using BinaryFields, MultilinearPoly, SHA, Random

# define dummy FS 
mutable struct FS
    ctx::SHA256_CTX
end

function FS(seed::Int)
    ctx = SHA.SHA256_CTX()
    buf = Vector{UInt8}(undef, 8)
    unsafe_store!(Ptr{UInt64}(pointer(buf)), UInt64(seed))
    SHA.update!(ctx, buf)
    return FS(ctx)
end

function absorb!(fs::FS, s::QuadraticEvals)
    SHA.update!(fs.ctx, reinterpret(UInt8, [s.e0, s.e1, s.e2]))
end

function squeeze_rng(fs::FS)::MersenneTwister
    digest = SHA.digest!(deepcopy(fs.ctx))
    return MersenneTwister(reinterpret(UInt32, digest))
end

function get_field(fs::FS, ::Type{T}) where T <: BinaryElem
    rng = squeeze_rng(fs)
    return rand(rng, T)
end

# define utils functions 
function random_poly(::Type{T}, k::Int) where T <: BinaryElem
    evals = rand(T, 2^k)
    return MultiLinearPoly(evals)
end

function inner(f::MultiLinearPoly{T}, b::MultiLinearPoly{T}, rs::Vector{T}) where T <: BinaryElem
    n1 = f.n
    n2 = b.n

    eval_pts = rs[1:n1 - n2]
    fp = partial_eval(f, eval_pts)
    @assert fp.n == b.n

    return sum(fp.evals .* b.evals)
end

function inner_from_running(prover::SumcheckProverInstance{T}, b::MultiLinearPoly{T}) where T <: BinaryElem
    @assert prover.f.n == b.n
    return sum(prover.f.evals .* b.evals)
end

k = 12
glues = [2, 5, 9]
bs = [random_poly(BinaryElem16, k - gi) for gi in glues]
separation_challenges = rand(BinaryElem16, length(glues))

f = random_poly(BinaryElem16, k)
b1 = random_poly(BinaryElem16, k)
h = inner(f, b1, rs)

# store randomness for partial evaluations here
rs = Vector{BinaryElem16}()

let 
    # === PROVER ===
    fs_prover = FS(1234)
    prover, s1 = SumcheckProverInstance(f, b1, h)
    absorb!(fs_prover, s1)

    folds = 0
    gl_idx = 1
    hs = BinaryElem16[]

    for i in 1:(k - 1)
        ri = get_field(fs_prover, BinaryElem16)
        si = fold!(prover, ri)
        absorb!(fs_prover, si)
        folds += 1

        if gl_idx <= length(glues) && folds == glues[gl_idx]
            bi = bs[gl_idx]
            # alpha = separation_challenges[gl_idx]

            hi = inner_from_running(prover, bi)
            push!(hs, hi)

            gl_i = introduce_new!(prover, bi, hi)
            absorb!(fs_prover, gl_i)

            alpha = get_field(fs_prover, BinaryElem16)
            glue!(prover, alpha)

            gl_idx += 1
        end
    end

    # === VERIFIER ===
    fs_verifier = FS(1234)
    folds = 0
    gl_idx = 1
    verifier, g1 = SumcheckVerifierInstance(b1, h, copy(prover.transcript))
    absorb!(fs_verifier, g1)

    for i in 1:(k - 1)
        ri = get_field(fs_verifier, BinaryElem16)
        push!(rs, ri)
        gi = fold!(verifier, rs[i])
        absorb!(fs_verifier, gi)
        folds += 1

        if gl_idx <= length(glues) && folds == glues[gl_idx]
            bi = bs[gl_idx]
            hi = hs[gl_idx]
            # alpha = separation_challenges[gl_idx]

            gl_i = introduce_new!(verifier, bi, hi)
            absorb!(fs_verifier, gl_i)

            alpha = get_field(fs_verifier, BinaryElem16)
            glue!(verifier, alpha)

            gl_idx += 1
        end
    end
    # Final check

    # emulate oracle access to f
    ri = get_field(fs_verifier, BinaryElem16)
    push!(rs, ri)
    f_eval = partial_eval(f, rs).evals[1]

    # perform final verification 
    ok = verify(verifier, rs[end], f_eval)
    @show ok
end 

# k = 5 

# randomness for partial evaluations
# rs = rand(BinaryElem16, k)

# f_evals = rand(BinaryElem16, 2^k)
# f = MultiLinearPoly(f_evals)

# g_evals = rand(BinaryElem16, 2^k)
# g = MultiLinearPoly(g_evals)

# h1 = sum(f.evals .* g.evals)

# # start the prover 
# prover = SumcheckProverInstance(f, g, h1)

# fold!(prover, rs[1])
# fold!(prover, rs[2])
# fold!(prover, rs[3])
# fold!(prover, rs[4])

# copy transcript
# tr = copy(prover.transcript)

# # start the verifier
# verifier = SumcheckVerifierInstance(g, h1, tr)
# fold!(verifier, rs[1])
# fold!(verifier, rs[2])
# fold!(verifier, rs[3])
# fold!(verifier, rs[4])

# # simulate oracle access to f
# f_eval = partial_eval(f, rs).evals[1]

# # verify the sumcheck claim 
# verify(verifier, rs[5], f_eval)