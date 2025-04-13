using BinaryFields
using MultilinearPoly

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
rs = rand(BinaryElem16, k)
glues = [2, 5, 9]
bs = [random_poly(BinaryElem16, k - gi) for gi in glues]
separation_challenges = rand(BinaryElem16, length(glues))

f = random_poly(BinaryElem16, k)
b1 = random_poly(BinaryElem16, k)
h = inner(f, b1, rs)

let 
    # === PROVER ===
    prover = SumcheckProverInstance(f, b1, h)
    folds = 0
    gl_idx = 1
    hs = BinaryElem16[]

    for i in 1:(k - 1)
        fold!(prover, rs[i])
        folds += 1

        if gl_idx <= length(glues) && folds == glues[gl_idx]
            bi = bs[gl_idx]
            alpha = separation_challenges[gl_idx]

            hi = inner_from_running(prover, bi)
            push!(hs, hi)

            introduce_new!(prover, bi, hi)
            glue!(prover, alpha)

            gl_idx += 1
        end
    end

    # === VERIFIER ===
    folds = 0
    gl_idx = 1
    verifier = SumcheckVerifierInstance(b1, h, copy(prover.transcript))

    for i in 1:(k - 1)
        fold!(verifier, rs[i])
        folds += 1

        if gl_idx <= length(glues) && folds == glues[gl_idx]
            bi = bs[gl_idx]
            hi = hs[gl_idx]
            alpha = separation_challenges[gl_idx]

            introduce_new!(verifier, bi, hi)
            glue!(verifier, alpha)

            gl_idx += 1
        end
    end

    # Final check

    # emulate oracle access to f
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