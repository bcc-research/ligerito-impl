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
rs = rand(BinaryElem16, k - 4) # for example stop after 8 folds 
glues = [2, 5]
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

    # prover folds rs - 1 times in total 
    for i in 1:(length(rs) - 1)
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

    # verifier folds rs - 1 times before the final check
    for i in 1:(length(rs) - 1)
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
    f_partial_eval = partial_eval(f, rs).evals

    # perform final verification 
    ok = verify_partial(verifier, rs[end], f_partial_eval)
    @show ok
end 