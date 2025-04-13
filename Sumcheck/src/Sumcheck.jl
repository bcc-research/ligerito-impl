module Sumcheck

export sumcheck_prover, sumcheck_verifier, double_sumcheck_prover, double_sumcheck_verifier
export SumcheckProverInstance, SumcheckVerifierInstance, fold!, glue!, introduce_new!, verify

using BinaryFields, MultilinearPoly

export QuadraticPoly

# aX^2 + bX + c
struct QuadraticPoly{T<:BinaryElem}
    a::T
    b::T
    c::T
end

# interpolates quadratic polynomial from evaluations at 0, 1 and x where we change each - to + in char 2 field
function quadratic_from_evals(at0::T, at1::T, atx::T, x = T(3)) where T <: BinaryElem
    numerator = atx + at0 + x * (at1 + at0)
    denominator = x*x + x
    a = numerator * inv(denominator)
    b = at1 + at0 + a
    return QuadraticPoly(a, b, at0)
end

function eval_quadratic(poly::QuadraticPoly{T}, r::T) where T <: BinaryElem
    return poly.a * r * r + poly.b * r + poly.c
end

function fold_quadratic(p1::QuadraticPoly{T}, p2::QuadraticPoly{T}, alpha::T) where T <: BinaryElem
    return QuadraticPoly(
        p1.a + alpha * p2.a,
        p1.b + alpha * p2.b,
        p1.c + alpha * p2.c,
    )
end

include("./ligeritho.jl")

#TODO: remove claimed_sum after testing
function sumcheck_prover(f::MultiLinearPoly{T}, claimed_sum::T) where T<:BinaryElem
    transcript = Vector{NTuple{2, T}}()
    ris = T[]
    current_poly = f
    current_sum = claimed_sum

    n = Int(log2(length(f.evals)))
    for _ in 1:n
        s0, s1 = eval_012(current_poly)
        @assert s0 + s1 == current_sum

        push!(transcript, (s0, s1))

        r_i = rand(T)
        push!(ris, r_i)

        current_poly = partial_eval(current_poly, [r_i])
        current_sum = sum(current_poly)
    end

    return transcript, ris
end

function sumcheck_verifier(transcript::Vector{Tuple{T,T}}, ris::Vector{T},
    claimed_sum::T, f::MultiLinearPoly{T}) where T<:BinaryElem

    h = claimed_sum
    for (i, (g0, g1)) in enumerate(transcript)
        @assert g0 + g1 == h
        h = g0 * (one(T) + ris[i]) + g1 * ris[i] # g = g0 * (1 - x) + g1 * x, which is + in char 2 field
    end

    f_eval = sum(partial_eval(f, ris))
    return f_eval == h
end


function double_sumcheck_prover(fp::MultiLinearPoly{T}, gp::MultiLinearPoly{T}) where T<:BinaryElem
    transcript = Vector{NTuple{3, T}}()
    ris = T[]
    f, g = fp, gp

    n = Int(log2(length(fp.evals)))
    for _ in 1:n
        s0, s1, s2 = eval_013_product(f, g)
        push!(transcript, (s0, s1, s2))

        r_i = rand(T)
        push!(ris, r_i)

        f = partial_eval(f, [r_i])
        g = partial_eval(g, [r_i])
    end

    return transcript, ris
end

function double_sumcheck_verifier(f::MultiLinearPoly{T}, g::MultiLinearPoly{T}, claimed_sum::T, transcript::Vector{NTuple{3, T}}, ris::Vector{T}) where T<:BinaryElem
    H = claimed_sum
    for (i, (gi0, gi1, gi2)) in enumerate(transcript)
        @assert gi0 + gi1 == H
        gi = quadratic_from_evals(gi0, gi1, gi2)
        H = eval_quadratic(gi, ris[i])
    end

    f_eval = sum(partial_eval(f, ris))
    g_eval = sum(partial_eval(g, ris))

    return f_eval * g_eval == H
end

end # module Sumcheck
