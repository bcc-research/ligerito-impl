using BinaryReedSolomon

mutable struct FoldData{T <: BinaryElem}
    n::Int
    sks_evaluations::Vector{Vector{T}}
    alpha::T
    beta::T
    rs::Vector{T}
end 

function FoldData(n::Int, queried_points::Vector{T}, alpha::U, sks_vks::AbstractVector{T}) where {U <: BinaryElem, T <: BinaryElem}
    sks_at_x = [evaluate_sks_at_x(2^n, sks_vks, x) for x in queried_points]
    FoldData{U}(n, sks_at_x, alpha, one(U), U[])
end

function register_separation_challenge!(fd::FoldData{T}, beta::T) where T<:BinaryElem
    fd.beta = beta
end

function add_r!(fd::FoldData{T}, r::T) where T<:BinaryElem
    push!(fd.rs, r)
end

function evaluate_poly(fd::FoldData{T}) where T<:BinaryElem 
    powers = Vector{T}(undef, length(fd.sks_evaluations))
    powers[1] = one(T)
    for i in 2:length(fd.sks_evaluations)
        powers[i] = powers[i-1] * fd.alpha
    end
    if fd.beta != one(T)
        powers .= fd.beta * powers
    end
    sum(tensor_rs_basis_with_eq(fd.sks_evaluations[i], fd.rs, powers[i]) for i in eachindex(fd.sks_evaluations))
end

mutable struct SumcheckVerifier{T}
    folds::Vector{FoldData{T}}
    sum::T
    transcript::Vector{NTuple{3, T}}
    ris::Vector{T}
    tr_reader::Int

    running_poly::Union{Nothing, QuadraticPoly{T}}
    to_glue::Union{Nothing, QuadraticPoly{T}}
end

function read_tr!(verifier::SumcheckVerifier{T}) where T <: BinaryElem
    @assert verifier.tr_reader ≤ length(verifier.transcript) "Transcript exhausted"
    g0, g1, g2 = verifier.transcript[verifier.tr_reader]
    verifier.tr_reader += 1
    return g0, g1, g2
end

function SumcheckVerifier(transcript::Vector{NTuple{3, T}}) where T <: BinaryElem
    SumcheckVerifier{T}([], zero(T), copy(transcript), T[], 1, nothing, nothing)
end 

function start!(verifier::SumcheckVerifier{T}, h::T) where T <: BinaryElem
    verifier.sum = h
    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == verifier.sum

    verifier.running_poly = quadratic_from_evals(g0, g1, g2)
    return QuadraticEvals(g0, g1, g2)
end

function new_fold!(verifier::SumcheckVerifier{U}, n::Int, queried_points::Vector{T}, sks_vks::AbstractVector{T}, alpha::U) where {U <: BinaryElem, T <: BinaryElem}
    fold = FoldData(n, queried_points, alpha, sks_vks)
    push!(verifier.folds, fold)
end

function add_r!(v::SumcheckVerifier{T}, r::T) where T <: BinaryElem
    push!(v.ris, r)
    for fold in v.folds
        add_r!(fold, r)
    end    
end


function fold!(verifier::SumcheckVerifier{T}, r::T) where T <: BinaryElem
    push!(verifier.ris, r)
    for fold in verifier.folds
        add_r!(fold, r)
    end
    verifier.sum = eval_quadratic(verifier.running_poly, r)

    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == verifier.sum

    verifier.running_poly = quadratic_from_evals(g0, g1, g2)
    return QuadraticEvals(g0, g1, g2)
end

function verifier_introduce_new!(verifier::SumcheckVerifier{T}, h::T) where T <: BinaryElem
    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == h

    verifier.to_glue = quadratic_from_evals(g0, g1, g2)
    return QuadraticEvals(g0, g1, g2)
end

function glue!(verifier::SumcheckVerifier{T}, beta::T) where T
    register_separation_challenge!(verifier.folds[end], beta)
    verifier.running_poly = fold_quadratic(verifier.running_poly, verifier.to_glue, beta)
end

function evaluate_basis_polys(verifier::SumcheckVerifier{T}, r::T) where T <: BinaryElem
    for fold in verifier.folds
        add_r!(fold, r)
    end
    return prod(evaluate_poly(fold) for fold in verifier.folds)
end

# function verify(verifier::SumcheckVerifier{T}, r::T, f_eval::T) where T <: BinaryElem
#     verifier.sum = eval_quadratic(verifier.running_poly, r)
#     basis_evals = evaluate_basis_polys(verifier, r)
#     return f_eval * basis_evals == verifier.sum
# end

function verify(verifier::SumcheckVerifier{T}, r::T, f_partial_eval::Vector{T}) where T <: BinaryElem
    k = Int(log2(length(f_partial_eval)))

    f_poly = MultiLinearPoly(f_partial_eval)
    final_rs = vcat([r], rand(T, k - 1))
    for ri in final_rs
        for fold in verifier.folds
            add_r!(fold, ri)
        end
        f_poly = partial_eval(f_poly, [ri])
    end
    verifier.sum = eval_quadratic(verifier.running_poly, r)
    basis_evals = evaluate_basis_polys(verifier, r)
    return f_poly.evals[1] * basis_evals == verifier.sum
end


export FoldData, add_r!, evaluate_poly, SumcheckVerifier, fold!, introduce_new!, glue!, verify, new_fold!, start!