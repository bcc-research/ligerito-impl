using BinaryFields
struct SumcheckProverInstance{T<:BinaryElem}
    f::MultiLinearPoly{T}
    basis_polys::Vector{MultiLinearPoly{T}}
    sum::T
    transcript::Vector{NTuple{3, T}}
    to_be_glued::Union{Nothing, MultiLinearPoly{T}}
end

function SumcheckProverInstance(f::MultiLinearPoly{T}, b1::MultiLinearPoly{T}, h1::T) where T <: BinaryElem
    basis_polys = [b1]
    transcript = NTuple{3, T}[]
    instance = SumcheckProverInstance(f, basis_polys, h1, transcript, nothing)
    _start!(instance)
    return instance
end

function _start!(inst::SumcheckProverInstance)
    s0, s1, s2 = eval_013_product(inst.f, inst.basis_polys[1])
    push!(inst.transcript, (s0, s1, s2))
    @assert s0 + s1 == inst.sum
end

function introduce_new!(inst::SumcheckProverInstance, bi::MultiLinearPoly{T}, h::T) where T <: BinaryElem
    s0, s1, s2 = eval_012_product(inst.f, bi)
    push!(inst.transcript, (s0, s1, s2))
    push!(inst.basis_polys, bi)
    @assert s0 + s1 == h
    inst.to_be_glued = bi
end

function glue!(inst::SumcheckProverInstance{T}, alpha::T) where T
    @assert inst.to_be_glued !== nothing "No polynomial to glue!"
    bi = inst.to_be_glued
    g_evals = [gi * alpha for gi in bi.evals]
    g = MultiLinearPoly(g_evals)
    push!(inst.basis_polys, g)
    inst.to_be_glued = nothing
end


function evals_at_x(evals_at_0::Vector{T}, evals_at_1::Vector{T}, x::T) where T <: BinaryElem
    return x .* evals_at_1 .+ (one(T) + x) .* evals_at_0
end

function eval_01x_product(inst::SumcheckProverInstance{T}, x::T = T(3)) where T <: BinaryElem
    f0 = partial_eval_at_0(inst.f).evals
    f1 = partial_eval_at_1(inst.f).evals
    f2 = evals_at_x(f0, f1, x)

    b0s_sum = partial_eval_at_0(inst.basis_polys[1]).evals
    b1s_sum = partial_eval_at_1(inst.basis_polys[1]).evals
    b2s_sum = evals_at_x(b0s_sum, b1s_sum, x)

    for i in 2:length(inst.basis_polys)
        b0i = partial_eval_at_0(inst.basis_polys[i]).evals
        b1i = partial_eval_at_1(inst.basis_polys[i]).evals
        b2i = evals_at_x(b0i, b1i, x)

        b0s_sum .= b0s_sum .+ b0i
        b1s_sum .= b1s_sum .+ b1i
        b2s_sum .= b2s_sum .+ b2i
    end

    s0 = sum(f0 .* b0s_sum)
    s1 = sum(f1 .* b1s_sum)
    s2 = sum(f2 .* b2s_sum)

    return s0, s1, s2
end

function fold!(inst::SumcheckProverInstance{T}, r::T) where T <: BinaryElem
    inst.f = partial_eval(inst.f, [r])
    for i in 1:length(inst.basis_polys)
        inst.basis_polys[i] = partial_eval(inst.basis_polys[i], [r])
    end

    s0, s1, s2 = eval_01x_product(inst)
    push!(inst.transcript, (s0, s1, s2))
end

struct SumcheckVerifierInstance{T}
    basis_polys::Vector{MultiLinearPoly{T}}
    separation_challenges::Vector{T}
    sum::T
    transcript::Vector{NTuple{3, T}}
    ris::Vector{T}
    tr_reader::Int

    running_poly::Union{Nothing, QuadraticPoly{T}}
    to_glue::Union{Nothing, QuadraticPoly{T}}
end


function read_tr!(verifier::SumcheckVerifierInstance{T}) where T <: BinaryElem
    @assert verifier.index ≤ length(verifier.transcript) "Transcript exhausted"
    g0, g1, g2 = verifier.transcript[verifier.index]
    verifier.index += 1
    return g0, g1, g2
end

function SumcheckVerifierInstance(b1::MultiLinearPoly{T}, h1::T, transcript::Vector{NTuple{3, T}}) where T <: BinaryElem
    verifier = SumcheckVerifierInstance([b1], [one(T)], h1, copy(transcript), T[], 1, nothing, nothing)
    start!(verifier)
    return verifier
end 

function start!(verifier::SumcheckVerifierInstance{T}) where T
    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == verifier.sum

    verifier.running_poly = quadratic_from_evals(g0, g1, g2)
end


function fold!(verifier::SumcheckVerifierInstance{T}, r::T) where T <: BinaryElem
    push!(verifier.ris, r)
    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == verifier.sum

    g = quadratic_from_evals(g0, g1, g2)
    verifier.sum = eval_quadratic(g, r)
end

function introduce_new!(verifier::SumcheckVerifierInstance{T}, bi::MultiLinearPoly{T}, h::T) where T
    g0, g1, g2 = read_tr!(verifier)
    @assert g0 + g1 == h

    push!(verifier.basis_polys, bi)
    verifier.to_glue = quadratic_from_evals(g0, g1, g2)
end

function glue!(verifier::SumcheckVerifierInstance{T}, alpha::T) where T
    push!(verifier.separation_challenges, alpha)
    verifier.running_poly = fold_quadratic(verifier.running_poly, verifier.to_glue, alpha)
end

function evaluate_basis_polys(verifier::SumcheckVerifierInstance{T}, r::T) where T
    push!(verifier.ris, r)
    b_eval = partial_eval(verifier.basis_polys[1], verifier.ris).evals[1]

    for i in 2:length(verifier.basis_polys)
        n = verifier.basis_polys[i].n
        eval_pts = verifier.ris[end - n + 1:end]
        bi_eval = partial_eval(verifier.basis_polys[i], eval_pts).evals[1]
        b_eval += verifier.separation_challenges[i] * bi_eval
    end

    return b_eval
end

function verify(verifier::SumcheckVerifierInstance{T}, r::T, f_eval::T) where T
    verifier.sum = eval_quadratic(verifier.running_poly, r)
    basis_evals = evaluate_basis_polys(verifier, r)
    return f_eval * basis_evals == verifier.sum
end
