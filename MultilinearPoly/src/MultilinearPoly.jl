module MultilinearPoly

export MultiLinearPoly, partial_eval_at_0, partial_eval_at_1, partial_eval, eval_012, eval_013_product, quadratic_from_evals, eval_quadratic, sum

using BinaryFields

struct MultiLinearPoly{T<:BinaryElem}
    evals::Vector{T}
    n::Int
end

function MultiLinearPoly(evals::Vector{T}) where T
    n = Int(log2(length(evals)))
    if 2^n != length(evals)
        error("Length of evals must be a power of 2")
    end
    MultiLinearPoly{T}(evals, n)
end

function Base.sum(p::MultiLinearPoly)
    return sum(p.evals)
end

# [Lemma 4.3. from: Proofs, Arguments, and Zero-Knowledge, Justin Thaler] 
function partial_eval(p::MultiLinearPoly{T}, r::Vector{T}) where T<:BinaryElem
    partial_evals = copy(p.evals)
    len = length(partial_evals)
    ONE = one(eltype(partial_evals))

    for ri in r
        len = div(len, 2)

        for i in 1:len
            left = partial_evals[i]
            right = partial_evals[i + len]
            partial_evals[i] = (ONE + ri) * left + ri * right # original formula is -, in char 2 it's +
        end
    end
    return MultiLinearPoly(partial_evals[1:len])
end

function partial_eval_at_0(p::MultiLinearPoly{T}) where T<:BinaryElem
    half = div(length(p.evals), 2)
    return MultiLinearPoly(p.evals[1:half])
end

function partial_eval_at_1(p::MultiLinearPoly{T}) where T<:BinaryElem
    half = div(length(p.evals), 2)
    return MultiLinearPoly(p.evals[half+1:end])
end

# [Lemma 4.4. from: Proofs, Arguments, and Zero-Knowledge, Justin Thaler] 
function eval_012(f::MultiLinearPoly{T}) where T<:BinaryElem
    f0 = partial_eval_at_0(f)
    f1 = partial_eval_at_1(f)

    return sum(f0), sum(f1)
end

# [Lemma 4.4. from: Proofs, Arguments, and Zero-Knowledge, Justin Thaler] 
function eval_013_product(f::MultiLinearPoly{T}, g::MultiLinearPoly{T}) where T<:BinaryElem
    f0 = partial_eval_at_0(f)
    g0 = partial_eval_at_0(g)
    s0 = sum(fi * gi for (fi, gi) in zip(f0.evals, g0.evals))

    f1 = partial_eval_at_1(f)
    g1 = partial_eval_at_1(g)
    s1 = sum(fi * gi for (fi, gi) in zip(f1.evals, g1.evals))

    # in char 2 field we can't use F(2), so we simply use 3, then 3 + 1 = 2 (but in the extension)
    f2_vals = [T(3) * f1i + (T(1) + T(3)) * f0i for (f1i, f0i) in zip(f1.evals, f0.evals)]
    g2_vals = [T(3) * g1i + (T(1) + T(3)) * g0i for (g1i, g0i) in zip(g1.evals, g0.evals)]
    s2 = sum(f2i * g2i for (f2i, g2i) in zip(f2_vals, g2_vals))

    return s0, s1, s2
end

end