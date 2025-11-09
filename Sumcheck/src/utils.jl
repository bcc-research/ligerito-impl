# Basis is tenzorizable as ⊗ (1, Bi(t_i))
# and rs is tenzorizable as ⊗ (1 + r_i, r_i)
function tensor_rs_basis_with_eq(basis::Vector{T}, rs::Vector{T}, alpha::T = one(T)) where {T<:BinaryElem}
    @assert length(basis) == length(rs) "Vectors must be of the same length"

    ONE = one(T)
    acc = ONE
    n = length(rs)
    @inbounds for i in 1:n
        idx = n - i + 1  # reverse rs
        if i == 1
            acc *= alpha * ((ONE + rs[idx]) + basis[i]*rs[idx])
            continue
        end
        acc *= ((ONE + rs[idx]) + basis[i]*rs[idx])
    end
    return acc
end 

function partial_tensor_rs_basis_with_eq(basis::Vector{T}, rs::Vector{T}, alpha::T = one(T)) where {T<:BinaryElem}
    @assert length(basis) >= length(rs) "Basis vector must be at least as long as rs vector"
    ONE = one(T)
    acc = ONE
    n = length(rs)
    @inbounds for i in 1:n
        idx = n - i + 1  # reverse rs
        if i == 1
            acc *= alpha * ((ONE + rs[idx]) + basis[i]*rs[idx])
            continue
        end
        acc *= ((ONE + rs[idx]) + basis[i]*rs[idx])
    end
    return acc
end 

export tensor_rs_basis_with_eq, partial_tensor_rs_basis_with_eq