# Basis is tenzorizable as ⊗ (1, Bi(t_i))
# and rs is tenzorizable as ⊗ (1 + r_i, r_i)
function tensor_rs_basis_with_eq(basis::Vector{T}, rs::Vector{T}, alpha::T = one(T)) where {T<:BinaryElem}
    @assert length(basis) == length(rs) "Vectors must be of the same length"

    ONE = one(T)
    acc = ONE
    @inbounds for i in eachindex(basis, rs)
        if i == 1
            acc *= alpha * ((ONE + rs[i]) + basis[i]*rs[i])
            continue
        end
        acc *= ((ONE + rs[i]) + basis[i]*rs[i])
    end
    return acc
end 

export tensor_rs_basis_with_eq