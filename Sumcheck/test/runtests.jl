using Test
using Sumcheck, BinaryReedSolomon, BinaryFields, MultilinearPoly

function expand_evaluation_point(rs::Vector{T}) where T <: BinaryElem
    one_elem = one(T)
    current_layer = [one_elem + rs[1], rs[1]]
    len = 2
    for i in 2:length(rs)
        next_layer_size = 2 * len
        next_layer = Vector{T}(undef, next_layer_size)

        ri_p_one = one_elem + rs[i]
        for j in 1:len
            next_layer[2*j - 1] = current_layer[j] * ri_p_one
            next_layer[2*j] = current_layer[j] * rs[i]
        end

        current_layer = next_layer
        len *= 2
    end

    return current_layer
end 

@testset "Kronecker basic test" begin
    T = BinaryElem32 
    n = 20

    r = rand(T, n)
    sks_vks = eval_sk_at_vks(2^n, T)

    t = T(rand(1:2^n - 1))
    sks_at_x = evaluate_sks_at_x(2^n, sks_vks, t)

    basis = expand_basis(sks_at_x)
    eq_r = expand_evaluation_point(r)

    inp = sum(basis .* eq_r)
    tensor_ip = tensor_rs_basis_with_eq(sks_at_x, r)
    @assert inp == tensor_ip "Kronecker product inner product does not match direct computation"
end

@testset "Separation challenges test" begin
    T = BinaryElem32 
    n = 20

    r = rand(T, n)
    sks_vks = eval_sk_at_vks(2^n, T)

    t1 = T(rand(1:2^n - 1))
    t2 = T(rand(1:2^n - 1))
    sks_at_t1 = evaluate_sks_at_x(2^n, sks_vks, t1)
    sks_at_t2 = evaluate_sks_at_x(2^n, sks_vks, t2)

    alpha1 = rand(T)
    alpha2 = rand(T)

    # so we want do check inner product of the form 
    # <(alpha1 * basis_1 + alpha2 * basis_2), eq_r>
    basis_1 = expand_basis(sks_at_t1)
    basis_2 = expand_basis(sks_at_t2)
    eq_r = expand_evaluation_point(r)

    combined_basis = alpha1 .* basis_1 .+ alpha2 .* basis_2
    inp = sum(combined_basis .* eq_r)

    tensor_ip1 = tensor_rs_basis_with_eq(sks_at_t1, r, alpha1)
    tensor_ip2 = tensor_rs_basis_with_eq(sks_at_t2, r, alpha2)

    @assert inp == tensor_ip1 + tensor_ip2 "Kronecker product inner product with separation challenges does not match direct computation"
end

@testset "Partial basis test" begin
    T = BinaryElem32 
    n = 20

    alpha = rand(T)
    r = rand(T, n)
    sks_vks = eval_sk_at_vks(2^n, T)

    t = T(rand(1:2^n - 1))
    sks_at_x = evaluate_sks_at_x(2^n, sks_vks, t)

    full_eval = tensor_rs_basis_with_eq(sks_at_x, r, alpha)

    # now let's partial evaluate on first 10 variables
    k = 10
    r_partial = r[1:k]
    basis_at_partial_r = partial_tensor_rs_basis_with_eq(sks_at_x[k+1:end], r_partial, alpha)

    rest_of_basis = expand_basis(sks_at_x[1:k])
    rest_of_basis .= basis_at_partial_r .* rest_of_basis

    poly = MultiLinearPoly(rest_of_basis)
    final_eval = partial_eval(poly, r[k+1:end])
    @assert length(final_eval.evals) == 1 "Final evaluation should be a single element"
    @assert full_eval == final_eval.evals[1] "Partial basis evaluation does not match full evaluation"
end