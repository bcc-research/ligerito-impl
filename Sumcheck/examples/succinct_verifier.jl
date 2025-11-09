using BinaryFields, MultilinearPoly, BinaryReedSolomon

function random_poly(::Type{T}, k::Int) where T <: BinaryElem
    evals = rand(T, 2^k)
    return MultiLinearPoly(evals)
end

k = 20 
sks_at_vks = eval_sk_at_vks(2^k, BinaryElem32)
x = BinaryElem32(rand(1:2^k - 1))
sks_at_x = evaluate_sks_at_x(2^k, sks_at_vks, x)

f = random_poly(BinaryElem32, k)

