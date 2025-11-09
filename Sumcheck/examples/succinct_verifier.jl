# using BinaryFields, MultilinearPoly, BinaryReedSolomon

# function random_poly(::Type{T}, k::Int) where T <: BinaryElem
#     evals = rand(T, 2^k)
#     return MultiLinearPoly(evals)
# end

# #= 
# We can start with a random polynomial that that we can partially evaluate at a few random points 
# =#

