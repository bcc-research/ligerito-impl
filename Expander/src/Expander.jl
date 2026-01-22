module Expander

using BinaryFields, StatsBase, Random
using Base.Threads: @threads

struct ExpanderMatrix{T<:BinaryElem}
    msg_len::Int
    constraints::Vector{Vector{Tuple{Int, T}}}
end

function ExpanderMatrix(num_of_constraints::Int, msg_len::Int, d::Int, ::Type{T}) where T <: BinaryElem
    constraints = Vector{Vector{Tuple{Int, T}}}(undef, num_of_constraints)

    for i in 1:num_of_constraints
        constraint_indices = sort(sample(1:msg_len, d; replace=false))
        constraint_values = rand(T, d)
        constraints[i] = [(idx, val) for (idx, val) in zip(constraint_indices, constraint_values)]
    end
    return ExpanderMatrix{T}(msg_len, constraints)
end

function encode!(c, expander::ExpanderMatrix{T}) where T <: BinaryElem
    for (i, constraint) in enumerate(expander.constraints)
        acc = zero(T)
        for (idx, val) in constraint
            acc += c[idx] * val
        end
        c[expander.msg_len + i] = acc
    end
end

export ExpanderMatrix, encode!


# const MIB = 1_048_576 # bytes 
# const L2_CACHE_SIZE_BYTES_PER_CORE = 3 * MIB
# const ELEMENT_SIZE_BYTES = 4 # BinaryElem32
# const CHUNK_SIZE = div(L2_CACHE_SIZE_BYTES_PER_CORE, ELEMENT_SIZE_BYTES)
# const CACHELINE = 128
# const STEP = CACHELINE ÷ ELEMENT_SIZE_BYTES

# include("structs.jl")

# @inline function prefetch(ptr::Ptr{T}) where T
#     Base.llvmcall
#     (
#         "declare void @llvm.prefetch(ptr <address>, i32 <rw>, i32 <locality>, i32 <cache type>)", 
#         Cvoid, 
#         (Ptr{T}, Cint, Cint, Cint), 
#         ptr, 0, 3, 1
#     )
# end

# function prefetch_span!(A::AbstractVector, start::Int, stop::Int)
#     baseptr = pointer(A, start)
#     n = stop - start + 1
#     @inbounds for off in 0:STEP:(n-1)
#         prefetch(Ptr{Cvoid}(baseptr + off))
#     end
#     return nothing
# end

# function to_bench_random(n::Int, d::Int)
#     indices = sample(1:n, d; replace=false)
#     values = rand(BinaryElem32, d)
#     x = rand(BinaryElem32, n)
#     return x, indices, values
# end 

# function to_bench_chunk(n::Int, d::Int, chunk_idx::Int)
#     start = (chunk_idx - 1) * CHUNK_SIZE + 1
#     stop = min(chunk_idx * CHUNK_SIZE, n)
#     indices = sample(start:stop, d; replace=false)
#     values = rand(BinaryElem32, d)
#     x = rand(BinaryElem32, n)
#     return x, indices, values
# end 

# function random_access(x::Vector{T}, indices::Vector{Int}, values::Vector{T}) where T <: BinaryElem
#     acc = zero(T) 
#     for i in eachindex(indices)
#         acc += x[indices[i]] * values[i]
#     end
#     return acc
# end

# function access_with_fetch(x::Vector{T}, indices::Vector{Int}, values::Vector{T}, chunk_idx::Int) where T <: BinaryElem
#     start = (chunk_idx - 1) * CHUNK_SIZE + 1
#     stop = min(chunk_idx * CHUNK_SIZE, length(x))
#     prefetch_span!(x, start, stop)
#     acc = zero(T) 
#     for i in eachindex(indices)
#         acc += x[indices[i]] * values[i]
#     end
#     return acc
# end

# # assumes that c is already allocated
# function encode!(c::Vector{T}, expander::ExpanderMatrix{T}) where T <: BinaryElem
#     num_of_chunks = length(expander.rows[1].chunks)
#     for chunk_idx in 1:num_of_chunks
#         start = (chunk_idx - 1) * CHUNK_SIZE + 1
#         stop = min(chunk_idx * CHUNK_SIZE, expander.msg_len)
#         prefetch_span!(c, start, stop)
#         encode_chunk!(c, expander, chunk_idx)
#     end
# end

# function encode_chunk!(c::Vector{T}, expander::ExpanderMatrix{T}, chunk_idx::Int) where T <: BinaryElem
#     msg_len = expander.msg_len
#     parity_len = length(expander.rows)

#     parity = @view c[msg_len+1 : msg_len+parity_len]
#     for i in 1:parity_len
#         constraint = expander.rows[i].chunks[chunk_idx]
#         for (idx, val) in constraint
#             parity[i] += c[idx] * val
#         end
#     end
# end

# function inp(x::Vector{T}, y::Vector{T}) where T <: BinaryElem
#     acc = zero(T) 
#     for i in eachindex(x)
#         acc += x[i] * y[i]
#     end
#     return acc
# end


# export L2_CACHE_SIZE_BYTES_PER_CORE, ELEMENT_SIZE_BYTES, CHUNK_SIZE
# export prefetch, encode!

# export random_access, to_bench_random, to_bench_chunk, access_with_fetch, inp

end # module Expander