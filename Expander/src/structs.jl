# struct ChunkedConstraint{T <:BinaryElem}
#     chunks::Vector{Vector{Tuple{Int, T}}} 
# end 

# @inline function determine_bucket(i::Int, bounds::AbstractVector{Int})
#     j = searchsortedfirst(bounds, i)
#     return j > length(bounds) ? length(bounds) : j
# end

# function ChunkedConstraint(msg_len::Int, num_of_chunks::Int, d::Int, bounds::AbstractVector{Int}, ::Type{T}) where T <: BinaryElem
#     chunks = [Vector{Tuple{Int, T}}() for _ in 1:num_of_chunks]
#     constraint_indices = sort(sample(1:msg_len, d; replace=false))
#     constraint_values = rand(T, d)

#     for (i, val) in zip(constraint_indices, constraint_values)
#         bucket = determine_bucket(i, bounds)
#         push!(chunks[bucket], (i, val))
#     end
#     return ChunkedConstraint{T}(chunks)
# end

# struct ExpanderMatrix{T <:BinaryElem}
#     msg_len::Int
#     rows::Vector{ChunkedConstraint{T}} 
# end

# function ExpanderMatrix(num_of_constraints::Int, msg_len::Int, num_of_chunks::Int, d::Int, chunk_size::Int, ::Type{T}) where T <: BinaryElem
#     bounds = [i * chunk_size for i in 1:num_of_chunks-1]
#     push!(bounds, msg_len)
#     rows = Vector{ChunkedConstraint{T}}(undef, num_of_constraints)

#     for i in 1:num_of_constraints
#         rows[i] = ChunkedConstraint(msg_len, num_of_chunks, d, bounds, T)
#     end
#     return ExpanderMatrix{T}(msg_len, rows)
# end

# struct SimpleExpander{T<:BinaryElem}
#     msg_len::Int
#     constraints::Vector{Vector{Tuple{Int, T}}}
# end

# function SimpleExpander(num_of_constraints::Int, msg_len::Int, d::Int, ::Type{T}) where T <: BinaryElem
#     constraints = Vector{Vector{Tuple{Int, T}}}(undef, num_of_constraints)

#     for i in 1:num_of_constraints
#         constraint_indices = sort(sample(1:msg_len, d; replace=false))
#         constraint_values = rand(T, d)
#         constraints[i] = [(idx, val) for (idx, val) in zip(constraint_indices, constraint_values)]
#     end
#     return SimpleExpander{T}(msg_len, constraints)
# end

# function encode_simple!(c::Vector{T}, expander::SimpleExpander{T}) where T <: BinaryElem
#     for (i, constraint) in enumerate(expander.constraints)
#         acc = zero(T)
#         for (idx, val) in constraint
#             acc += c[idx] * val
#         end
#         c[expander.msg_len + i] = acc
#     end
# end

# using Base.Threads: @threads
# function encode_simple_threads!(c::Vector{T}, expander::SimpleExpander{T}) where {T<:BinaryElem}
#     m = length(expander.constraints)
#     msg_len = expander.msg_len
#     @threads for i in 1:m
#         @inbounds for (idx, val) in expander.constraints[i]
#             c[msg_len + i] += c[idx] * val
#         end
#     end
# end

# function encode_msgs_with_threads!(msgs::Vector{Vector{T}}, expander::SimpleExpander{T}) where T <: BinaryElem
#     @threads for msg in msgs
#         encode_simple!(msg, expander)
#     end
# end

# function encode_msgs_manual_chunks!(msgs::Vector{Vector{T}}, expander::SimpleExpander{T}) where T <: BinaryElem
#     nt = Threads.nthreads()
#     n = length(msgs)
#     chunk_size = div(n, nt)
    
#     Threads.@sync for t in 1:nt
#         Threads.@spawn begin
#             start_idx = (t - 1) * chunk_size + 1
#             stop_idx = t == nt ? n : t * chunk_size
            
#             @inbounds for i in start_idx:stop_idx
#                 encode_simple!(msgs[i], expander)
#             end
#         end
#     end
# end

# function encode_msgs!(msgs::Vector{Vector{T}}, expander::SimpleExpander{T}) where T <: BinaryElem
#     for msg in msgs
#         encode_simple_threads!(msg, expander)
#     end
# end

# function encode_msgs_simple!(msgs::Vector{Vector{T}}, expander::SimpleExpander{T}) where T <: BinaryElem
#     for msg in msgs
#         encode_simple_threads!(msg, expander)
#     end
# end

# function encode_msgs_expander!(msgs::Vector{Vector{T}}, expander::ExpanderMatrix{T}) where T <: BinaryElem
#     @threads for msg in msgs
#         encode!(msg, expander)
#     end
# end

# export ChunkedConstraint, ExpanderMatrix, SimpleExpander, encode_simple!, encode_simple_threads!, encode_msgs!, encode_msgs_simple!, encode_msgs_expander!
# export encode_msgs_with_threads!, encode_msgs_manual_chunks!