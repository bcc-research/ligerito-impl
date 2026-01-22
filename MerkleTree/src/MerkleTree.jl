module MerkleTree

using SHA, BCCSHA
import Base: sizeof

export build_merkle_tree, get_root, get_depth, build_merkle_tree_fast, build_merkle_tree_odd
export CompleteMerkleTree, MerkleRoot, BatchedMerkleProof
export prove, verify


# for now simply assume that isbitstype(T) = true
hash_leaf(leaf) = SHA.sha256(reinterpret(UInt8, leaf))
hash_siblings(left_hash::Vector{UInt8}, right_hash::Vector{UInt8}) = SHA.sha256(vcat(left_hash, right_hash))

function is_power_of_two(n::Int)
    return n > 0 && (n & (n - 1)) == 0
end

struct CompleteMerkleTree
    tree::Vector{Vector{UInt8}}  # Each layer is flat 32*n bytes
end

function reset_ctx!(ctx::BCCSHA.SHA256_CTX)
    copyto!(ctx.state, BCCSHA.SHA2_256_initial_hash_value)
    ctx.bytecount = 0
    fill!(ctx.buffer, 0x00)
    ctx.used = false
end

function hash_array(x::AbstractArray)
    n = length(x)
    xh = Vector{UInt8}(undef, 32 * n)

    xh_ptr = pointer(xh)

    nt = Threads.nthreads() 
    chunk_size = ceil(Int, n / nt)
    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            ctx = BCCSHA.SHA256_CTX()
            start_idx = (t - 1) * chunk_size + 1
            end_idx = min(t * chunk_size, n)

            @inbounds for i in start_idx:end_idx
                BCCSHA.update!(ctx, reinterpret(UInt8, x[i]))
                BCCSHA.digest_inplace!(xh_ptr + (i-1)*32, ctx)
                reset_ctx!(ctx)
            end
        end
    end

    return xh
end

function build_merkle_tree(leaves::AbstractArray)
    @assert is_power_of_two(length(leaves))
    if isempty(leaves)
        return []
    end

    n = length(leaves)
    lh = hash_array(leaves)
    layers = [lh]

    current_layer = lh
    current_n = n

    while current_n > 1
        next_n = (current_n + 1) ÷ 2  # ceil div by 2
        next_layer = Vector{UInt8}(undef, 32 * next_n)
        next_ptr = pointer(next_layer)
        curr_ptr = pointer(current_layer)

        nt = Threads.nthreads()
        chunk_size = ceil(Int, next_n / nt)

        Threads.@sync for t in 1:nt
            Threads.@spawn begin
                ctx = BCCSHA.SHA256_CTX()
                buf64 = Vector{UInt8}(undef, 64)
                buf_prt = pointer(buf64)

                start_idx = (t-1)*chunk_size + 1
                end_idx = min(t*chunk_size, next_n)

                @inbounds for i in start_idx:end_idx
                    src = curr_ptr + (2*(i-1))*32
                    unsafe_copyto!(buf_prt, src, 64)

                    BCCSHA.update!(ctx, buf64)
                    BCCSHA.digest_inplace!(next_ptr + (i-1)*32, ctx)
                    reset_ctx!(ctx)
                end
            end
        end

        push!(layers, next_layer)

        current_layer = next_layer
        current_n = next_n
    end

    return CompleteMerkleTree(layers)
end

function build_merkle_tree_odd(leaves::AbstractArray)
    if isempty(leaves)
        return []
    end

    n = length(leaves)
    lh = hash_array(leaves)
    layers = [lh]

    current_layer = lh
    current_n = n

    # Dummy hash - zeros
    dummy_hash = zeros(UInt8, 32)

    while current_n > 1
        next_n = (current_n + 1) ÷ 2  # ceil div by 2
        next_layer = Vector{UInt8}(undef, 32 * next_n)
        next_ptr = pointer(next_layer)
        curr_ptr = pointer(current_layer)

        # Check if we need to pad with dummy
        is_odd = (current_n % 2 == 1)

        nt = Threads.nthreads()
        chunk_size = ceil(Int, next_n / nt)

        Threads.@sync for t in 1:nt
            Threads.@spawn begin
                ctx = BCCSHA.SHA256_CTX()
                buf64 = Vector{UInt8}(undef, 64)
                buf_prt = pointer(buf64)

                start_idx = (t-1)*chunk_size + 1
                end_idx = min(t*chunk_size, next_n)

                @inbounds for i in start_idx:end_idx
                    # Check if this is the last element and layer is odd
                    if is_odd && i == next_n
                        # Hash the last node with dummy
                        src = curr_ptr + (2*(i-1))*32
                        unsafe_copyto!(buf_prt, src, 32)
                        unsafe_copyto!(buf_prt + 32, pointer(dummy_hash), 32)
                    else
                        src = curr_ptr + (2*(i-1))*32
                        unsafe_copyto!(buf_prt, src, 64)
                    end

                    BCCSHA.update!(ctx, buf64)
                    BCCSHA.digest_inplace!(next_ptr + (i-1)*32, ctx)
                    reset_ctx!(ctx)
                end
            end
        end

        push!(layers, next_layer)

        current_layer = next_layer
        current_n = next_n
    end

    return CompleteMerkleTree(layers)
end

struct MerkleRoot
    root::Vector{UInt8}
end

sizeof(x::MerkleRoot) = sizeof(x.root)
get_root(c::CompleteMerkleTree) = MerkleRoot(isempty(c.tree) ? nothing : c.tree[end])
get_depth(c::CompleteMerkleTree) = length(c.tree) - 1

struct FlatLayerView
    data::Vector{UInt8}
end

function Base.getindex(flv::FlatLayerView, idx::Int)
    start = (idx - 1) * 32 + 1
    return @view flv.data[start:start+31]
end

function ith_layer!(current_layer, queries_len, queries, proof)
    next_queries_len = 0
    i = 1

    while i <= queries_len
        query = queries[i]
        sibling = query ⊻ 1

        next_queries_len += 1
        queries[next_queries_len] = query >> 1

        if i == queries_len
            push!(proof, current_layer[sibling + 1])
            break
        end

        if query % 2 != 0
            push!(proof, current_layer[sibling + 1])
            i += 1
        else
            if queries[i + 1] != sibling
                push!(proof, current_layer[sibling + 1])
                i += 1
            else
                i += 2
            end
        end
    end

    return next_queries_len
end

struct BatchedMerkleProof
    proof::Vector{Vector{UInt8}}
end

sizeof(x::BatchedMerkleProof) = sum(length, x.proof)

function prove(tree::CompleteMerkleTree, queries)
    proof = Vector{Vector{UInt8}}()
    depth = length(tree.tree) - 1

    queries_buff = copy(queries)

    queries_buff .-= 1
    queries_cnt = length(queries)
    for i in 1:depth
        current_layer = FlatLayerView(tree.tree[i])
        queries_cnt = ith_layer!(current_layer, queries_cnt, queries_buff, proof)
    end

    return BatchedMerkleProof(proof)
end

function verify_ith_layer!(layer, queries, curr_cnt, proof, proof_cnt)
    next_cnt = 0
    i = 1

    while i <= curr_cnt
        query = queries[i]
        sibling = query ⊻ 1

        next_cnt += 1
        queries[next_cnt] = query >> 1

        if i == curr_cnt
            proof_cnt += 1
            pp = proof[proof_cnt]
            if query % 2 != 0
                layer[next_cnt] = hash_siblings(pp, layer[i])
            else
                layer[next_cnt] = hash_siblings(layer[i], pp)
            end
            break
        end

        if query % 2 != 0
            proof_cnt += 1
            pp = proof[proof_cnt]
            layer[next_cnt] = hash_siblings(pp, layer[i])
            i += 1
        else
            if queries[i + 1] != sibling
                proof_cnt += 1
                pp = proof[proof_cnt]
                layer[next_cnt] = hash_siblings(layer[i], pp)
                i += 1
            else
                layer[next_cnt] = hash_siblings(layer[i], layer[i + 1])
                i += 2
            end
        end
    end

    return (next_cnt, proof_cnt)
end

function verify(root::MerkleRoot, batched_proof::BatchedMerkleProof; depth, leaves, leaf_indices)    
    proof = copy(batched_proof.proof)
    layer = [hash_leaf(leaf) for leaf in leaves]
    queries = copy(leaf_indices)

    queries .-= 1

    curr_cnt = length(queries)
    proof_cnt = 0

    for _ in 1:depth
        (curr_cnt, proof_cnt) = verify_ith_layer!(layer, queries, curr_cnt, proof, proof_cnt)
    end

    return layer[1] == root.root
end

end # module MerkleTree