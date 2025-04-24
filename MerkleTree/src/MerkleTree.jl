module MerkleTree

using SHA
import Base: sizeof

export build_merkle_tree, get_root, get_depth, build_merkle_tree_fast
export CompleteMerkleTree, MerkleRoot, BatchedMerkleProof
export prove, verify

function threaded_each(nt::Int, xs::AbstractVector, f::Function; threshold::Int=256)
    n = length(xs)

    if n <= threshold || nt == 1
        @inbounds for i in 1:n
            f(xs, i)
        end
        return
    end

    chunk_size = ceil(Int, n / nt)
    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            start_idx = (t - 1) * chunk_size + 1
            end_idx = min(t * chunk_size, n)

            @inbounds for i in start_idx:end_idx
                f(xs, i)
            end
        end
    end
end

# for now simply assume that isbitstype(T) = true
hash_leaf(leaf) = sha256(reinterpret(UInt8, leaf))
hash_siblings(left_hash::Vector{UInt8}, right_hash::Vector{UInt8}) = sha256(vcat(left_hash, right_hash))

function is_power_of_two(n::Int)
    return n > 0 && (n & (n - 1)) == 0
end

struct CompleteMerkleTree
    tree::Vector{Vector{Vector{UInt8}}}
end

function build_merkle_tree_fast(leaves::AbstractVector)
    @assert is_power_of_two(length(leaves))
    if isempty(leaves)
        return []
    end

    nt = Threads.nthreads()
    n = length(leaves)

    layer0 = Vector{Vector{UInt8}}(undef, n)
    f = (curr, i) -> curr[i] = hash_leaf(leaves[i])
    threaded_each(nt, layer0, f)

    tree = [layer0]

    while length(tree[end]) > 1
        parent_layer = tree[end]
        parent_len = length(parent_layer)
        child_len = parent_len ÷ 2
        child_layer = Vector{Vector{UInt8}}(undef, child_len)

        if child_len <= 64 || nt == 1
            for i in 1:child_len
                left = parent_layer[2i - 1]
                right = parent_layer[2i]
                child_layer[i] = hash_siblings(left, right)
            end
        else 
            chunk_size = ceil(Int, child_len / nt)
            Threads.@sync for t in 1:nt
                Threads.@spawn begin
                    start_idx = (t - 1) * chunk_size + 1
                    end_idx = min(t * chunk_size, child_len)
        
                    @inbounds for i in start_idx:end_idx
                        left = parent_layer[2i - 1]
                        right = parent_layer[2i]
                        child_layer[i] = hash_siblings(left, right)
                    end
                end
            end
        end

        push!(tree, child_layer)
    end

    return CompleteMerkleTree(tree)
end


function build_merkle_tree(leaves::AbstractVector)
    @assert is_power_of_two(length(leaves))
    if isempty(leaves)
        return []
    end

    # current_layer = [hash_leaf(leaf) for leaf in leaves]
    current_layer = Vector{Vector{UInt8}}(undef, length(leaves))
    Threads.@threads for i in eachindex(leaves)
        current_layer[i] = hash_leaf(leaves[i])
    end
    tree = [current_layer]

    while length(current_layer) > 1
        next_layer_size = length(current_layer) ÷ 2
        next_layer = [Vector{UInt8}(undef, 32) for _ in 1:next_layer_size]

        threshold = 1024 # i'm just hardcoding this for now

        if next_layer_size < threshold
            for i in 1:next_layer_size
                left = current_layer[2i - 1]
                right = current_layer[2i]
                next_layer[i] = hash_siblings(left, right)
            end
        else
            Threads.@threads for i in 1:next_layer_size
                left = current_layer[2i - 1]
                right = current_layer[2i]
                next_layer[i] = hash_siblings(left, right)
            end
        end

        push!(tree, next_layer)
        current_layer = next_layer
    end

    return CompleteMerkleTree(tree)
end

struct MerkleRoot
    root::Vector{UInt8}
end

sizeof(x::MerkleRoot) = sizeof(x.root)
get_root(c::CompleteMerkleTree) = MerkleRoot(isempty(c.tree) ? nothing : c.tree[end][1])
get_depth(c::CompleteMerkleTree) = length(c.tree) - 1

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
        queries_cnt = ith_layer!(tree.tree[i], queries_cnt, queries_buff, proof)
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