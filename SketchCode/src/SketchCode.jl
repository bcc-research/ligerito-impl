module SketchCode

using BinaryFields, Expander, BinaryReedSolomon, MerkleTree

function sketch_encode!(c::AbstractArray{T}, expander::ExpanderMatrix{T}, rs) where T <: BinaryElem
    Expander.encode!(c, expander) 
    BinaryReedSolomon.encode_non_systematic!(rs, c[expander.msg_len+1:end])
end

function poly2mat(poly::Vector{T}, m::Int, n::Int, expander_compression::Int) where {T <: BinaryElem}
    reed_solomon_rate = 1/2
    m_target = m + Int(div(m, expander_compression * reed_solomon_rate))
    mat = zeros(T, m_target, n)

    nt = Threads.nthreads()
    chunk_size = ceil(Int, n / nt)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            start_col = (t-1)*chunk_size + 1
            end_col = min(t*chunk_size, n)

            @inbounds for j in start_col:end_col
                for i in 1:m
                    mat[i, j] = poly[(j-1)*m + i]
                end
            end
        end
    end

    return mat
end

function encode_poly!(poly_mat::AbstractMatrix{T}, expander::ExpanderMatrix{T}, rs; parallel=true) where T <: BinaryElem
    n = size(poly_mat, 2)

    # first we encode each column with sketched code
    if parallel
        nt = Threads.nthreads()
        chunk_size = ceil(Int, n / nt)

        Threads.@sync for t in 1:nt
            Threads.@spawn begin
                start_col = (t-1)*chunk_size + 1
                end_col = min(t*chunk_size, n)

                @inbounds for j in start_col:end_col
                    sketch_encode!(@view(poly_mat[:, j]), expander, rs)
                end
            end
        end
    else
        @inbounds for j in 1:n
            sketch_encode!(@view(poly_mat[:, j]), expander, rs)
        end
    end

    leaves = eachrow(poly_mat)
    @time tree = build_merkle_tree_odd(leaves)

    return tree
end

function sketch_code(::Type{T}, m::Int, expander_compression::Int) where T <: BinaryElem
    inv_rs_rate = 2
    expander_syndrome_len = div(m, expander_compression)
    D = 8

    expander = ExpanderMatrix(expander_syndrome_len, m, D, T)
    rs = reed_solomon(T, expander_syndrome_len, expander_syndrome_len * inv_rs_rate)
    return (expander, rs)
end

function get_things()
    n = 2^30 
    m = 2^23 
    k = 2^7 

    expander_compression = 8 
    T = BinaryElem32 

    poly = rand(T, n)
    (expander, rs) = sketch_code(T, m, expander_compression)

    poly_mat = poly2mat(poly, m, k, expander_compression)
    return (poly_mat, expander, rs)
end

export poly2mat, sketch_code, encode_poly!
export get_things

end # module SketchCode
