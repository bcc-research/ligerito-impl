module RAA

using BinaryFields, Base, Random
export encode, get_row, mat_a_t!, dot

function repeat!(message::Vector{T}, q::Int, repeated::Vector{T}) where T <: BinaryElem
    n = length(message)
    @assert length(repeated) == n * q "Length of repeated array must be n * q"

    for i in 1:n
        for j in 1:q
            repeated[q * (i - 1) + j] = message[i]
        end
    end
end

function mat_a!(codeword::Vector{T}) where T <: BinaryElem
    for i in 2:length(codeword)
        codeword[i] += codeword[i - 1]
    end
end

function encode(message::Vector{T}, q::Int, pi1::Vector{Int}, pi2::Vector{Int}) where T <: BinaryElem
    codeword = Vector{T}(undef, length(message) * q)

    repeat!(message, q, codeword)
    Base.permute!!(codeword, pi1) 
    mat_a!(codeword)
    Base.permute!!(codeword, pi2)
    mat_a!(codeword)

    return codeword
end

function mat_a_t!(row::Vector{Bool})
    for i in length(row)-1:-1:1
        row[i] ⊻= row[i + 1]
    end
end

function Fr_t!(row::Vector{Bool}, q::Int)
    n = length(row) ÷ q
    for i in 1:n
        row[i] = row[q * (i - 1) + 1]
        for j in 2:q
            row[i] ⊻= row[q * (i - 1) + j]
        end
    end
end

function a_t_pi1_composed(n::Int, q::Int, i::Int, pi1::Vector{Int})
    row = fill(true, n * q)
    for j in 1:(n * q)
        if pi1[j] > i
            row[j] = false
        end
    end 
    
    return row
end

function get_row(n::Int, q::Int, i::Int, pi1::Vector{Int}, pi2::Vector{Int})
    row = a_t_pi1_composed(n, q, i, pi2)
    mat_a_t!(row)
    Base.permute!!(row, pi1)
    Fr_t!(row, q)

    return row[1:n]
end

function dot(x::Vector{T}, row::Vector{Bool}) where T <: BinaryElem
    result = zero(T)
    n = length(row)
    for i in 1:n
        if row[i]
            result += x[i]
        end
    end

    return result
end

export test
function test(n::Int, q::Int)
    x = rand(BinaryElem32, n)

    pi1 = shuffle(1:n * q)
    pi2 = shuffle(1:n * q)

    pi1_inv = invperm(pi1)
    pi2_inv = invperm(pi2)

    codeword = encode(x, q, pi1, pi2)

    i = rand(1:n*q)
    @time row = get_row(n, q, i, pi1_inv, pi2_inv)

    @assert codeword[i] == dot(x, row) "Codeword and dot product do not match"
end

end