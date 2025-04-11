module BinaryReedSolomon

using BinaryFields

# Utility functions
is_pow_2(n) = 2^(round(Int, log2(n))) == n


export fft!, ifft!, compute_twiddles, compute_twiddles!, eval_sk_at_vks, compute_pis, evaluate_basis
include("./binaryfft.jl")

export reed_solomon, encode, encode!, encode_non_systematic!, message_length, block_length, log_message_length, log_block_length
include("./reedsolomon.jl")

end # module BinaryReedSolomon
