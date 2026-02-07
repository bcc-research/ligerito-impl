using Test, Expander, BinaryFields
using BenchmarkTools

@testset "Emulate Encoding Benchmark" begin
    msg_len = 2^23
    d = 8
    
    # Prepare benchmark data
    (constraint, data) = emulate_expander(msg_len, d)
    out = zeros(UInt32, 128)
    
    # Run benchmark
    println("Benchmarking emulate_encoding with msg_len=$msg_len, d=$d")
    @btime emulate_encoding($constraint, $data, $out)
    
    # Test correctness
    @test length(out) == 128
    @test all(isa(x, UInt32) for x in out)
end

@testset "Prefetch L2 chunk" begin
    num_of_constraints = 2^20 
    msg_len = 2^23 
    d = 32
    num_of_chunks = div(msg_len, CHUNK_SIZE)

    # matrix = ExpanderMatrix(num_of_constraints, msg_len, num_of_chunks, d, CHUNK_SIZE, BinaryElem32)
    # msg = rand(BinaryElem32, msg_len + num_of_constraints)
    # @time encode!(msg, matrix)

    simple_expander = SimpleExpander(num_of_constraints, msg_len, d, BinaryElem32)
    msg2 = rand(BinaryElem32, msg_len + num_of_constraints)
    @time encode_simple_threads!(msg2, simple_expander)

    # nt = Threads.nthreads()
    # msgs = [rand(BinaryElem32, msg_len + num_of_constraints) for _ in 1:nt]
    # Threads.@sync for t in 1:nt
    #     Threads.@spawn begin
    #         # @time encode_simple!(msgs[t], simple_expander)
    #         @time encode!(msgs[t], matrix)
    #     end
    # end

    # @time begin
    #     for i in 1:CHUNK_SIZE
    #         prefetch(pointer(arr, i))
    #     end
    # end
end