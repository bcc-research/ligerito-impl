using BinaryFields
using MultilinearPoly

k = 5 

# randomness for partial evaluations
rs = rand(BinaryElem16, k)

f_evals = rand(BinaryElem16, 2^k)
f = MultiLinearPoly(f_evals)

g_evals = rand(BinaryElem16, 2^k)
g = MultiLinearPoly(g_evals)

h1 = sum(f.evals .* g.evals)

# start the prover 
prover = SumcheckProverInstance(f, g, h1)

fold!(prover, rs[1])
fold!(prover, rs[2])
fold!(prover, rs[3])
fold!(prover, rs[4])

tr = copy(prover.transcript)

verifier = SumcheckVerifierInstance(g, h1, tr)
fold!(verifier, rs[1])
fold!(verifier, rs[2])
fold!(verifier, rs[3])
fold!(verifier, rs[4])

f_eval = partial_eval(f, rs).evals[1]
verify(verifier, rs[5], f_eval)