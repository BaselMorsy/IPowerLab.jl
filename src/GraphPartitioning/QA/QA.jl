using PyCall

@pyinclude("src/GraphPartitioning/QA/dwave_call.py")

sample_qubo = pyimport("src.GraphPartitioning.QA.dwave_call").sample_qubo

qubo_matrix = [[-1 0 3]; [0 2 1]; [ 0 0 -0.5]]

sample_qubo(qubo_matrix)