# In order to run this file, you need to follow a few steps:

# 1. Install the Julia package PyCall. When adding PyCall, it will by default use a miniconda environment
# You need to install dwave-ocean-sdk into that environment

# 2. (Only necessary if you need the clous services) Create an account on D-waves  website and get an API key
# You should then save the api key in a file named api_key.json in the same directory as this file. The file should have the following format:
# {
#     "api_key": "<your api key>",
# }

using PyCall

@pyinclude("src/GraphPartitioning/QA/dwave_call.py")

sample_qubo = pyimport("src.GraphPartitioning.QA.dwave_call").sample_qubo

qubo_matrix = [[-1 0 3]; [0 2 1]; [ 0 0 -0.5]]

sample_qubo(qubo_matrix)