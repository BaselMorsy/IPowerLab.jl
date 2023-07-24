import dimod 
import json
from dwave.system import DWaveSampler, LeapHybridSampler, EmbeddingComposite

def sample_qubo(
        qubo_matrix,
        num_reads = 10,
        solver = "simulated_annealing",
        return_type = "list"
):
    bqm = dimod.BinaryQuadraticModel(
        qubo_matrix,
        "BINARY"
    )

    if solver == "qpu":
        sampler = EmbeddingComposite(DWaveSampler(
            token = json.load(open("src/GraphPartitioning/QA/api_key.json", "r"))["api_key"]
        ))
    
    elif solver == "hybrid":
        sampler = LeapHybridSampler(
            token = json.load(open("src/GraphPartitioning/QA/api_key.json", "r"))["api_key"]
        )
    
    elif solver == "simulated_annealing":
        sampler = dimod.SimulatedAnnealingSampler()
    
    elif solver == "random":
        sampler = dimod.RandomSampler()
    
    else:
        raise ValueError(f"Invalid solver {solver}")
    
    if solver != "hybrid":

        sample_set = sampler.sample(
            bqm,
            num_reads = num_reads,
        ).aggregate()
    
    else:
        sample_set = sampler.sample(
            bqm,
        ).aggregate()

    if return_type == "list":
        return  [(sample, energy, counts) for sample, energy, counts in zip(sample_set.record.sample, sample_set.record.energy, sample_set.record.num_occurrences)]
    else:
        return sample_set