#!/usr/bin/env python3

import os
import msprime

def sim(rr, n_reps, seed):
    return msprime.simulate(
                sample_size=100, Ne=10000, length=1e6,
                mutation_rate=1e-8, recombination_rate=rr,
                num_replicates=n_reps, random_seed=seed)

def dump(tslist, path):
    os.makedirs(path)
    for i, ts in enumerate(tslist):
        ts.dump(f"{path}/{i}.ts")

if __name__ == "__main__":
    n_reps = 5000
    dump(sim(1e-8, n_reps, 31415), "data/a")
    dump(sim(2e-8, n_reps, 27182), "data/b")
