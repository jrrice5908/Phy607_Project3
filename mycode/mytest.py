import numpy as np
from mymcmc import propose_new_state, mcmc_update_params, mcmc_step
from wmysystem import System 

params = System(100, 10, 0.2, 0.04, 0.5, 2.0, 0.1, 1.0, 0.1, 0.05, 0.01, 0.005)

for step in range(100):
	for i in range(params.num_parts):
		mcmc_step(i)
	mcmc_update_params()
	if step % 10 == 0:
		print(f"Step {step}: align_strength = {align_strength:.3f}, attraction_strength = {attraction_strength:.3f}")

