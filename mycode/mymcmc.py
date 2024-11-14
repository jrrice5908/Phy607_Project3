import numpy as np
from wmysystem import System

params = System(100, 10, 0.2, 0.04, 0.5, 2.0, 0.1, 1.0, 0.1, 0.05, 0.01, 0.005)

positions, orientations = params.init_states()
align_strength, attract_strength = params.init_align_attract()

def propose_new_state(i):
	new_orientation = (orientations[i] + params.orientation_noise * np.random.uniform(-1, 1)) % (2 * np.pi)
	displacement = params.init_speeds(i) * np.array([np.cos(new_orientation), np.sin(new_orientation)])
	new_position = (positions[i] + displacement) % System.length
	return new_position, new_orientation

def mcmc_update_params():
# Propose new parameters
	new_align_strength = align_strength + np.random.normal(0, 0.01)
	new_attraction_strength = attraction_strength + np.random.normal(0, 0.001)

# Calculate the prior probabilities for the current and proposed parameters
	current_prior = (np.exp(-0.5 * ((align_strength - params.align_avg) / params.align_stdev) ** 2) *
			np.exp(-0.5 * ((attraction_strength - params.attract_avg) / params.attract_stdev) ** 2))
	new_prior = (np.exp(-0.5 * ((new_align_strength - params.align_avg) / params.align_stdev) ** 2) *
			np.exp(-0.5 * ((new_attraction_strength - params.attract_avg) / params.attract_stdev) ** 2))

# Calculate the likelihood for the current and proposed parameters
	current_likelihood = sum(params.align_E(i, orientations[i]) + params.attract_E(i, positions[i]) for i in range(N))
	new_likelihood = sum(params.align_E(i, orientations[i]) + params.attract_E(i, positions[i]) for i in range(N))

# Metropolis-Hastings criterion
	acceptance_ratio = np.exp(-System.beta * (new_likelihood - current_likelihood)) * (new_prior / current_prior)
	if np.random.rand() < acceptance_ratio:
		align_strength = new_align_strength
		attract_strength = new_attract_strength

def mcmc_step(i):
	current_energy = params.align_E(i, orientations[i]) + params.attract_E(i, positions[i])
	new_position, new_orientation = propose_new_state(i)
	proposed_energy = params.align_E(i, new_orientation) + params.attract_E(i, new_position)

	if np.random.rand() < np.exp(-params.beta * (proposed_energy - current_energy)):
		positions[i] = new_position
		orientations[i] = new_orientation
