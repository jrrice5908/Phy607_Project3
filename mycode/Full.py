import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
N = 100                 # Number of rods
L = 100                 # Size of 2D box (L x L)
mean_speed = 0.2        # Mean speed of rods
speed_variability = 0.04 # Variability in speed
density_radius = 0.5     # Radius within which density affects alignment
attraction_radius = 2.0  # Radius within which rods are attracted to each other
theta_noise = 0.1        # Random noise added to orientation at each step
beta = 1.0               # Inverse temperature for MCMC

# Initialize positions, orientations (angles), and speeds of rods
positions = np.random.uniform(0, L, (N, 2))
orientations = np.random.uniform(0, 2 * np.pi, N)
speeds = np.random.normal(mean_speed, speed_variability, N)

# Initialize priors for alignment and attraction strength
align_strength_prior_mean = 0.1
align_strength_prior_std = 0.05
attraction_strength_prior_mean = 0.01
attraction_strength_prior_std = 0.005

# Initialize the parameters (alignment and attraction strength)
align_strength = np.random.normal(align_strength_prior_mean, align_strength_prior_std)
attraction_strength = np.random.normal(attraction_strength_prior_mean, attraction_strength_prior_std)

# Distance calculation function
def distance(p1, p2):
    delta = np.abs(p1 - p2)
    delta = np.where(delta > L / 2, L - delta, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

# Define energy functions for the likelihood
def alignment_energy(idx, orientation):
    distances = distance(positions[idx], positions)
    neighbors = distances < density_radius
    mean_theta = np.arctan2(
        np.sin(orientations[neighbors]).mean(),
        np.cos(orientations[neighbors]).mean()
    )
    angle_diff = (orientation - mean_theta + np.pi) % (2 * np.pi) - np.pi
    return align_strength * (angle_diff ** 2)

def attraction_energy(idx, position):
    distances = distance(position, positions)
    neighbors = (distances < attraction_radius) & (distances > 0)
    attraction_energy = (distances[neighbors] ** 2).sum()
    return attraction_strength * attraction_energy

# Propose new state for MCMC updates
def propose_new_state(idx):
    new_orientation = (orientations[idx] + theta_noise * np.random.uniform(-1, 1)) % (2 * np.pi)
    displacement = speeds[idx] * np.array([np.cos(new_orientation), np.sin(new_orientation)])
    new_position = (positions[idx] + displacement) % L
    return new_position, new_orientation

# MCMC for Bayesian parameter inference
def mcmc_update_params():
    global align_strength, attraction_strength
    
    # Propose new parameters
    new_align_strength = align_strength + np.random.normal(0, 0.01)
    new_attraction_strength = attraction_strength + np.random.normal(0, 0.001)
    
    # Calculate the prior probabilities for the current and proposed parameters
    current_prior = (np.exp(-0.5 * ((align_strength - align_strength_prior_mean) / align_strength_prior_std) ** 2) *
                     np.exp(-0.5 * ((attraction_strength - attraction_strength_prior_mean) / attraction_strength_prior_std) ** 2))
    new_prior = (np.exp(-0.5 * ((new_align_strength - align_strength_prior_mean) / align_strength_prior_std) ** 2) *
                 np.exp(-0.5 * ((new_attraction_strength - attraction_strength_prior_mean) / attraction_strength_prior_std) ** 2))

    # Calculate the likelihood for the current and proposed parameters
    current_likelihood = sum(alignment_energy(i, orientations[i]) + attraction_energy(i, positions[i]) for i in range(N))
    new_likelihood = sum(alignment_energy(i, orientations[i]) + attraction_energy(i, positions[i]) for i in range(N))
    
    # Metropolis-Hastings criterion
    acceptance_ratio = np.exp(-beta * (new_likelihood - current_likelihood)) * (new_prior / current_prior)
    if np.random.rand() < acceptance_ratio:
        align_strength = new_align_strength
        attraction_strength = new_attraction_strength

# MCMC update for rods
def mcmc_step(idx):
    global positions, orientations

    current_energy = alignment_energy(idx, orientations[idx]) + attraction_energy(idx, positions[idx])
    new_position, new_orientation = propose_new_state(idx)
    proposed_energy = alignment_energy(idx, new_orientation) + attraction_energy(idx, new_position)

    if np.random.rand() < np.exp(-beta * (proposed_energy - current_energy)):
        positions[idx] = new_position
        orientations[idx] = new_orientation

# Visualization function
def plot_rods():
    plt.figure(figsize=(6, 6))
    plt.xlim(0, L)
    plt.ylim(0, L)
    for i in range(N):
        x, y = positions[i]
        dx = np.cos(orientations[i]) * 0.3
        dy = np.sin(orientations[i]) * 0.3
        plt.arrow(x, y, dx, dy, head_width=0.1, head_length=0.2, fc='blue', ec='blue')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# Main loop
for step in range(100):
    for i in range(N):
        mcmc_step(i)
    mcmc_update_params()
    if step % 10 == 0:
        print(f"Step {step}: align_strength = {align_strength:.3f}, attraction_strength = {attraction_strength:.3f}")
        plot_rods()

