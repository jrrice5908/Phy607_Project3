import matplotlib as plt
import numpy as np
from .system import System #Importing from system.py
from .particles import Particle #Importing from particle.py

# Simulation parameters
N = 1000                 # Number of rods
length = 100                 # Size of 2D box (L x L)
mean_speed = 0.2        # Mean speed of rods
speed_variability = 0.04 # Variability in speed
align_factor = 0.1       # Alignment factor, controls influence of neighbors
density_radius = 0.5     # Radius within which density affects alignment
attraction_radius = 2.0  # Radius within which rods are attracted to each other
attraction_strength = 0.01  # Strength of attraction force
theta_noise = 0.1        # Random noise added to orientation at each step
direction_change_prob = 0.1  # Probability of a direction flip per update
