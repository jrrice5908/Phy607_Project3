import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from .system import System

class Particle:
	'''
	This defines how the particles in the systen behave. 
    
	Parameters:
	-----------
	- position: list for storing a particles position
	- speed: list for storing a particles 2D speed
	- orientation: angle of particles alignment
	'''

	def __init__(self, position, speed, orientation):
		self.position = np.array(position)
		self.speed = np.array(speed)
		self.orientation = orientation
		self.local_density_list = []
		self.particle_energy_list = []
	
        
	def find_distance(self, other):
		delta = np.abs(self - other)
		delta = np.where(delta > L / 2, L - delta, delta)
		return np.sqrt((delta ** 2).sum(axis=-1))
    
	def local_density(self):
		distances = distance(positions[idx], positions)
		neighbors = distances < density_radius
		return np.sum(neighbors) - 1  # exclude the rod itself
    
	def update_alignment(rod_orientation, mean_theta):
		#Aligning each particle to either mean_theta or 180 apart
		# Calculate the angular difference
		angle_diff = (mean_theta - rod_orientation + np.pi) % (2 * np.pi) - np.pi
		if np.abs(angle_diff) <= np.pi / 2:
			# Align directly to mean_theta if it's within 90 degrees
			return mean_theta
		else:
			# Otherwise, align 180 degrees away
			return (mean_theta + np.pi) % (2 * np.pi)
        
	def mean_orientation(idx):
		"""Compute mean orientation of neighbors within density radius."""
		distances = distance(positions[idx], positions)
		neighbors = distances < density_radius
		if np.sum(neighbors) > 1:
			return np.arctan2(
				np.sin(orientations[neighbors]).mean(),
				np.cos(orientations[neighbors]).mean()
			)
		return orientations[idx]
    
def attraction_force(idx):
	"""Compute the attractive force vector on a rod given by index `idx`."""
	distances = distance(positions[idx], positions)
	neighbors = (distances < attraction_radius) & (distances > 0)
    
	if np.sum(neighbors) > 0:
		# Calculate the displacement vectors toward each neighbor within attraction radius
		displacement_vectors = positions[neighbors] - positions[idx]
        
		# Apply periodic boundary adjustments to the displacement vectors
		displacement_vectors = np.where(displacement_vectors > L / 2, displacement_vectors - L, displacement_vectors)
		displacement_vectors = np.where(displacement_vectors < -L / 2, displacement_vectors + L, displacement_vectors)
        
		# Calculate the mean attraction vector
		attraction_vector = displacement_vectors.mean(axis=0)
        
		# Scale by attraction strength
		return attraction_strength * attraction_vector
	return np.array([0.0, 0.0])
