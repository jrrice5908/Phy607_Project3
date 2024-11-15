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
	- local_density_list: list for storing total denst=ity of each mcmc step
	- particle_energy_list: list for storing energies of particles for each mcmc step
	- new_pos: initializing xy positions as 0,0
	- new_orientation: initializing orientation as 0
	'''

	def __init__(self, position, orientation, speed):
		self.position = np.array(position)
		self.speed = np.array(speed)
		self.orientation = orientation
		self.local_density_list = []
		self.particle_energy_list = []
		self.new_pos = np.zeros(2)
		self.new_orientation = np.zeros(1)
        
	def find_distance(self, other):
		delta = np.sqrt((self.position[0] - other.position[0])**2+ (self.position[1] - other.position[1])**2)
		return delta

	def find_new_distance(self, other):
		delta = np.sqrt((self.new_pos[0] - other.new_pos[0])**2+ (self.new_pos[1] - other.new_pos[1])**2)
		return delta
