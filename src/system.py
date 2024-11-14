import numpy as np
#from .particle import Particle
import particle as pp
class System:
	'''
	This defines the box/boundary in which we are simulating. 
	
	Parameters:
	-----------
	- length (float): length of one side of the square box.
	- particles (list):list of all particles contained in box.
	'''
	
	def __init__(self, length, part_num, temp):             
		self.length = length
		self.area = length ** 2
		self.part_num = part_num
		self.particles = []
		for part in range(self.part_num):
			part_i = pp.Particle(
				np.random.uniform(0,length,2),
				np.random.uniform(0, 2*np.pi),
				np.random.uniform(0, 0.5))
			self.particles.append(part_i)
		self.temp = temp
	def local_denisty_list(self):
		self.local_density_list = []
		for particle in self.particles:
			partx_den_list = []
			for particle2 in self.particles:
				distance_x = particle.find_distance(particle2)
				partx_den_list.append(distance_x)
			neighbors = distances < density_radius
			total_density = np.sum(neighbors) -1
			self.local_density_list.append(total_density)
	
	def mean_orientation_list(self):
		self.mean_orientation_list = []
		for particle in self.particles:
			part_or_list = []
			for particle2 in self.particles:
				distances = particle.find_distance(particle2)
				if (distances < density_radius):
					part_or_list.append(particle2.orientation)
			self.mean_orientation_list.append(np.mean(part_or_list))
				
	def compute_energy(self, alpha = 1.0, beta = 0.2):
		self.total_energy = 0
		for i in self.particles:
			for j in self.particles:
				d = i.find_distance(j)
				orientation_diff = np.cos(i.orientation - j.orientation)
				pair_energy = - alpha * orientation_diff - beta / (1 + d **2)
				self.total_energy += pair_energy / 2
		return self.total_energy

	def compute_new_energy(self, alpha = 1.0, beta = 0.2):
		self.new_total_energy = 0
		for i in self.particles:
			for j in self.particles:
				d = i.find_new_distance(j)
				orientation_diff = np.cos(i.new_orientation - j.new_orientation)
				pair_energy = - alpha * orientation_diff - beta / (1 + d ** 2)
				self.new_total_energy += pair_energy / 2
		return self.new_total_energy

	def mcmc_proposal(self, position_step=0.1, orientation_step=0.1):
		prior_energy = self.compute_energy()
		for part in self.particles:
			if np.random.rand() < 1:
				delta_pos = position_step * (np.random.rand(2) -0.5)
				part.new_pos = (part.position + delta_pos) % self.length
			if np.random.rand() < 1:
				delta_orientation = orientation_step * (np.random.rand() -0.5)
				part.new_orientation = (part.orientation + delta_orientation) % (2 * np.pi)
			
			if np.random.rand() < 0.1:
				part.new_orientation = (part.new_orientation + np.pi) % (2 * np.pi)
		
				
	
	def mcmc_step(self):
		self.mcmc_proposal()
		prior_energy = self.compute_energy()
		prop_energy = self.compute_new_energy()
		delta_e = prop_energy - prior_energy
		if delta_e < 0 or np.random.rand() < np.exp(-delta_energy/ self.temp):
			for i in self.particles:
				i.position = i.new_pos
				i.orientation = i.new_orientation

system = System(10, 100,0.1)
system.mcmc_step()

#    def add_particle(self, particle):               #System.add_particle, adds particle to class list
#        #Add a particle to the simulation box
#        self.particles.append(particle)
        
#    def initialize_positions(self, N, length):      #System.initialize_positions, returns [position,[vx,vy]]
#        particles = []
#        for i in range(N):
#            position = np.random.uniform(0, length, (N, 2))
#            particle = Particle(position, np.zeros(2))    #Particles are initialized with speed=0
#            particles.append(particle)
#        return particles
        
#    def initialize_speeds(self, mean_speed, speed_variability):     #System.initialize_speeds
#        #Initialize particle velocities by sampling from Gaussian distribution centered around mean
#        for particle in self.particles:
#            particle.speed = np.random.normal(mean_speed, speed_variability, N)
            
#    def intitialize_orientations(self, N):              #System.initialize_orientations
#        for particle in self.particles:
#            particle.orientation = np.random.uniform(0, 2 * np.pi, N)




            
