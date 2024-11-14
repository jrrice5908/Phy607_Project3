import numpy as np
from .particle import Particle

class System:
    '''
    This defines the box/boundary in which we are simulating. 
    
    Parameters:
    -----------
    - length (float): length of one side of the square box.
    - particles (list):list of all particles contained in box.
    '''
    
	def __init__(self, length, part_num):             
		self.length = length
		self.area = length ** 2
		self.part_num = part_num
		self.particles = []
		for part in range(self.part_num):
			part_i = Particle(
				np.random.uniform(0,length,2),
				np.random.uniform(0, 0.5),
				np.random.uniform(0,2 * np.pi))
			self.particles.append(part_i)
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
				
	def compute_energy(self)
		self.total_energy = 0
		for i in particles:
			for j in particles:
				d = i.distance(j)
				orientation_diff = np.cos(i.orientation - j.orientation)
				pair_energy = - alpha * orientation_diff - beta / (1 + d **2)
				self.total_energy += pair_energy / 2
	def mcmc_proposal(self):
		pass
	
	def mcmc_step(self):
		pass
			
        
	
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




            
