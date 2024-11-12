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
    
    def __init__(self, length):             
        self.length = length
        self.area = length ** 2
        self.particles = []
        
    def add_particle(self, particle):               #System.add_particle, adds particle to class list
        #Add a particle to the simulation box
        self.particles.append(particle)
        
    def initialize_positions(self, N, length):      #System.initialize_positions, returns [position,[vx,vy]]
        particles = []
        for i in range(N):
            position = np.random.uniform(0, length, (N, 2))
            particle = Particle(position, np.zeros(2))    #Particles are initialized with speed=0
            particles.append(particle)
        return particles
        
    def initialize_speeds(self, mean_speed, speed_variability):     #System.initialize_speeds
        #Initialize particle velocities by sampling from Gaussian distribution centered around mean
        for particle in self.particles:
            particle.speed = np.random.normal(mean_speed, speed_variability, N)
            
    def intitialize_orientations(self, N):              #System.initialize_orientations
        for particle in self.particles:
            particle.orientation = np.random.uniform(0, 2 * np.pi, N)
            
