import numpy as np
import matplotlib.pyplot as plt
import particle as pp

class System:
	'''
	This defines the box/boundary in which we are simulating. 
	
	Parameters:
	-----------
	- length (float): length of one side of the square box.
	- particles (list): list of all particles contained in box.
	- part_num = (integer) number of particles
	- temp = (float) temperature of system
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
				
	def compute_energy(self, alpha = 1, beta = 25):
		#Compute energy at current step
		self.total_energy = 0
		for i in self.particles:
			for j in self.particles:
				d = i.find_distance(j)
				orientation_diff = np.cos(i.orientation - j.orientation)
				pair_energy = - alpha * orientation_diff**2 - beta / (1 + d **2)
				self.total_energy += pair_energy / 2
		return self.total_energy

	def compute_new_energy(self, alpha = 1, beta = 25):
		#Compute energy of proposed system
		self.new_total_energy = 0
		for i in self.particles:
			for j in self.particles:
				d = i.find_new_distance(j)
				orientation_diff = np.cos(i.new_orientation - j.new_orientation)
				pair_energy = - alpha * orientation_diff**2 - beta / (1 + d ** 2)
				self.new_total_energy += pair_energy / 2
		return self.new_total_energy

	def mcmc_proposal(self, position_step=2, orientation_step=2):
		#Proposed step for each particle position and orientation
		prior_energy = self.compute_energy()
		for part in self.particles:
			if np.random.rand() < 1:
				delta_pos = position_step * (np.random.rand(2) -0.5)
				part.new_pos = (part.position + delta_pos) % self.length
			if np.random.rand() < 1:
				delta_orientation = orientation_step * (np.random.rand() -0.5)
				part.new_orientation = (part.orientation + delta_orientation) % (2 * np.pi)
			
			if np.random.rand() < 0.0:
				part.new_orientation = (part.new_orientation + np.pi) % (2 * np.pi)				
	
	def metropolis_step(self):
		#Taking one Metropolis Hastings step
		#Use prposed energy to determine whether to take step
		self.mcmc_proposal()
		prior_energy = self.compute_energy()
		prop_energy = self.compute_new_energy()
		delta_e = prop_energy - prior_energy
		if delta_e < 0 or np.random.rand() < np.exp(-delta_e/ self.temp):
#			print("yes")
			for i in self.particles:
				i.position = i.new_pos
				i.orientation = i.new_orientation
				i.orientation_list.append(i.new_orientation)
				i.position_list.append(i.new_pos)
			return True
		else:
#			print("nah fam")
			return False

	def run_mcmc(self, iterations, burn_in_fraction =0.2):
		
		burn_in = int(iterations * burn_in_fraction)
		acceptance = 0
		E_list = []
		
		part_list = []
		
		for i in range(iterations):
			if self.metropolis_step():
				acceptance += 1
			E_list.append(self.total_energy)
			part_list.append(self.particles)
		if i >= burn_in:
			pass
			#part_list.append(self.particles)
		return E_list, part_list
					
	def plot_system(self):
		#Plotting system of N particles
		plt.figure(figsize= (6,6))
		plt.xlim(0,self.length)
		plt.ylim(0, self.length)
		x = []
		y = []
		dx = []
		dy = []
		for i in self.particles:
			x = (i.position[0])
			y = (i.position[1])
			dx = (np.cos(i.orientation) * 0.3)
			dy = (np.sin(i.orientation) * 0.3)
			plt.arrow(x,y,dx,dy, head_width =  -0.1, head_length = 0.2, fc = 'blue', ec = 'blue')
		plt.gca().set_aspect('equal', adjustable = 'box')
		plt.show()
			
	def compute_prior(self, particle):
		#Compute the prior probability for a particle's position and orientation
		pos_prior = np.exp(-np.sum((particle.position - self.length / 2) ** 2) / (2 * (self.length / 4) ** 2))
		orientation_prior = 1  # Uniform prior on orientation (default)
		return pos_prior * orientation_prior

	def compute_likelihood(self):
        #Compute likelihood using energy, lower energy configs have higher liklihood
		energy = self.compute_energy()
		return -energy / self.temp

	def compute_system_prior(self, particles):
		#Computing the prior for the system/configm, based on the spread of particles
		positions = [p.position for p in particles]
		mean_position = np.mean(positions, axis=0)
		spread = np.sum([np.linalg.norm(p.position - mean_position) for p in particles])
		#Using Gaussian dist because we want a smaller spread (want things clumped together
		return -spread

	def metropolis_step_with_prior(self):
		#Metropolis-Hastings step with priors
		prior_old = self.compute_system_prior(self.particles)
		likelihood_old = self.compute_likelihood()
		
		self.mcmc_proposal()
        	
		prior_energy = self.compute_energy()
		prop_energy = self.compute_new_energy()

		prior_new = self.compute_system_prior(self.particles)
		likelihood_new = self.compute_likelihood()

		#Want energy change to include priors
		delta_e = (prop_energy - prior_energy) + (prior_new - prior_old)

		if delta_e < 0 or np.random.rand() < np.exp(-delta_e / self.temp):
			for i in self.particles:
				i.position = i.new_pos
				i.orientation = i.new_orientation
				i.orientation_list.append(i.orientation)
				i.position_list.append(i.position)
			return True
		else:
			return False


	def run_mcmc_with_prior(self, iterations, burn_in_fraction=0.2):
        #Run MCMC with priors for a number of iterations
		burn_in = int(iterations * burn_in_fraction)
		acceptance = 0
		energy_list = []
		part_list = []
		for i in range(iterations):
			if self.metropolis_step_with_prior():
				acceptance += 1
			energy_list.append(self.total_energy)
			part_list.append(self.particles)
		print(f"Acceptance rate: {acceptance / iterations:.2%}")
		return energy_list, part_list

	def plot_system_with_prior(system):
        
    # Plot the particles
		for particle in system.particles:
			plt.scatter(particle.position[0], particle.position[1], color="red", s=30)
			plt.arrow(
				particle.position[0], 
				particle.position[1], 
				np.cos(particle.orientation) * 0.3, 
				np.sin(particle.orientation) * 0.3, 
				head_width=0.1, 
				head_length=0.15, 
				fc="blue", 
				ec="blue"
			)
    
		plt.xlim(0, system.length)
		plt.ylim(0, system.length)
		plt.title("Particle System with Prior Distribution")
		plt.xlabel("x")
		plt.ylabel("y")
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

	def autocorrelation_length(self, param):
		#Finding autocorrelation, will tell us about convergence
		param = np.array(param)
		n = len(param)
		lag = n // 2
		
		mean = np.mean(param)
		var = np.var(param)
	

		autocorr = []
		for i in range(lag):
			corr = np.sum((param[:n-lag] - mean) * (param[lag:] - mean)) / var / (n - lag)
			autocorr.append(corr)
	
		autocorrelation_len = 1 + 2 * np.sum(np.array(autocorr)[:1])

		return autocorrelation_len, autocorr


#Acutally running full sim:

system = System(length=10, part_num=100, temp=0.1)

system.plot_system_with_prior()
E,p = system.run_mcmc_with_prior(100000)
system.plot_system_with_prior()
auto_len, auto = system.autocorrelation_length(E)
print(auto_len)
plt.plot(system.particles[0].position_list)
plt.plot(system.particles[0].orientation_list)
plt.show()

plt.plot(E)
plt.show()

plt.plot(auto)
plt.show()
