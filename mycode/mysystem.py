import numpy as np

class System:
   
    def __init__(self,num_parts, length, avg_speed, speed_var, density_rad, attract_rad, orientation_noise, beta, align_avg,align_stdev, attract_avg, attract_stdev):
        self.num_parts = num_parts
        self.length = length
        self.avg_speed = avg_speed
        self.speed_var = speed_var
        self.density_rad = density_rad
        self.attract_rad = attract_rad
        self.orientation_noise = orientation_noise
        self.beta = beta
        self.align_avg = align_avg
        self.align_stdev = align_stdev
        self.attract_avg = attract_avg
        self.attract_stdev = attract_stdev
                
    def init_align_attract(self):
        align_strength = np.random.normal(self.align_avg, self.align_stdev)
        attract_strength = np.random.normal(self.attract_avg, self.attract_stdev)
        return align_strength, attract_strength

    def init_states(self):
        positions = np.random.uniform(0, self.length, (self.num_parts, 2))
        orientations = np.random.uniform(0, 2 * np.pi, self.num_parts)
        return positions, orientations

    def init_speeds(self):
        speeds = np.random.normal(self.avg_speed, self.speed_var, self.num_parts)
        return speeds

    def dist_between(self, point1, point2):
        change = np.abs(point1 - point2)
        change = np.where(change > self.length / 2, self.length - change, change)
        return np.sqrt((change ** 2).sum(axis=-1))
         
    def align_E(self, i, orientation):
        positions, orientations = self.init_states()
        align_strength, attract_strength = self.init_align_attract()
        distances = self.dist_between(positions[i], positions)
        neighbors = distances < self.density_rad
        avg_theta = np.arctan2(
                        np.sin(orientations[neighbors]).mean(),
                        np.cos(orientations[neighbors]).mean()
                        )
        angle_diff = (orientation - avg_theta + np.pi) % (2 * np.pi) - np.pi
        return align_strength * (angle_diff ** 2)
    
    def attract_E(self, i, position):
        positions, orientations = self.init_states()
        align_strength, attract_strength = self.init_align_attract()
        distances = self.dist_between(position, positions)
        neighbors = (distances < self.attract_rad) & (distances > 0)
        attract_E = (distances[neighbors] ** 2).sum()
        return attract_strength * attract_E
