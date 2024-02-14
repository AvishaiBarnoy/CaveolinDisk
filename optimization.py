import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.optimize import minimize
import sys

class Optimization:
    def __init__(self, vertices, ideal_distances, ideal_angles, k_angle, k_edges, max_steps=200, dt=0.1,
            optimizer="gd", use_adam_optimizer=False,
            log_energy=None, log_geom=None, n_geom=10):
        self.vertices = np.array(vertices, dtype=float)
        self.num_vertices = len(vertices)
        self.ideal_distances = np.array(ideal_distances, dtype=float)
        self.ideal_angles = np.array(ideal_angles, dtype=float)
        self.k_angle = k_angle
        self.k_edges = k_edges
        self.max_steps = max_steps
        self.dt = dt
        self.optimizer = optimizer
        self.current_step = 0
        self.log_energy = log_energy
        self.log_geom = log_geom
        self.n_geom = n_geom
        self.use_adam_optimizer = use_adam_optimizer
        if use_adam_optimizer:
            if self.optimizer.lower() != "gd":
                sys.exit("Adam optimizer requires gradient descent")
            self.beta1 = 0.9  # Exponential decay rates for moment estimates
            self.beta2 = 0.999
            self.epsilon = 1e-8  # Small constant to avoid division by zero
            self.m = np.zeros_like(vertices)  # First moment estimate
            self.v = np.zeros_like(vertices)  # Second moment estimate
            self.t = 0  # Time step

    def calculate_energy(self, vetrices):
        # Ensure vertices are reshaped properly
        vertices = np.array(self.vertices).reshape((self.num_vertices, 2))

        energy = 0.0
        # Calculate energy from spring distances
        distances = np.linalg.norm(np.roll(vertices, -1, axis=0) - vertices, axis=1)
        energy += self.k_edges * np.sum((distances - self.ideal_distances) ** 2)

        # Calculate energy from spring angles
        for i in range(self.num_vertices):
            prev_vertex = (i - 1) % self.num_vertices
            next_vertex = (i + 1) % self.num_vertices
            angle = self.calculate_angle(vertices[prev_vertex], vertices[i], vertices[next_vertex])
            ideal_angle = self.ideal_angles[i]
            energy += self.k_angle * (angle - ideal_angle) ** 2
        
        return energy

    def calculate_angle(self, prev_point, current_point, next_point):
        vector1 = prev_point - current_point
        vector2 = next_point - current_point
        dot_product = np.dot(vector1, vector2)
        magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        angle = np.arccos(dot_product / magnitude_product)
        return angle

    def bfgs_optimizer(self):
        result = minimize(self.calculate_energy, self.vertices.flatten(), method='BFGS', jac=False, options={'disp':True})
        # self.vertices = result.x.reshape((self.num_vertices, 2))
        self.vertices = result.x.reshape((self.num_vertices, 2))
        return self.vertices 

    def calculate_gradient(self):
        gradient = np.zeros_like(self.vertices)
        for i in range(self.num_vertices):
            for j in range(2):  # Loop over x and y components
                delta = 0.0001  # Small perturbation for numerical differentiation
                self.vertices[i, j] += delta
                energy_plus = self.calculate_energy(self.vertices)
                self.vertices[i, j] -= 2 * delta
                energy_minus = self.calculate_energy(self.vertices)
                self.vertices[i, j] += delta  # Restore original value
                gradient[i, j] = (energy_plus - energy_minus) / (2 * delta)
        return gradient

    def calculate_forces(self):
        forces = np.zeros_like(self.vertices)
        # Calculate forces from spring distances
        for i in range(self.num_vertices):
            next_vertex = (i + 1) % self.num_vertices
            distance_vector = self.vertices[next_vertex] - self.vertices[i]
            distance = np.linalg.norm(distance_vector)
            ideal_distance = self.ideal_distances[i]
            spring_force = self.k_edges * (distance - ideal_distance)
            forces[i] += spring_force * distance_vector / distance
            forces[next_vertex] -= spring_force * distance_vector / distance
        
        # Calculate forces from spring angles
        for i in range(self.num_vertices):
            prev_vertex = (i - 1) % self.num_vertices
            next_vertex = (i + 1) % self.num_vertices
            angle = self.calculate_angle(self.vertices[prev_vertex], self.vertices[i], self.vertices[next_vertex])
            ideal_angle = self.ideal_angles[i]
            diff_angle = angle - ideal_angle
            prev_vector = self.vertices[i] - self.vertices[prev_vertex]
            next_vector = self.vertices[next_vertex] - self.vertices[i]
            cross_product = np.cross(prev_vector, next_vector)
            torque = self.k_angle * diff_angle
            force1 = -torque * (prev_vector / np.linalg.norm(prev_vector)) * (1 / np.linalg.norm(prev_vector))
            force2 = torque * (next_vector / np.linalg.norm(next_vector)) * (1 / np.linalg.norm(next_vector))
            forces[i] += force1
            forces[next_vertex] += force2
        
        return forces
    

    def update_position(self):
        # print initial energy
        if self.current_step == 1:
            initial_energy = self.calculate_energy(self.vertices)
            sys.stdout.write(f"Initial energy: {initial_energy}\n")

        if self.optimizer.lower() == "gd":
            gradient = self.calculate_gradient()
            if self.use_adam_optimizer:
                self.t += 1
                self.m = self.beta1 * self.m + (1 - self.beta1) * gradient
                self.v = self.beta2 * self.v + (1 - self.beta2) * (gradient ** 2)
                m_hat = self.m / (1 - self.beta1 ** self.t)
                v_hat = self.v / (1 - self.beta2 ** self.t)
                self.vertices -= self.dt * m_hat / (np.sqrt(v_hat) + self.epsilon)
            else:
                self.vertices -= gradient * self.dt
        elif self.optimizer.lower() == "bfgs":
            self.vertices = self.bfgs_optimizer()
        elif optimizer.lower() == "forces":
            forces = self.calculate_forces()
            self.vertices += forces * self.dt
        else:
            sys.exit("Illegal optimizer")

        self.current_step += 1
        energy = round(self.calculate_energy(self.vertices), 8)
        grad = None
        # also log gradient
        if self.optimizer.lower()  == 'gd':
            grad = sum(self.calculate_gradient())
        # Log to file every n_geom steps, default 100
        if self.log_geom and self.current_step % self.n_geom == 0:
            with open(self.log_geom, 'a') as f:
                f.write(f"n={self.current_step}\t e={energy}\t grad={grad}\n")
                np.savetxt(f, self.vertices)
        # log step and energy
        if self.log_energy:
            with open(self.log_energy, 'a') as f:
                f.write(f"n={self.current_step}\t e={energy}\t grad={grad}\n")

        # print last energy value
        if self.current_step == self.max_steps:
            sys.stdout.write(f"Final energy: {energy}\n")

    def should_stop(self):
        if self.current_step >= self.max_steps:
            return True
        return False

    def log_std(self):
        if self.optimizer.lower() == "gd":
            if self.use_adam_optimizer:
                sys.stdout.write("Using Adam gradient descent\n")
            else:
                sys.stdout.write("Minimizing using gradient descent\n")
        elif self.optimizer.lower() == "bfgs":
            sys.stdout.write("Minimizing using BFGS optimizer\n")
        else:
            sys.stdout.write("Minimizing using force simulation\n")


class Animation:
    def __init__(self, sim):
        self.sim = sim
        self.fig, self.ax = plt.subplots()
        self.lines, = self.ax.plot([], [], 'bo-')
        self.energy_text = self.ax.text(0.02, 0.95, '', transform=self.ax.transAxes)
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-100, 100)
        self.ax.set_ylim(-100, 100)

    def init_animation(self):
        self.lines.set_data([], [])
        self.energy_text.set_text('')
        return self.lines, self.energy_text

    def animate(self, i):
        annotations = []
        for j, vertex in enumerate(self.sim.vertices):
            annotation = self.ax.annotate(str(j), (vertex[0], vertex[1]), xytext=(5,5),
                    textcoords="offset points", ha='center', color='black')
            annotations.append(annotation)

        if self.sim.should_stop():
            return [self.lines, self.energy_text] + annotations

        self.sim.update_position()
        self.lines.set_data(np.append(self.sim.vertices[:, 0], self.sim.vertices[0, 0]),
                            np.append(self.sim.vertices[:, 1], self.sim.vertices[0, 1]))

        energy_text = 'Step: {:d}   Energy: {:.6f}'.format(self.sim.current_step, self.sim.calculate_energy(
            self.sim.vertices))
        self.energy_text.set_text(energy_text)

        return [self.lines, self.energy_text] + annotations

    def start_animation(self):
        animation = FuncAnimation(self.fig, self.animate, init_func=self.init_animation, frames=self.sim.max_steps,
                interval=50, blit=True)
        plt.show()


if __name__ == "__main__":
    # Example Usage:
    vertices = [[-5, -5], [3, -6], [3, 2], [-2, 3]]  # Example vertices
    ideal_distances = [10, 4, 6, 4]  # Example ideal distances
    ideal_angles = [np.deg2rad(90), np.deg2rad(90), np.deg2rad(120), np.deg2rad(60)]  # Example ideal angles

    k_edges = 1
    k_angle = 1
    max_steps = 100
    dt = 0.05
    # optimizer options: gradient descent (gd), bfgs, forces
    #       gd can turn on adam optimizer flag
    #       TODO: if Adam without gd raise flag
    # optimizer = "bfgs"
    optimizer = "gd"
    use_adam_optimizer = False 
    log_energy = "results/log_energy.txt"
    log_geom = "results/log_geom.txt"
    n_geom = 10
    
    sim = Optimization(vertices, ideal_distances, ideal_angles, k_angle, k_edges, max_steps, dt,
            optimizer, use_adam_optimizer, log_energy, log_geom, n_geom)
    sim.log_std()

    # Toggle on/off the visualization by setting visualize=True or False
    visualize = True 

    if visualize:
        animation = Animation(sim)
        animation.start_animation()
    else:
        while not sim.should_stop():
            sim.update_position()
