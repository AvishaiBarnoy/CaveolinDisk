import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import distance
import sys
import re

class GeometryOptimizer:
    def __init__(self, num_vertices, vertices, ideal_distances, ideal_angles, k_edges, k_angle):
        self.num_vertices = num_vertices
        self.vertices = vertices
        self.ideal_distances = ideal_distances
        self.ideal_angles = ideal_angles
        self.k_edges = k_edges
        self.k_angle = k_angle
        self.current_step = 0

    def calculate_energy(self, vertices):
        vertices = np.array(vertices).reshape((self.num_vertices, 2))

        energy = 0.0
        for i in range(self.num_vertices):
            next_vertex = (i + 1) % self.num_vertices
            distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
            ideal_distance = self.ideal_distances[i]
            energy += self.k_edges * (distance - ideal_distance) ** 2
        
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

    def optimize_geometry(self):
        result = minimize(self.calculate_energy, self.vertices, method='cg', options={'disp': True, 'maxiter': 25000})
        optimized_vertices = result.x.reshape((-1, 2))
        return optimized_vertices, result.fun

def calculate_angles(points):
    angles = []
    for i in range(len(points)):
        prev_point = points[i - 1]
        current_point = points[i]
        if i == len(points) - 1:
            next_point = points[0]
        else:
            next_point = points[i + 1]
        vector1 = prev_point - current_point
        vector2 = next_point - current_point
        dot_product = np.dot(vector1, vector2)
        magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        angle = np.arccos(dot_product / magnitude_product)
        angles.append(angle)
    return angles


def calculate_edges(points):
    distances = []
    for i in range(len(points)):
        current_point = points[i]
        if i == len(points) - 1:
            next_point = points[0]
        else:
            next_point = points[i + 1]
        dist = distance.euclidean(current_point, next_point)
        distances.append(dist)
    return distances

def calculate_circle_points(mod_length) -> list[tuple[float, float]]:
    min_radius = max(mod_length) / 2
    max_radius = sum(mod_length) / 2

    if min_radius > max_radius:
        return []
        raise ValueError("Segments are too long to fit in a circle")

    x_0 = (min_radius + max_radius) / 2
    d_0 = ((max_radius - min_radius) / 2) * 0.9

    while d_0 > 1e-9:
        current_angle = check_circle(mod_length, x_0)
        if current_angle < 2 * np.pi:
            x_0 -= d_0
        else:
            x_0 += d_0
        d_0 /= 1.8
    # print(f"The smallest radius that can fit all segments is {x_0}")
    angle = 0.0
    points = []
    for length in mod_length:
        x = x_0 * np.cos(angle)
        y = x_0 * np.sin(angle)
        points.append((x, y))
        angle += 2 * np.arcsin(length / (2 * x_0))

    return points

def check_circle(mod_length, radius: float) -> float:
    sum_of_angles = 0
    for length in mod_length:
        angle = 2 * np.arcsin(length / (2 * radius))
        sum_of_angles += angle
    return sum_of_angles

def get_ideal_dist(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    pattern = r"\[(.*?)\]"
    match = re.search(pattern, first_line)
    if match:
        extracted_list = match.group(1).split(', ')
        extracted_list = [int(x) for x in extracted_list]
        # print("Extracted list:", extracted_list)
    else:
        print("No list found in the string.")
    return extracted_list

def calc_ideal_angle(L, xi=2, R=7):
    '''
    f_param = Delta_epsilon / sqrt(k*k_t) * h / a
    D_epsilon - energy diff of lipids on-top of protein and in membrane
    h - monolayer thickness ~ 2 nm, a - lipid length ~ 1 nm
    k - monolayer rigidiy ~ 1e-9 J, k_t - tilt modulus ~ 3 mJ/m
    '''
    f_param = 0.3
    return f_param * 1 / (2/np.tanh(7/xi) + 1/np.tanh(L/xi))

if __name__ == "__main__":
    # Example usage:
    filename=sys.argv[1]
    initial_geometry = np.loadtxt(filename)

    k_edges = 100.0
    k_angle = 1.0

    extracted_list = get_ideal_dist(filename)
    print("Combination:", extracted_list)
    ideal_distances = extracted_list

    L = 2 # half-distance between proteins
    id_angle = np.pi - calc_ideal_angle(L, xi=2, R=7)

    ideal_angles = [id_angle] * len(ideal_distances)
    num_vertices = len(initial_geometry)
    print(f"Ideal angle: {id_angle}") 
    
    optimizer = GeometryOptimizer(num_vertices, initial_geometry.flatten(), ideal_distances, ideal_angles, k_edges, k_angle)
    optimized_vertices, min_energy = optimizer.optimize_geometry()
    optimized_vertices = optimized_vertices.reshape((num_vertices, 2))
    
    print("final angles:\n",calculate_angles(optimized_vertices))
    final_edges = calculate_edges(optimized_vertices)
    print(f"final edges ({len(final_edges)}):\n", final_edges)
    
    output_file = "geom_opt.txt"
    with open(f"{output_file}", "w") as f:
        first_line = f"# combination: {extracted_list}\n"  
        f.write(first_line)
        np.savetxt(f, optimized_vertices)
    
    visualize = False
    if visualize:
        plt.plot(*zip(*initial_geometry), 'o-', c='r', label='initia')
        plt.plot(*zip(*optimized_vertices), 'o-', c='b', label='final')
        plt.legend()
        plt.show()
