import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import distance
import sys
import re
import argparse

class GeometryOptimizer:
    def __init__(self, num_vertices, vertices, ideal_distances, ideal_angles, k_edges, k_angle, Lid,
                 repulsion=False, conserve_membrane=False, save=True):
        self.num_vertices = num_vertices
        self.vertices = vertices
        self.ideal_distances = ideal_distances
        self.ideal_angles = ideal_angles
        self.k_edges = k_edges
        self.k_angle = k_angle
        self.current_step = 0
        self.repulsion = repulsion
        self.Lid = Lid
        self.conserve_membrane = conserve_membrane
        self.save = save

    def calculate_energy(self, vertices):
        vertices = np.array(vertices).reshape((self.num_vertices, 2))

        energy = 0.0

        if self.conserve_membrane == False:
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]
                energy += self.k_edges * (distance - ideal_distance) ** 2

        elif self.conserve_membrane == True:
            total_membrane = self.calc_total_membrane_area(self.ideal_distances)
            current_membrane = 0
            for i in range(self.num_vertices):
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]

                if i % 2 != 0:
                    current_membrane += distance
                elif i % 2 == 0:
                    energy += self.k_edges * (distance - ideal_distance) ** 2

            energy += self.k_edges * (current_membrane - total_membrane) ** 2

        if self.repulsion:
            for i in range(self.num_vertices):
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]

                if i % 2 != 0:
                    energy += self.repulsion_k(Lid=self.Lid) * (distance - self.Lid) ** 2

        for i in range(self.num_vertices):
            prev_vertex = (i - 1) % self.num_vertices
            next_vertex = (i + 1) % self.num_vertices
            angle = self.calculate_angle(vertices[prev_vertex], vertices[i], vertices[next_vertex])
            ideal_angle = self.ideal_angles[i]
            energy += self.k_angle * (angle - ideal_angle) ** 2

        return energy

    def calc_total_membrane_area(self, initial_distances):
        total_membrane = sum(initial_distances[1::2])
        return total_membrane

    def repulsion_k(self, Lid, R=7):
        h = 2 # nm
        a = 0.7 # nm
        xi = 2 # nm
        k = 0.8e-19 # J
        kt = 30e-3 # N/m
        depsilon = 4 # kT/nm

        A = (h/a)**2 * depsilon/np.sqrt(k*kt/1e18) * 4.11e-21
        k = 2*A* ( (np.sinh(self.Lid/xi)*np.cosh(self.Lid/xi)) * (1/np.tanh(Lid/xi) + 2/np.tanh(R/xi)) - (np.sinh(Lid/xi))**4 )\
                /(xi**2 * (1/np.tanh(Lid/xi)) + 2/np.tanh(R/xi))
        return k

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

# TODO: organize all these functions
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
    sys.stdout.write(f"The smallest radius that can fit all segments is {x_0}")
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
        extracted_list = [float(x) for x in extracted_list]
    else:
        sys.stdout.write("No list found in the string.")
        sys.stdout.write("Please look at input file structure.")
        sys.exit(1)
    return extracted_list

def calc_ideal_angle(L, R=7, xi=2):
    '''
    f_param = Delta_epsilon / sqrt(k*k_t) * h / a
    D_epsilon - energy diff of lipids on-top of protein and in membrane ~ 0.3
    h - monolayer thickness ~ 2 nm, a - lipid length ~ 0.7 nm
    k - monolayer rigidiy ~ 1e-9 J, k_t - tilt modulus ~ 30 mJ/m
    '''
    h = 2 # nm
    a = 0.7 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    depsilon = 4 # kT/nm
    f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
    return np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))

def caveolin_radius(L, R=7, xi=2):
    phi = calc_ideal_angle(L=L, R=R)
    R_c = (R+L*np.cos(np.pi-phi))/np.sin(np.pi-phi)
    return R_c

def calc_n_disks(L, R):
    Rc = caveolin_radius(L=L, R=R)
    circumference = 2*np.pi*Rc
    n_disks = circumference/(2*R + 2*L)
    return round(n_disks)

def calc_k(L, R=7, xi=2):
    h = 2 # nm
    a = 0.7 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def main(geometry_file, L, R, output_file, save=True, conserve_membrane=False, repulsion=False):
    n_disks = calc_n_disks(L=L, R=r_disk)

    k_edges = 100.0

    if repulsion:
        k_angle = calc_k(L=L, R=R, xi=2)
    elif not repulsion:
        k_angle = 1.0

    initial_geometry = np.loadtxt(args.inputfile)
    extracted_list = get_ideal_dist(geometry_file)
    sys.stdout.write(f"Initial configuration: {extracted_list}\n")
    ideal_distances = extracted_list

    id_angle = calc_ideal_angle(L, xi=2, R=7)

    ideal_angles = [id_angle] * len(ideal_distances)
    num_vertices = len(initial_geometry)
    sys.stdout.write(f"\nIdeal angle: {id_angle}")

    sys.stdout.write(f"""
N proteins: {n_disks}
Protein radius: {r_disk} nm
Half-distance: {L} nm
Estimated caveolin radius {(2*r_disk+L)*n_disks / (2 * np.pi)} nm\n""")
    if conserve_membrane:
        sys.stdout.write("Membrane treatment: changing total area conserved.\n")
    elif not conserve_membrane:
        # TODO: find a way to discern between case of uniform distribution and not conserving membrane
        sys.stdout.write("""Membrane treatment: static membrane.
    total area might be conserved but uniformly distributed, check generate_geom.py
    if that option was used.""")

    optimizer = GeometryOptimizer(num_vertices, initial_geometry.flatten(), ideal_distances, ideal_angles,
                                  k_edges, k_angle, Lid=L, repulsion=repulsion, conserve_membrane=conserve_membrane)
    optimized_vertices, min_energy = optimizer.optimize_geometry()
    optimized_vertices = optimized_vertices.reshape((num_vertices, 2))

    # TODO: maybe remove this print of angles
    sys.stdout.write(f"final angles:\n {calculate_angles(optimized_vertices)}\n")

    # important to see that edges don't break
    final_edges = calculate_edges(optimized_vertices)
    sys.stdout.write(f"final edges ({len(final_edges)}):\n {final_edges}\n")
    sys.stdout.write(f"total membrane eare: {sum(final_edges[1::2])}")
    if save == True:
        sys.stdout.write(f"writing final geometry to: {output_file}\n")
        with open(f"{output_file}", "w") as f:
            first_line = f"# combination: {final_edges}\n"
            f.write(first_line)
            np.savetxt(f, optimized_vertices)

    sys.exit(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="2D caveolin simulation optimizer",
        description="optimizer script for geometry",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-i", "--inputfile", required=True, help="path to input file")
    parser.add_argument("-o", "--outputfile", help="path to output file", default="geom_opt.txt")
    parser.add_argument("-s", "--save", action="store_true", help="default is to not save final geometry")
    parser.add_argument("-c", "--conserve_membrane", action="store_true", help="conserve membrane and allow distance between proteins to change")
    parser.add_argument("-r", "--repulsion", action="store_true", help="minimizes with protein repulsion")
    args = parser.parse_args()

    r_disk = 7
    L = 2 # half-distance between proteins

    main(geometry_file=args.inputfile, L=L, R=r_disk, save=args.save, conserve_membrane=args.conserve_membrane,
         output_file=args.outputfile, repulsion=args.repulsion)
