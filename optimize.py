import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial import distance
import sys
import re
import argparse

class GeometryOptimizer:
    def __init__(self, vertices, ideal_distances, k_edges, k_angle, Lid,
                 ideal_angles=[], repulsion=False, conserve_membrane=False, save=True):
        # TODO: work with numpy arrays instead of lists
        self.vertices = vertices
        self.num_vertices = len(self.vertices)//2
        self.ideal_distances = ideal_distances

        def initiate_ideal_angles(geometry):
            """takes initial geometry and initializes ideal angles list"""
            h = 2; a = 0.7; depsilon = 4; xi = 2; k =0.8e-19; kt = 30e-3; R = 7
            f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
            geometry = np.array(geometry).reshape((self.num_vertices, 2))
            cyclic_geometry = np.vstack([geometry, geometry[0]])
            angle_lst = []
            for i, j in enumerate(geometry):
                if i % 2 != 0:
                    L = distance.euclidean(cyclic_geometry[i], cyclic_geometry[i+1])
                    ideal_angle = np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
                    angle_lst.append(ideal_angle)
            return angle_lst
        self.ideal_angles = initiate_ideal_angles(self.vertices)

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

        # TODO: each energy step ideal angles list is recalculated and then used
        #       for energy calculation
        if self.conserve_membrane == False:
            # conserve_membrane name is misleading since membrane might be conserved (depending
            # on the options in generate_geom.py) it's just that distance between proteins/disks
            # is not allowed to change. # TODO: change conserve_membrane name
            # 
            # TODO: just use the initially calculated ideal_angle
            next_vertex = (i + 1) % self.num_vertices
            distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
            ideal_distance = self.ideal_distances[i]
            # make sure edges don't change
            energy += self.k_edges * (distance - ideal_distance) ** 2

        elif self.conserve_membrane == True:
            # if distance between proteins/disks allowed to change
            total_membrane = self.calc_total_membrane_area(self.ideal_distances)

            self.ideal_angles = self.update_ideal_angles(vertices)

            current_membrane = 0
            for i in range(self.num_vertices):
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]

                if i % 2 != 0:
                    current_membrane += distance
                elif i % 2 == 0:
                    energy += self.k_edges * (distance - ideal_distance) ** 2

            # constraint on total membrane, k_edges should be big 
            energy += self.k_edges * (current_membrane - total_membrane) ** 2

        if self.repulsion:
            # misleading name, takes into accound elastic energy
            for i in range(self.num_vertices):
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]

                if i % 2 != 0:
                    energy += self.elastic_energy(distance)

        for i in range(self.num_vertices):
            # iterates num_vertices to make sure num angles is computed correctly
            prev_vertex = (i - 1) % self.num_vertices
            next_vertex = (i + 1) % self.num_vertices
            angle = self.calculate_angle(vertices[prev_vertex], vertices[i], vertices[next_vertex])
            energy += self.k_angle * (angle - self.ideal_angles[i]) ** 2

        return energy

    def update_ideal_angles(self, geometry):
        # updates ideal angles from geometry, since id_angle is a function of L
        h = 2; a = 0.7; depsilon = 4; xi = 2; k =0.8e-19; kt = 30e-3; R = 7
        f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
        cyclic_geometry = np.vstack([geometry, geometry[0]])
        angle_lst = []
        for i, _ in enumerate(cyclic_geometry):
            if i % 2 != 0:
                L = distance.euclidean(cyclic_geometry[i], cyclic_geometry[i+1])
                ideal_angle = np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
                angle_lst.extend([ideal_angle, ideal_angle])
        return angle_lst

    def calc_total_membrane_area(self, initial_distances):
        """called one time to get total membrane area"""
        total_membrane = sum(initial_distances[1::2])
        return total_membrane

    def elastic_energy(self, L, R=7, h=2, a=0.7, xi=2, k=0.8e-19, kt=30e-3, depsilon=4):
        """elastic energy depending on L"""
        A = (h/a)**2 * depsilon/np.sqrt(k*kt/1e18) * 4.11e-21
        F = - A / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
        return F

    def calculate_angle(self, prev_point, current_point, next_point):
        """calculates single angle from geometry for deviation energy from ideal angle"""
        vector1 = prev_point - current_point
        vector2 = next_point - current_point
        dot_product = np.dot(vector1, vector2)
        magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        angle = np.arccos(dot_product / magnitude_product)
        return angle

    def optimize_geometry(self):
        """where the magic happens"""
        # TODO: add option for optimizer method and even optimization plan
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

def calc_ideal_angle(L, R=7, xi=2, h=2, a=0.7, k=0.8e-19, kt=30e-3, depsilon=4):
    # TODO: move into optimize
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

    k_edges = 100.0 # to keep proteins/disks rigid 

    # TODO: think if it makes sense to give k_angle during optimization
    if repulsion:
        k_angle = calc_k(L=L, R=R, xi=2)
    elif not repulsion:
        k_angle = 1.0

    initial_geometry = np.loadtxt(args.inputfile)
    extracted_list = get_ideal_dist(geometry_file)
    sys.stdout.write(f"Initial configuration: {extracted_list}\n")
    ideal_distances = extracted_list # TODO: rename and this line is redundant 

    # num_vertices = len(initial_geometry)

    sys.stdout.write(f"""
N proteins: {n_disks}
Protein radius: {r_disk} nm
Half-distance: {L} nm
Estimated caveolin radius {(2*r_disk+L)*n_disks / (2 * np.pi)} nm\n""")
    if conserve_membrane:
        sys.stdout.write("""Membrane treatment:
distance between proteins is allowed to change but total membrane is conserved\n""")
    elif not conserve_membrane:
        sys.stdout.write("""Membrane treatment:
total membrane is conserved but excesss membrane is distributed uniformly between proteins.\n""")

    # def __init__(self, vertices, ideal_distances, k_edges, k_angle, Lid, repulsion=False, conserve_membrane=False, save=True):
    # TODO: why is initial_geometry flattened and then reshaped?
    optimizer = GeometryOptimizer(vertices=initial_geometry.flatten(), ideal_distances=ideal_distances,
                                  k_edges=k_edges, k_angle=k_angle, Lid=L, repulsion=repulsion, conserve_membrane=conserve_membrane)
    optimized_vertices, min_energy = optimizer.optimize_geometry()
    optimized_vertices = optimized_vertices.reshape((len(initial_geometry), 2))

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
        description="optimizer script for geometry",
        epilog="Not all those who physics are lost"
    )
    io = parser.add_argument_group("input/output options")
    io.add_argument("-i", "--inputfile", required=True, help="path to input file")
    io.add_argument("-o", "--outputfile", help="path to output file", default="geom_opt.txt")
    io.add_argument("-s", "--save", action="store_true", help="default is to not save final geometry")

    options = parser.add_argument_group("optimization options")
    options.add_argument("-c", "--conserve_membrane", action="store_true", help="conserve membrane and allow distance between proteins to change")
    options.add_argument("-r", "--repulsion", action="store_true", help="minimizes with protein repulsion")

    args = parser.parse_args()

    r_disk = 7
    L = 2 # half-distance between proteins

    main(geometry_file=args.inputfile, L=L, R=r_disk, save=args.save, conserve_membrane=args.conserve_membrane,
         output_file=args.outputfile, repulsion=args.repulsion)
