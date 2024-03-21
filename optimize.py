import numpy as np
from scipy.optimize import minimize
from scipy.spatial import distance
import sys
import re
import argparse

class GeometryOptimizer:
    def __init__(self, vertices, ideal_distances, k_edges, k_angle,
                 optimizer, n_steps=5000,
                 energy_method='new',
                 ideal_angles=[], repulsion=False, conserve_membrane=False, save=True):
        # TODO: decide which parameter are input and which are initialized inside
        # TODO: work with numpy arrays instead of lists
        self.vertices = vertices
        self.num_vertices = len(self.vertices)//2
        self.ideal_distances = ideal_distances

        self.energy_method = energy_method.lower() # TODO: implement in the future use choices

        # assumes initial geoemtry is close to optimal one
        if self.energy_method == "old":
            def initiate_ideal_angles(geometry):
                # TODO: move to energy calculation for energy_method!
                """takes initial geometry and initializes ideal angles list"""
                h = 2; a = 0.7; depsilon = 4; k = 0.4e-19; kt = 20e-3; R = 7
                xi = np.sqrt(k/kt) * 1e9 # J / nm
                f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
                geometry = np.array(geometry).reshape((self.num_vertices, 2))
                cyclic_geometry = np.vstack([geometry, geometry[0]])
                angle_lst = []
                for i, _ in enumerate(geometry):
                    if i % 2 != 0:
                        L = distance.euclidean(cyclic_geometry[i], cyclic_geometry[i+1])
                        ideal_angle = np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
                        angle_lst.append(ideal_angle)
                return angle_lst

            # allows for external input of ideal angle
            if ideal_angles:
                self.ideal_angles = ideal_angles
            elif not ideal_angles:
                self.ideal_angles = initiate_ideal_angles(self.vertices)

        self.k_edges = k_edges
        self.k_angle = k_angle
        self.repulsion = repulsion
        self.conserve_membrane = conserve_membrane
        self.save = save


        # TODO: add optimization plan -> mv todo to input file
        self.opt_method = optimizer # TODO: implement more methods, currently only: cg, bfgs, l-bfgs-b
        self.n_steps = n_steps

    def calculate_energy(self, vertices):
        # TODO: add documentation

        # TODO: implement differnet energy calculation methods
        #   1. write energy functions for new method
        #   2. identify which function calls are for both methods
        #   3. write new logic

        vertices = np.array(vertices).reshape((self.num_vertices, 2))
        current_angles = self.calculate_angles(vertices)
        energy = 0.0

        if self.energy_method == 'new': # default should be Misha's new energy function
            # iterate over length and angles
            # for each length check two angles
            L_lst = self.calc_L_lst(vertices)
            for i,_ in enumerate(L_lst):
                energy += self.calc_new_energy(L=L_lst[i], phi=current_angles[i])

            for i in range(self.num_vertices):
                # TODO: use this logic for constraints and edge calculation, maybe refactor
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]

                if i % 2 == 0:
                    # constraint only on proteins/disks
                    energy += self.k_edges * (distance - ideal_distance) ** 2

            # constraint on total membrane, k_edges should be big
            total_membrane = self.calc_total_membrane_area(self.ideal_distances)
            current_membrane = sum(L_lst[0::2])
            energy += self.k_edges * (current_membrane - total_membrane) ** 2

        elif self.energy_method == 'old':
            if selfconserve_membrane == False:
                # conserve_membrane name is misleading since membrane might be conserved (depending
                # on the options in generate_geom.py) it's just that distance between proteins/disks
                # is not allowed to change. # TODO: change conserve_membrane name
                next_vertex = (i + 1) % self.num_vertices
                distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i]
                # make sure edges don't change
                energy += self.k_edges * (distance - ideal_distance) ** 2

            elif self.conserve_membrane == True:
                # logic if distance between proteins/disks allowed to change
                # TODO: mv total membrane to parameters and calculate it there if input has relevant flags
                total_membrane = self.calc_total_membrane_area(self.ideal_distances)

                self.ideal_angles = self.update_ideal_angles(vertices)

                current_membrane = 0
                for i in range(self.num_vertices):
                    # TODO: use this logic for constraints and edge calculation, maybe refactor
                    next_vertex = (i + 1) % self.num_vertices
                    distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                    ideal_distance = self.ideal_distances[i]

                    if i % 2 != 0:
                        # sums for total membrane conservation
                        current_membrane += distance
                    elif i % 2 == 0:
                        # constraint only on proteins/disks
                        energy += self.k_edges * (distance - ideal_distance) ** 2

                # constraint on total membrane, k_edges should be big
                energy += self.k_edges * (current_membrane - total_membrane) ** 2

            if self.repulsion:
                # misleading name, takes into accound elastic energy # TODO: rename
                for i in range(self.num_vertices):
                    next_vertex = (i + 1) % self.num_vertices
                    distance = np.linalg.norm(vertices[next_vertex] - vertices[i])
                    ideal_distance = self.ideal_distances[i]

                    if i % 2 != 0:
                        energy += self.elastic_energy(distance)

            # sums energy price due to deviation from ideal angle
            for i in range(self.num_vertices):
                prev_vertex = (i - 1) % self.num_vertices
                next_vertex = (i + 1) % self.num_vertices

                # TODO: call from current_angles array
                angle = self.calculate_angle(vertices[prev_vertex], vertices[i], vertices[next_vertex])
                energy += self.k_angle * (angle - self.ideal_angles[i]) ** 2

        return energy

    def calc_L_lst(self, geometry):
        L_lst = []
        for i,j in enumerate(geometry):
            if i % 2 != 0:
                next_vertex = (i + 1) % len(geometry)
                distance = np.linalg.norm(geometry[next_vertex] - geometry[i])
                L_lst.extend([distance,distance])
        return np.array(L_lst)

    def calc_new_energy(self, L, phi, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
        """
        energy calculation for each L and phi
        """
        # constants
        xi = np.sqrt(k/kt) * 1e9 # J / nm
        r = 1/np.tanh(R/xi)
        K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
        e = h/a * depsilon # kT

        # dependency on L
        l = 1/np.tanh(L/xi)

        # convert phi according to convention
        phi = np.pi - phi

        f1 = 0.5 * K * l * (2*r + l) / (r + l) * phi**2
        f2 = l / (r + l) * e * phi
        f3 = 0.5 / (r + l) * e**2 / K

        F = f1 - f2 - f3
        return F

    def calculate_angles(self, points):
        angles = []
        for i in range(len(points)):
            # TODO: change to % i and enumerate for better readability
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
        return np.array(angles)

    def update_ideal_angles(self, geometry):
        """
        updates ideal angles from geometry, since id_angle is a function of L
        """
        h = 2; a = 0.7; depsilon = 4; k = 0.4e-19; kt = 20e-3; R = 7
        xi = np.sqrt(k/kt) * 1e9 # J / nm
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
        """
        called one time to get total membrane area
        """
        total_membrane = sum(initial_distances[1::2])
        return total_membrane

    def elastic_energy(self, L, R=7, h=2, a=0.7, k=0.4e-19, kt=20e-3, depsilon=4):
        """
        elastic energy depending on L
        """
        xi = np.sqrt(k/kt) * 1e9 # J / nm
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
        result = minimize(self.calculate_energy, self.vertices, method=self.opt_method, options={'disp': True, 'maxiter': self.n_steps})
        optimized_vertices = result.x.reshape((-1, 2))
        return optimized_vertices, result.fun

# TODO: organize all these functions

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

def calc_ideal_angle(L, R=7, h=2, a=0.7, k=0.4e-19, kt=20e-3, depsilon=4):
    # TODO: move into optimize
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
    return np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))

def caveolin_radius(L, R=7,k=0.4e-19, kt=20e-3):
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    phi = calc_ideal_angle(L=L, R=R)
    R_c = (R+L*np.cos(np.pi-phi))/np.sin(np.pi-phi)
    return R_c

def calc_n_disks(L, R):
    Rc = caveolin_radius(L=L, R=R)
    circumference = 2*np.pi*Rc
    n_disks = circumference/(2*R + 2*L)
    return round(n_disks)

def calc_k(L, R=7):
    h = 2 # nm
    a = 0.7 # nm
    k = 0.4e-19 # J
    kt = 20e-3 # N/m
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def main(geometry_file, L, R, output_file, save=True, conserve_membrane=False, repulsion=False,
         optimizer='cg', n_steps=5000, energy_method='new'):
    n_disks = calc_n_disks(L=L, R=r_disk)

    k_edges = 100.0 # to keep proteins/disks rigid 

    if repulsion:
        k_angle = calc_k(L=L, R=R)
    elif not repulsion:
        # if elastic energy is not taken into account optimization works only on angles
        #   and the value of k_angle is not important just needs to be much smaller than k_edges
        k_angle = 1.0

    initial_geometry = np.loadtxt(args.inputfile)
    ideal_distances = get_ideal_dist(geometry_file) # TODO: rename, name misleading due to exapnsion of optimization options
    sys.stdout.write(f"Initial configuration: {ideal_distances}\n")

    sys.stdout.write(f"""
N proteins: {n_disks}
Protein radius: {r_disk} nm
Half-distance: {L} nm
Estimated caveolin radius {(2*r_disk+L)*n_disks / (2 * np.pi)} nm

minimization method: {optimizer}
energy calculation method: {energy_method}\n""")

    if conserve_membrane:
        sys.stdout.write("""Membrane treatment:
distance between proteins is allowed to change but total membrane is conserved\n""")
    elif not conserve_membrane:
        sys.stdout.write("""Membrane treatment:
total membrane is conserved but excesss membrane is distributed uniformly between proteins.\n""")

    # TODO: why is initial_geometry flattened and then reshaped? -> input into scipy.minimize should be flat but rather than that...
    optimizer = GeometryOptimizer(vertices=initial_geometry.flatten(), ideal_distances=ideal_distances,
                                  optimizer=optimizer, n_steps=n_steps, energy_method='new',
                                  k_edges=k_edges, k_angle=k_angle, repulsion=repulsion, conserve_membrane=conserve_membrane)
    optimized_vertices, _ = optimizer.optimize_geometry() # _ is final energy but is unassigned 'cause it's unused
    optimized_vertices = optimized_vertices.reshape((len(initial_geometry), 2))

    # TODO: inisde optimize print angles after minimize is called, maybe in optimize_geometry()
    # sys.stdout.write(f"final angles:\n {calculate_angles(optimized_vertices)}\n")

    # important to see that edges don't break
    final_edges = calculate_edges(optimized_vertices)
    sys.stdout.write(f"final edges ({len(final_edges)}):\n {final_edges}\n")
    sys.stdout.write(f"total membrane area: {sum(final_edges[1::2])}\n")
    if save == True:
        sys.stdout.write(f"""writing final geometry to: {output_file}\n""")
        with open(f"{output_file}", "w") as f:
            first_line = f"# combination: {final_edges}\n"
            f.write(first_line)
            np.savetxt(f, optimized_vertices)

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
    options.add_argument("-opt", "--optimizer", choices=['cg', 'bfgs', 'l-bfgs-b'], default="cg", help="which optimizer to use")
    options.add_argument("-n", "--n_steps", type=int, default="25000", help="N steps before optimization stops")
    options.add_argument('-e', "--energy_method", type=str, choices=['old', 'new'], help="""old assumes that in initial geometry the
angles areclose to ideal angle and in the new one they aren't, the methods use different energy functions albeit a little similar.""")
    args = parser.parse_args()

    r_disk = 7
    L = 2 # half-distance between proteins

    main(geometry_file=args.inputfile, L=L, R=r_disk, save=args.save, conserve_membrane=args.conserve_membrane,
         output_file=args.outputfile, repulsion=args.repulsion, n_steps=args.n_steps, optimizer=args.optimizer,
         energy_method=args.energy_method)

    sys.exit(0)
