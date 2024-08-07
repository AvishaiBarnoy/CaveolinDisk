import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import euclidean
import sys
import re
import argparse

class GeometryOptimizer:
    def __init__(self, vertices, ideal_distances, k_edges, k_angle,
                 optimizer, n_steps=5000,
                 energy_method='new',
                 ideal_angles=[], repulsion=False, conserve_membrane=False, save=True, L=None,
                 c0=0.3):
        # TODO: decide which parameter are input and which are initialized inside

        # TODO: change vertices to not be flattend, only faltten when enter minimize
        self.vertices = vertices
        self.num_vertices = len(self.vertices)//2
        self.ideal_distances = ideal_distances

        self.energy_method = energy_method.lower() # TODO: implement in the future use choices

        # assumes initial geoemtry is close to optimal one
        # TODO: remove in future versions
        if self.energy_method == "old":
            def initiate_ideal_angles(geometry, R=7, h=2, a=0.7, depsilon=4, k=0.4e-19, kt=20e-3):
                """takes initial geometry and initializes ideal angles list"""
                # TODO: move to energy calculation for energy_method!
                xi = np.sqrt(k/kt) * 1e9 # J / nm
                f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
                geometry = np.array(geometry).reshape((self.num_vertices, 2))
                cyclic_geometry = np.vstack([geometry, geometry[0]])
                angle_lst = []
                for i, _ in enumerate(geometry):
                    if i % 2 != 0:
                        dist = euclidean(cyclic_geometry[i], cyclic_geometry[i+1])
                        ideal_angle = np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(dist/xi))
                        angle_lst.append(ideal_angle)
                return angle_lst

            # allows for external input of ideal angle
            if ideal_angles:
                self.ideal_angles = ideal_angles
            elif not ideal_angles:
                self.ideal_angles = initiate_ideal_angles(self.vertices)

        self.k_edges = k_edges
        self.k_angle = k_angle
        self.repulsion = repulsion # check if this is needed
        self.c0 = c0 # cholesterol concentration in reservoir 
        self.num_vertices = len(self.vertices)//2
        self.conserve_membrane = conserve_membrane
        self.save = save

        if not L:
            # TODO: placeholder to use this
            self.L = 2
            # self.num_vertices = len(self.vertices)//2
            self.total_membrane = self.L * self.num_vertices

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
            # R_lst = [val / 2 for val in self.ideal_distances[0::2] for _ in range(2)]

            for i in range(self.num_vertices):
                # TODO: change to % i and enumerate for better readability
                prev_point = vertices[i - 1]
                current_point = vertices[i]
                if i == self.num_vertices - 1:
                    next_point = vertices[0]
                else:
                    next_point = vertices[i + 1]

                if i % 2 == 0:
                    R_i = euclidean(next_point, current_point) / 2
                    L_i = euclidean(prev_point, current_point) / 2
                    energy += self.calc_new_energy(L=L_i, phi=current_angles[i], R=R_i)
                if i % 2 != 0:
                    R_i = euclidean(prev_point, current_point) / 2
                    L_i = euclidean(next_point, current_point) / 2
                    energy += self.calc_new_energy(L=L_i, phi=current_angles[i], R=R_i)

            for i,_ in enumerate(L_lst):
            # for i in range(self.num_vertices):
                # TODO: use this logic for constraints and edge calculation, maybe refactor
                next_vertex = (i + 1) % self.num_vertices
                protein_length = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i] # TODO: rename

                if i % 2 == 0:
                    # constraint only on proteins/disks
                    energy += self.k_edges * (protein_length - ideal_distance) ** 2

            # constraint on total membrane, k_edges should be big
            total_membrane = self.calc_total_membrane_area(self.ideal_distances)
            current_membrane = sum(L_lst[0::2]) # no L_st/2 because we need full membrane area
            energy += self.k_edges * (current_membrane - total_membrane) ** 2

        elif self.energy_method == "cholesterol":
            # iterate over length and angles
            # for each length check two angles
            L_lst = self.calc_L_lst(vertices)
            # R_lst = [val / 2 for val in self.ideal_distances[0::2] for _ in range(2)]

            for i in range(self.num_vertices):
                # TODO: change to % i and enumerate for better readability
                prev_point = vertices[i - 1]
                current_point = vertices[i]
                if i == self.num_vertices - 1:
                    next_point = vertices[0]
                else:
                    next_point = vertices[i + 1]

                if i % 2 == 0:
                    R_i = euclidean(next_point, current_point) / 2
                    L_i = euclidean(prev_point, current_point) / 2
                    # TODO: can sum switch in energy method here and not complete same block of code
                    energy += self.calc_cholesterol_energy(L=L_i, phi=current_angles[i], R=R_i)
                if i % 2 != 0:
                    R_i = euclidean(prev_point, current_point) / 2
                    L_i = euclidean(next_point, current_point) / 2
                    # TODO: can sum switch in energy method here and not complete same block of code
                    energy += self.calc_cholesterol_energy(L=L_i, phi=current_angles[i], R=R_i)

            for i,_ in enumerate(L_lst):
            # for i in range(self.num_vertices):
                # TODO: use this logic for constraints and edge calculation, maybe refactor
                next_vertex = (i + 1) % self.num_vertices
                protein_length = np.linalg.norm(vertices[next_vertex] - vertices[i])
                ideal_distance = self.ideal_distances[i] # TODO: rename

                if i % 2 == 0:
                    # constraint only on proteins/disks
                    energy += self.k_edges * (protein_length - ideal_distance) ** 2

            # constraint on total membrane, k_edges should be big
            total_membrane = self.calc_total_membrane_area(self.ideal_distances)
            current_membrane = sum(L_lst[0::2]) # no L_st/2 because we need full membrane area
            energy += self.k_edges * (current_membrane - total_membrane) ** 2

        elif self.energy_method == 'old':
            if self.conserve_membrane == False:
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
                        energy += self.elastic_energy(distance/2)

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

    def calc_cholesterol_energy(self, L, phi, zeta=1/3, c0=0.3, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
        """
        energy calculation for each L and phi
        """
        # constants
        xi = np.sqrt(k/kt) * 1e9 # J / nm
        K = np.sqrt(k*kt) * 1e-9 / 4.11e-21 # kT/nm
        e = h/a * depsilon # kT

        # in this forumlation l, r are defined as the parameter inside the hyperbolic functions and not the functions themselves 
        l = L/xi
        rho = R/xi

        # convert phi according to convention
        phi = np.pi - phi

        Ein =   K * (1/np.tanh(rho) - 0.5 * k/4.11e-21 * zeta**2 * a*c0*(1-c0) * (rho + np.sinh(rho)*np.cosh(rho))/(np.sinh(rho)**2))
        Eout =  K * (1/np.tanh(l)   - 0.5 * k/4.11e-21 * zeta**2 * a*c0*(1-c0) * (l + np.sinh(l)*np.cosh(l))/(np.sinh(l)**2))

        F = (2*Ein+Eout)/(2*Ein+2*Eout)*Eout*phi**2 - Eout/(Ein+Eout)*e*phi - e**2 / (2*Ein+2*Eout)
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

    def update_ideal_angles(self, geometry, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
        """
        updates ideal angles from geometry, since id_angle is a function of L
        """
        xi = np.sqrt(k/kt) * 1e9 # J / nm
        f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
        cyclic_geometry = np.vstack([geometry, geometry[0]])
        angle_lst = []
        for i, _ in enumerate(cyclic_geometry):
            if i % 2 != 0:
                dist = euclidean(cyclic_geometry[i], cyclic_geometry[i+1])
                ideal_angle = np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(dist/xi))
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
        sys.stdout.write(f"final angles:\n {self.calculate_angles(optimized_vertices)}\n")
        print("Finished optimizing")
        # TODO: mv priniting out of main into here, maybe use an auxillary function for printing
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
        dist = euclidean(current_point, next_point)
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

def calc_k(L, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7):
    """
    calculates spring constant for quadratic approximation for deviation angle from ideal angle
    """
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def main(geometry_file, L_i, R, output_file, save=True, conserve_membrane=False, repulsion=False,
         optimizer='cg', n_steps=5000, energy_method='new', cholesterol=0.3):
    # TODO: L_i is not used!
    k_edges = 100.0 # to keep proteins/disks rigid 

    if repulsion:
        k_angle = calc_k(L=L_i, R=R)
    elif not repulsion:
        # if elastic energy is not taken into account optimization works only on angles
        #   and the value of k_angle is not important just needs to be much smaller than k_edges
        k_angle = 1.0

    initial_geometry = np.loadtxt(args.inputfile)
    ideal_distances = get_ideal_dist(geometry_file) # TODO: rename, name misleading due to exapnsion of optimization options
    sys.stdout.write(f"Initial configuration: {ideal_distances}\n")

    sys.stdout.write(f"""
Protein radius: {r_disk} nm
Half-distance: {L_i} nm

minimization method: {optimizer}
energy calculation method: {energy_method}\n""")

    if conserve_membrane:
        sys.stdout.write("""Membrane treatment:
distance between proteins is allowed to change but total membrane is conserved\n""")
    elif not conserve_membrane:
        sys.stdout.write("""Membrane treatment:
total membrane is conserved but excesss membrane is distributed uniformly between proteins.\n""")

    optimizer = GeometryOptimizer(vertices=initial_geometry.flatten(), ideal_distances=ideal_distances,
                                  optimizer=optimizer, n_steps=n_steps, energy_method='new',
                                  k_edges=k_edges, k_angle=k_angle, repulsion=repulsion, conserve_membrane=conserve_membrane)
    optimized_vertices, _ = optimizer.optimize_geometry() # _ is final energy but is unassigned 'cause it's unused
    optimized_vertices = optimized_vertices.reshape((len(initial_geometry), 2))

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
    # TODO: change conserve_membrane so there are four scenarios:
    #           1. conserved and initial geometry membrane equals conserved membrane
    #           2. don't conserve and initial geometry membrane equals conserved membrane
    #           3. conserved and initial geoemtry membrane isn't equal conserved membrane
    #           4. don't conserve and initial geometry membrane isn't equals conserved membrane
    options.add_argument("-c", "--conserve_membrane", action="store_true", help="conserve membrane and allow distance between proteins to change")
    options.add_argument("-r", "--repulsion", action="store_true", help="minimizes with protein repulsion")
    options.add_argument("-opt", "--optimizer", choices=['cg', 'bfgs', 'l-bfgs-b'], default="cg", help="which optimizer to use, default: conjugate gradient")
    options.add_argument("-n", "--n_steps", type=int, default="25000", help="N steps before optimization stops")
    options.add_argument('-e', "--energy_method", type=str, choices=['old', 'new', 'cholesterol'], default="new", help="""old assumes that in initial geometry the angles areclose to ideal angle and in the new one they aren't, the methods use different energy functions albeit a little similar.
default: new energy method""")
    options.add_argument("-C", "--cholesterol", type=float, default=0.3, help="cholesterol concentration in reservoir, default: 0.3")
    options.add_argument("-L", type=float, help="half-distance")
    args = parser.parse_args()

    r_disk = 7
    L = 2 # half-distance between proteins
    # TODO: calculate total membrane from initial L

    main(geometry_file=args.inputfile, L_i=L, R=r_disk, save=args.save, conserve_membrane=args.conserve_membrane,
         cholesterol=args.cholesterol, output_file=args.outputfile, repulsion=args.repulsion, n_steps=args.n_steps, optimizer=args.optimizer,
         energy_method=args.energy_method)

    sys.exit(0)
