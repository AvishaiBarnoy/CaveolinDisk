import os
import numpy as np
import argparse
from scipy.spatial.distance import euclidean

def calculate_angles(points):
    """
    calculates angles from a set of [x,y] coordiantes
    """
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
    return np.array(angles)

def angles_energy(file_path, k=1, ideal_angle=2.930):
    """
    calculates energy of angles deviation from ideal angle
    """
    # Your calculation logic here, this function should return a float
    geometry = np.loadtxt(file_path)
    angles = calculate_angles(geometry)
    last_value = 0.5 * k * np.sum(np.power(angles - ideal_angle, 2))
    return last_value

def calculate_new_energy(file_path):
    points = np.loadtxt(file_path)
    geometry_energy = 0
    for i in range(len(points)):
        prev_point = points[i - 1]
        current_point = points[i]
        if i == len(points) - 1:
            next_point = points[0]
        else:
            next_point = points[i + 1]

        # get phi_i
        vector1 = prev_point - current_point
        vector2 = next_point - current_point
        dot_product = np.dot(vector1, vector2)
        magnitude_product = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        phi_i= np.arccos(dot_product / magnitude_product)

        # get L_i and R_i
        if i % 2 == 0:
            R_i = euclidean(next_point, current_point) / 2
            L_i = euclidean(prev_point, current_point) / 2
        elif i % 2 != 0:
            R_i = euclidean(prev_point, current_point) / 2
            L_i = euclidean(next_point, current_point) / 2

        geometry_energy += new_energy(L=L_i, phi=phi_i, R=R_i)
    return geometry_energy

def new_energy(L, phi, R, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
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

def find_lowest_value_directory(root_dir, k, ideal_angle, energy_method='old'):
    """
    loop through all directories to calculate angle from structures
    returns lowest energy for each structure
    """
    lowest_value = float('inf')
    lowest_dir = None
    structures_energy = []
    for dir_name in os.listdir(root_dir):
        dir_path = os.path.join(root_dir, dir_name)
        if os.path.isdir(dir_path):
            optimized_geometry = 'geom_opt.txt'
            last_geometry_file = os.path.join(dir_path, optimized_geometry)
            if os.path.exists(last_geometry_file):
                if energy_method == 'old':
                    value = angles_energy(last_geometry_file, k, ideal_angle)
                elif energy_method == 'new':
                    value = calculate_new_energy(last_geometry_file)
                structures_energy.append(value)
                if value < lowest_value:
                    lowest_value = value
                    lowest_dir = dir_name
    return lowest_value, lowest_dir


def calc_ideal_angle(L, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7, depsilon=4):
    '''
    f_param = Delta_epsilon / sqrt(k*k_t) * h / a
    D_epsilon - energy diff of lipids on-top of protein and in membrane ~ 0.3
    h - monolayer thickness ~ 2 nm, a - lipid length ~ 0.7 nm
    k - monolayer rigidiy ~ 1e-9 J, k_t - tilt modulus ~ 30 mJ/m
    '''
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
    return np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="for each ΔN finds the lowest optimized combination",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument('-e', "--energy_method", type=str, choices=['old', 'new'], default="new", help="""old assumes that in
initial geometry the angles areclose to ideal angle and in the new one they aren't, the methods use different energy functions
albeit a little similar.""")

    args = parser.parse_args()
    k = 1
    L = 5
    ideal_angle = calc_ideal_angle(L=L, R=7)

    root_directory = os.getcwd()
    lowest_value, lowest_dir = find_lowest_value_directory(root_directory, k=k, ideal_angle=ideal_angle, energy_method=args.energy_method)

    # can probably remove after automation works
    if args.energy_method == 'old':
        print(f"Ideal angle: {ideal_angle}, with kappa_tilde: {k}")
    print(f"Lowest value: {lowest_value}, found in directory: {lowest_dir}")

    # copy geometry to analysis folder  
    os.system(f"cp {lowest_dir}/geom_opt.txt ../analysis/{root_directory[-4:]}.txt")

    # append energy to energy.py
    os.system(f"echo {lowest_value} >> ../analysis/e_values.txt")
