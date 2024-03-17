import os
import numpy as np
import re
import argparse

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

def angles_energy(file_path, ideal_angle, L, k, R=7):
    """
    calculates energy of angles deviation from ideal angle
    """
    geometry = np.loadtxt(file_path)
    angles = calculate_angles(geometry)
    last_value = 0.5 * k * np.sum(np.power(angles - ideal_angle, 2))
    return last_value

def calc_k(L, R=7, xi=2):
    h = 2 # nm
    a = 0.7 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    K_const = np.sqrt(k*kt/1e18) * (2/np.tanh(R/xi) + 1/np.tanh(L/xi))*1/np.tanh(L/xi) / (1/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return K_const / 4.11e-21

def repulsion_energy(file_path, L, R=7):
    h = 2 # nm
    a = 0.7 # nm
    xi = 2 # nm
    k = 0.8e-19 # J
    kt = 30e-3 # N/m
    depsilon = 4 # kT/nm

    A = (h/a)**2 * depsilon/np.sqrt(k*kt/1e18) * 4.11e-21
    F = - A / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))
    return F

def total_energy(file_path, ideal_angle, k, L):
    e_angles = angles_energy(file_path, ideal_angle, k=k, L=L)
    e_repulson = repulsion_energy(file_path, L)
    return e_angles + e_repulson

def find_lowest_value_directory(root_dir, L, k, ideal_angle):
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
                value = total_energy(last_geometry_file, ideal_angle=ideal_angle, L=L, k=k)
                structures_energy.append(value)
                if value < lowest_value:
                    lowest_value = value
                    lowest_dir = dir_name
    return lowest_value, lowest_dir


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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="lowest_dn",
        description="extract lowest energy structure for each ΔN",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-o", "--outputfile", help="path to output file", default="e_values.txt")
    parser.add_argument("-r", "--repulsion", action="store_true", help="""takes protein's repulsion excilictly when checking for structure with minimum energy,
CATUION: returned energy includes repulsion so make sure you use appropriate plotting scripts not to double counting energy""")
    args = parser.parse_args()

    L = 1
    if args.repulsion:
        k = calc_k(L, R=7, xi=2)
    elif not args.repulsion:
        k = 1
    ideal_angle = calc_ideal_angle(L, xi=2, R=7)

    root_directory = os.getcwd()
    lowest_value, lowest_dir = find_lowest_value_directory(root_dir=root_directory, ideal_angle=ideal_angle, L=L, k=k)

    # can probably remove after automation works
    print(f"Ideal angle: {ideal_angle}, k: {k}")
    print(f"Lowest value: {lowest_value}, found in directory: {lowest_dir}")

    # copy geometry to analysis folder  
    os.system(f"cp {lowest_dir}/geom_opt.txt ../analysis/{root_directory[-4:]}.txt")

    # append energy to energy.py
    os.system(f"echo {lowest_value} >> ../analysis/{args.outputfile}")
