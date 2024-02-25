import os
import numpy as np
import re

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
    # print(last_value)
    return last_value

def find_lowest_value_directory(root_dir, k, ideal_angle):
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
                value = angles_energy(last_geometry_file, k, ideal_angle)
                structures_energy.append(value)
                if value < lowest_value:
                    lowest_value = value
                    lowest_dir = dir_name
    return lowest_value, lowest_dir

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
    k = 1
    L = 2
    ideal_angle = np.pi - calc_ideal_angle(L, xi=2, R=7) 
    
    root_directory = os.getcwd()
    lowest_value, lowest_dir = find_lowest_value_directory(root_directory, k=k, ideal_angle=ideal_angle)
    
    # can probably remove after automation works
    print(f"Ideal angle: {ideal_angle}, with kappa_tilde: {k}")
    print(f"Lowest value: {lowest_value}, found in directory: {lowest_dir}")
    
    # copy geometry to analysis folder  
    os.system(f"cp {lowest_dir}/geom_opt.txt ../analysis/{root_directory[-4:]}.txt")

    # append energy to energy.py
    os.system(f"echo {lowest_value} >> ../analysis/energies.txt") 
