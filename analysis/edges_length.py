import numpy as np
import sys

def calculate_edge_lengths(geometry_file):
    with open(geometry_file, 'r') as file:
        lines = file.readlines()

    # Extract coordinates from the file
    coordinates = []
    for line in lines:
        if not line.startswith('#'):
            parts = line.strip().split()
            x = float(parts[0])
            y = float(parts[1])
            coordinates.append([x, y])
    coordinates.append(coordinates[0])
    # Convert coordinates to NumPy array for easier manipulation
    coordinates = np.array(coordinates)

    # Calculate the differences between consecutive coordinates
    differences = np.diff(coordinates, axis=0)

    # Calculate the Euclidean distances (edge lengths)
    edge_lengths = np.linalg.norm(differences, axis=1)

    return np.round(edge_lengths)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python edges_length.py <first_geometry> <last_geometry>")
        sys.exit(1)

    geometry_file = sys.argv[1]
    edge_lengths = calculate_edge_lengths(geometry_file)
    print("Edge first lengths:", edge_lengths)
    
    geometry_file = sys.argv[2] 
    edge_lengths = calculate_edge_lengths(geometry_file)
    print("Edge last lengths:", edge_lengths)
