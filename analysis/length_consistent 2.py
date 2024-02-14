import sys
import ast

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

final_edges = calculate_edge_lengths("last_geometry.txt")

string_list = sys.argv[1]
initial_edges = ast.literal_eval(string_list)

print(f"initial: {initial_edges}")
print(f"final: {final_edges}")
