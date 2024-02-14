import sys
import ast
import numpy as np
import re

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

geom_file = sys.argv[1]
with open(geom_file, 'r') as f:
    first_line = f.readline().strip()
    pattern = r"\[(.*?)\]"
    match = re.search(pattern, first_line)
initial_edges = np.array(ast.literal_eval('['+match[1]+']'))

# print(f"initial: {initial_edges}")
# print(f"final: {final_edgese")
# print(f"diff: {abs(final_edges - initial_edges)}")
diff = np.count_nonzero(final_edges - initial_edges)
# print(diff)
if diff == 0:
    sys.stdout.write(f"{sys.argv[1]} edges consistent\n")
else: 
    sys.stdout.write(f"{sys.argv[1]} edges NOT consistent\n")
