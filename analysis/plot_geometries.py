# plot_geometries.py

import matplotlib.pyplot as plt

def read_geometry(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        geometry = []
        for line in lines:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            x = float(parts[0])
            y = float(parts[1])
            geometry.append((x, y))
        return geometry

def plot_geometries(first_geometry_file, last_geometry_file):
    first_geometry = read_geometry(first_geometry_file)
    last_geometry = read_geometry(last_geometry_file)
    
    # add first point to end for closed form plotting
    first_geometry.append(first_geometry[0])
    last_geometry.append(last_geometry[0])

    # Plotting the geometries
    plt.figure(figsize=(8, 6))
    plt.plot(*zip(*first_geometry), label='First Geometry')
    plt.plot(*zip(*last_geometry), label='Last Geometry')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('First and Last Geometries')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python plot_geometries.py <first_geometry_file> <last_geometry_file>")
        sys.exit(1)
    first_geometry_file = sys.argv[1]
    last_geometry_file = sys.argv[2]
    plot_geometries(first_geometry_file, last_geometry_file)

