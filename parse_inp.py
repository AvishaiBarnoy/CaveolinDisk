import sys
import argparse
# TODO: add default input filename
# TODO: after reading inputfile should adjust analysis scripts accordingly
# TODO: replace sys with argparse

def read_inputfile():
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python read_inputfile.py <input_file>")
        sys.exit(1)  # Exit the program with an error code

    # Extract the input file name from the command-line arguments
    input_file = sys.argv[1]
    try:
        # initialize parameters
        params = {"r_disk": None, "L": None}

        with open(input_file, 'r') as file:
            # Perform operations on the input file
            print(f"Reading input from file: {input_file}")

            # Process each line
            for line in file:
                # Strip whitespace characters from the beginning and end of the line
                line = line.strip()

                # Ignore empty lines and comment lines
                if not line or line.startswith('#'):
                    continue

                # Process non-empty, non-comment lines
                # converts input into appropriate types, raises errors if fails
                line = line.split('=')
                if line[0].upper() == 'L':
                    try:
                        params['L'] = float(line[1])
                    except ValueError:
                        print("Error: wrong type for half-distance between disks (L) value")
                        print(f"'{line[1]}' should be an float or an integer")
                        sys.exit(1)

                if line[0].lower() == "r":
                    try:
                        params['r_disk'] = float(line[1])
                    except ValueError:
                        print("Error: wrong type for caveolin radius (R) value")
                        print(f"'{line[1]}' should be an float or an integer")
                        sys.exit(1)

        if None in params.values():
            print(f"Error: missing parameter values")


        # print(f"parameters {params}")
        return params

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)  # Exit the program with an error code

if __name__ == "__main__":
    inps = read_inputfile()
    print(f"parameters {inps}")
