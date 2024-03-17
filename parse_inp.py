import argparse
import sys
# TODO: after reading inputfile should adjust scripts accordingly

def read_inputfile(inputfile):
    try:
        # initialize parameters
        params = {"r_disk": None, "L": None, "excess_membrane": 0}

        with open(inputfile, 'r') as file:
            # Perform operations on the input file
            print(f"Reading input from file: {inputfile}")

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

                if line[0].lower() == "excess_membrane":
                    try:
                        params['excess_membrane'] = int(line[1])
                        if params['excess_membrane'] == 0:
                            # no redistribution os excess membrane
                            print("adjust scripts to not conserve membrane")
                            print("don't edit any file\n")
                        elif params['excess_membrane'] == 1:
                            # uniform distribution between membranes
                            print("adjust scripts to uniformly distribute excess membrane")
                            print("edit in generate_geom.py\n")
                        elif params['excess_membrane'] == 2:
                            # initially distribute uniformly but allow to change during
                            print("adjust scripts to uniformly distribute excess membrane")
                            print("then allow distance between proteins to change")
                            print("edit generate_geom.py and run_simulation.sh for optimzie.py cli options\n")
                    except ValueError:
                        print("Error: wrong type for optin for dealing with excess membrane")
                        print(f"'{line[1]}' should be an integer")
                        sys.exit(1)

        if None in params.values():
            print(f"Error: missing parameter values")

        # print(f"parameters {params}")
        sys.stdout.write("Input file read correctly")
        return params

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)  # Exit the program with an error code

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="parse_inp",
        description="script for input file parsing",
        epilog="Not all those who physics are lost"
    )
    parser.add_argument("-i", "--inputfile", help="path to input file with simulation instructions", default="inputfile.txt")
    args = parser.parse_args()

    inp_params = read_inputfile(inputfile=args.inputfile)

    sys.stdout.write(f"Final parameters {inp_params}\n")
    sys.exit(0)
