#!/bin/zsh

# Check if the number of arguments is correct
if [ $# -ne 1 ]; then
    echo "Usage: $0 <new_L_value>"
    exit 1
fi

# Assign the new value for L from the command-line argument
new_L_value=$1

# TODO: if all values equal to new_L_value skip and don't sed
#           grep "L\ =" *.py | awk '{print $4}'
# Iterate over each Python file containing "L = "
for file in $(grep -l "L = " *.py); do
    # Replace the value of L with the new value using sed
    # This version should work on both macOS and Linux
    # sed -i '' -e "s/L = [0-9]\+/L = $new_L_value/" "$file"
    sed -i -e "s/L\ \= [0-9]/L\ \= $new_L_value/g" $file 
    echo "Updated $file"
done

sed -i -e "s/L\ \= [0-9]/L\ \= $new_L_value/g" $file 
