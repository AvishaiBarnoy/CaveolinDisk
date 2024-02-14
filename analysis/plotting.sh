#!/bin/bash

# Define the Bash function
plot_geometries() {
    # Check if a filename is provided
    if [ -z "$1" ]; then
        echo "Usage: plot_geometries <filename>"
        return 1
    fi
    
    line=`head -1 $1`
    
    count=$(echo "$line" | grep -o '\[.*\]' | awk -F ',' '{print NF}')
    echo "Number of elements between brackets: $count"

    # Extract first and last geometries
    head -n "$((count + 2))" "$1" | tail -n "$count" > first_geometry.txt
    tail -n "$count" "$1" > last_geometry.txt

    # Call Python script to plot the geometries
    python plot_geometries.py first_geometry.txt last_geometry.txt

    # Clean up temporary files
   # rm first_geometry.txt last_geometry.txt
}

# Call the function with the provided filename as argument
plot_geometries "$1"
