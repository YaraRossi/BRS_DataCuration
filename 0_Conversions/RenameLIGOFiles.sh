#!/bin/bash

# Check if the user provided a directory path
if [[ -z "$1" ]]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

# Assign the directory path from the command-line argument
directory="$1"

# Check if the provided path is a valid directory
if [[ ! -d "$directory" ]]; then
    echo "Error: '$directory' is not a valid directory."
    exit 1
fi

# Loop through all CSV files in the specified directory
for file in "$directory"/*.csv; do
    # Ensure we don't process non-existing files (e.g., if no CSV files are present)
    [[ -e "$file" ]] || continue

    # Extract the filename without the path
    filename=$(basename "$file")

    # Split the filename into segments using underscore as delimiter
    IFS='_' read -r segment1 segment2 segment3 segment4 segment5 segment6 rest <<< "$filename"

    # Determine Station based on segment1
    if [[ "$segment1" == "H1" ]]; then
        Station="HAN"
    elif [[ "$segment1" == "L1" ]]; then
        Station="LIV"
    else
        continue  # Skip files that don't match expected segment1 values
    fi

    # Set Network value
    Network="LI"

    # Check if segment2 starts with "PEM"
    if [[ "$segment2" == PEM* ]]; then
        Location="VA"
        Channel="DH${segment6}"
    else
        # Determine Location based on segment4
        if [[ ${#segment4} -eq 4 ]]; then
            Location="${segment4:0:1}${segment4:3:1}"
        else
            Location=""
        fi

        # Determine Channel based on segment5
        if [[ ${#segment5} -eq 2 ]]; then
            Channel="DJ${segment5:1:1}"
        else
            Channel="DH${segment5:0:1}"
        fi
    fi

    # Construct new filename
    new_filename="${Network}_${Station}_${Location}_${Channel}.csv"

    # Rename the file
    mv "$file" $1/"$new_filename"

    echo "Renamed: $file -> $directory/$new_filename"

done
