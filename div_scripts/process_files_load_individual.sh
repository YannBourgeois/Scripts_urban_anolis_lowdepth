#!/bin/bash

# Output file
output_file="mutation_load_individual.txt"

# Clear the output file if it already exists
> $output_file

# Loop through each file matching the pattern
for file in Results_for_individual_load/*.ml; do
    # Extract the filename without the path
    filename=$(basename "$file")
    
    # Extract the base name and the effect
    base_name=$(echo "$filename" | awk -F'_mdup' '{print $1 "_mdup"}')
    effect=$(echo "$filename" | awk -F'.bam_' '{print $2}' | awk -F'.ml' '{print $1}')
    
    # Read the content of the file and append it to the output
    awk -v base_name="$base_name" -v effect="$effect" '{
        print base_name "\t" effect "\t" $0
    }' "$file" >> $output_file
done
