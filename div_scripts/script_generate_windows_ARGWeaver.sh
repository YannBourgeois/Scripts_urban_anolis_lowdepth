#!/bin/bash

# Define the start and end coordinates
start=121500000
end=123500000
chr="chr4"
interval=50000

# Open a file to write the results
output_file="array_arg.txt"
> "$output_file"

# Loop over the range in 50,000 bp intervals
for ((i=start; i<end; i+=interval)); do
    from=$i
    to=$((i+interval))
    if [ $to -gt $end ]; then
        to=$end
    fi
    echo -e "$chr:$from-$to" >> "$output_file"
done

echo "Intervals have been generated and saved to $output_file"
