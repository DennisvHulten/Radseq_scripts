#!/usr/bin/env python3

import sys

def load_clonal_groups(clonal_groups_file):
    clonal_groups = []
    with open(clonal_groups_file, 'r') as f:
        for idx, line in enumerate(f, start=1):  # Start numbering from 1
            group, members = line.strip().split(":")  # Strip to remove newline issues
            members = [m.strip().split("(")[0].strip() for m in members.split(",")]  # Remove anything after "("
            group_number = idx
            clonal_groups.append((group.strip(), members, group_number))
    return clonal_groups

def process_popfile(popfile, clonal_groups):
    output = []
    
    with open(popfile, 'r') as f:
        for line in f:
            sample_id, pop = line.strip().split("\t")
            clonal_group_number = None
            
            # Check if sample_id is in any of the clonal groups
            for group, members, group_number in clonal_groups:
                if sample_id in members:
                    clonal_group_number = group_number
                    break
            
            # If sample_id is found in a group, add the clonal group number
            if clonal_group_number:
                output.append(f"{sample_id}, {pop}, {clonal_group_number}")
            else:
                output.append(f"{sample_id}, {pop}, None")  # In case no group matches
    
    return output

def main(popfile, clonal_groups_file, output_file):
    clonal_groups = load_clonal_groups(clonal_groups_file)
    updated_popfile = process_popfile(popfile, clonal_groups)
    
    # Write the updated information to the output file
    with open(output_file, 'w') as f:
        for line in updated_popfile:
            f.write(line + "\n")

if __name__ == "__main__":
    # Ensure the right number of arguments are passed
    if len(sys.argv) != 4:
        print("Usage: python script.py <popfile> <clonal_groups_file> <output_file>")
        sys.exit(1)
    
    # Run the script with the provided arguments
    popfile = sys.argv[1]
    clonal_groups_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(popfile, clonal_groups_file, output_file)
