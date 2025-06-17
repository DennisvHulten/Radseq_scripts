import argparse
import numpy as np
import os

def assign_population(file_path, outfile_path):
    groups = {}

    # First pass: Read the file and determine group sizes
    with open(file_path, 'r') as file:
        for line in file:
            cols = line.strip().split(',')
            sample_id = cols[0]
            values = np.array(cols[2:], dtype=float)
            loc = cols[1]
            pop = np.where(values == np.max(values))[0][0] + 1  # Assign population
            if pop not in groups:
                groups[pop] = set()
            groups[pop].add(sample_id)

    # Determine major, minor, and other clades
    sorted_clades = sorted(groups, key=lambda k: len(groups[k]), reverse=True)
    clade_labels = {sorted_clades[0]: "Major_Clade"}  # Assign major clade
    if len(sorted_clades) > 1:
        clade_labels[sorted_clades[1]] = "Minor_Clade"  # Assign minor clade
    for i, clade in enumerate(sorted_clades[2:], start=3):  # Assign all others as "Clade_X"
        clade_labels[clade] = f"Clade_{i}"

    # Second pass: Write results with clade assignments
    with open(file_path, 'r') as file, open(outfile_path, 'w') as outfile:
        for line in file:
            cols = line.strip().split(',')
            sample_id = cols[0]
            values = np.array(cols[2:], dtype=float)
            loc = cols[1]
            pop = np.where(values == np.max(values))[0][0] + 1  # Assign population
            
            # Assign clade based on sorted ranking
            clade = clade_labels.get(pop, f"Clade_{pop}")  # Default to "Clade_X" if not found
            outfile.write(f"{sample_id}\t{loc}\tpop_{pop}\t{clade}\n")

def compare_population_call(outfile_paths):
    sample_groups = {}
    inconsistent_samples = {}
    missing_samples = {}

    for file_path in outfile_paths:
        with open(file_path, 'r') as file:
            for line in file:
                sample_id, loc, population, clade = line.strip().split()

                if sample_id not in sample_groups:
                    sample_groups[sample_id] = {}

                sample_groups[sample_id][file_path] = clade  # Use clade directly from file

        print(f"Processed file: {file_path}")

    # Check for inconsistencies
    inconsistent_samples = {
        sample_id: file_groups
        for sample_id, file_groups in sample_groups.items()
        if len(set(file_groups.values())) > 1  # Sample assigned to different clades across files
    }

    # Check for missing samples
    all_files = set(outfile_paths)
    missing_samples = {
        sample_id: list(all_files - set(file_groups.keys()))
        for sample_id, file_groups in sample_groups.items()
        if len(file_groups) < len(outfile_paths)  # Sample missing in some files
    }

    # Output results
    with open('inconsistent_samples_popcall_d1.txt', 'w') as outfile:
        outfile.write("Inconsistent Samples:\n")
        for sample_id, file_groups in inconsistent_samples.items():
            outfile.write(f"Sample {sample_id}: {file_groups}\n")

    with open('missing_samples_popcall_d1.txt', 'w') as outfile:
        outfile.write("Missing Samples (Not present in every file):\n")
        for sample_id, files in missing_samples.items():
            outfile.write(f"Sample {sample_id} missing in: {files}\n")

def main():
    parser = argparse.ArgumentParser(description="Compare population assignments from multiple CSV files.")
    parser.add_argument('csv_files', nargs='+', help="List of CSV files to process")
    args = parser.parse_args()
    
    output_files = []
    for csv_file in args.csv_files:
        outfile = os.path.splitext(csv_file)[0] + "_pop_assignments.txt"
        assign_population(csv_file, outfile)
        output_files.append(outfile)
    
    compare_population_call(output_files)
    print("Comparison complete. Check inconsistent_samples_popcall_d1.txt and missing_samples_popcall_d1.txt for results.")

if __name__ == "__main__":
    main()
