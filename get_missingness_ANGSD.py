#!/usr/bin/env python3
import sys

def load_sample_names(name_file):
    with open(name_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def calculate_missingness(geno_file, names):
    with open(geno_file, 'r') as f:
        total = []
        missing = []
        for line in f:
            fields = line.strip().split()
            genos = fields[4:]  # skip chr, pos, major, minor
            if not total:
                total = [0] * len(genos)
                missing = [0] * len(genos)
            for i, g in enumerate(genos):
                total[i] += 1
                if g == '-1':
                    missing[i] += 1

    for i in range(len(total)):
        name = names[i] if i < len(names) else f"Individual_{i+1}"
        miss_pct = 100 * missing[i] / total[i]
        print(f"{name}\t{miss_pct:.2f}\t({missing[i]}/{total[i]})")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_missingness_named.py <geno_file.txt> <sample_names.txt>")
        sys.exit(1)

    geno_path = sys.argv[1]
    name_path = sys.argv[2]

    sample_names = load_sample_names(name_path)
    calculate_missingness(geno_path, sample_names)