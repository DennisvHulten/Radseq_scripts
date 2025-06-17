#!/usr/bin/env python3

import sys

def process_line(line):
    # Extract the part before the colon
    prefix = line.split(":")[0]
    
    # Split the prefix by "-" and extract the first 3 characters of each part
    parts = prefix.split("-")
    shortened_parts = [part[:3] for part in parts]
    
    # Check if the first three characters are not the same
    if len(shortened_parts) > 1 and shortened_parts[0] != shortened_parts[1]:
        return line
    return None

def main(file_path):
    # Open and read the input file
    with open(file_path, 'r') as f:
        for line in f:
            result = process_line(line.strip())
            if result:
                print(result)

if __name__ == "__main__":
    # Ensure a file argument is passed
    if len(sys.argv) != 2:
        print("Usage: python filter_lines.py <file>")
        sys.exit(1)
    
    # Run the script with the provided file
    main(sys.argv[1])
