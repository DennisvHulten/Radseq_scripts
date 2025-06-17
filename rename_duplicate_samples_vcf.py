import sys

def rename_vcf_samples(vcf_file, output_file):
    sample_counts = {}
    with open(vcf_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#CHROM'):  # The header line with sample names
                parts = line.strip().split('\t')
                header = parts[:9]  # The first 9 columns are fixed VCF columns
                samples = parts[9:]  # Sample names
                renamed_samples = []
                for sample in samples:
                    if sample not in sample_counts:
                        sample_counts[sample] = 1
                        renamed_samples.append(sample)
                    else:
                        sample_counts[sample] += 1
                        renamed_samples.append(f"{sample}_x")
                outfile.write('\t'.join(header + renamed_samples) + '\n')
            else:
                outfile.write(line)  # Copy the rest of the file as is

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rename_vcf_samples.py <input.vcf> <output.vcf>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    rename_vcf_samples(input_vcf, output_vcf)
