#!/usr/bin/env python3
"""
Script to update GFF file by making gene IDs (Soltu, PGSC, Sotub) unique with haplotype information
"""

import re
import sys

def extract_haplotype(chr_name):
    """Extract haplotype number from chromosome name (e.g., chr_01_3 -> 3)"""
    match = re.search(r'_(\d+)$', chr_name)
    return match.group(1) if match else None

def update_gene_ids(attributes, haplotype):
    """Update gene IDs (Soltu, PGSC, Sotub) in the attributes string with haplotype prefix"""
    if not haplotype:
        return attributes

    # Pattern to match gene IDs starting with Soltu, PGSC, or Sotub adapt for rice AGIS_
    gene_pattern = r'((?:Soltu|PGSC|Sotub|AGIS_)[^;,\s]+)'

    def replace_gene_id(match):
        original_id = match.group(1)
        return f"{original_id}.Hap{haplotype}"

    return re.sub(gene_pattern, replace_gene_id, attributes)

def process_gff_file(input_file, output_file):
    """Process GFF file and update gene IDs (Soltu, PGSC, Sotub) with haplotype information"""

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    outfile.write(line + '\n')
                    continue

                # Split GFF line into components
                parts = line.split('\t')
                if len(parts) != 9:
                    print(f"Warning: Line {line_num} doesn't have 9 columns, skipping: {line}")
                    outfile.write(line + '\n')
                    continue

                chr_name = parts[0]
                attributes = parts[8]

                # Extract haplotype from chromosome name
                haplotype = extract_haplotype(chr_name)

                if haplotype:
                    # Update gene IDs in attributes
                    updated_attributes = update_gene_ids(attributes, haplotype)
                    parts[8] = updated_attributes

                # Write updated line
                outfile.write('\t'.join(parts) + '\n')

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python update_gff_haplotypes.py <input_gff> <output_gff>")
        print("Example: python update_gff_haplotypes.py input.gff3 output.gff3")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Processing {input_file}...")
    process_gff_file(input_file, output_file)
    print(f"Updated GFF file saved as {output_file}")

if __name__ == "__main__":
    main()
