#!/usr/bin/env python3
"""
Promoter Region Extractor

This script extracts promoter regions from genomic data based on GFF annotations.
It can optionally remove overlaps with coding sequences.

Usage:
    python promoter_extractor.py --gff <gff_file> --fasta <fasta_file> \
                                 --output <output_file> --length <promoter_length> \
                                 [--tss-downstream <tss_downstream_length>] \
                                 [--remove-overlaps] \
                                 [--discard_short_promoters]
"""

import argparse
import os
import sys
import pyranges as pr
import pandas as pd
import pyfaidx
from pathlib import Path


def extract_promoters(gff_file, fasta_file, output_file, promoter_length=2000, tss_downstream_length=0, downstream_length = 0, remove_overlaps=False, discard_short_promoters=False):
    """
    Extract promoter regions from GFF file based on CDS features.

    Args:
        gff_file (str): Path to the GFF file
        fasta_file (str): Path to the FASTA file
        output_file (str): Path to the output GFF file
        promoter_length (int): Length of the promoter region to extract
        tss_downstream_length (int): Length to extend the promoter region into the gene body after the TSS
        downstream_length (int): Length of the downstream region to extract
        remove_overlaps (bool): Whether to remove overlaps with coding sequences
        discard_short_promoters (bool): Whether to discard promoters shorter than the specified length due to overlaps
    """
    print(f"Reading GFF file: {gff_file}")
    ann = pr.read_gff3(gff_file)
    print(ann)
    # Select CDS features
    print(f"Extracting CDS/exon features")
    selector = (ann.Feature == 'exon')
    cds = ann[selector]

   

    # Build promoter regions directly in genomic coordinates around the TSS.
    # This avoids transcript-spliced coordinates when adding --tss-downstream.
    promoter_total_length = promoter_length + tss_downstream_length
    print(f"Extending promoter sequence through TSS by {tss_downstream_length}bp (total promoter length: {promoter_total_length}bp)")
    print(f"Building promoter intervals in genomic coordinate space")

    exon_df = cds.df.copy()
    promoter_records = []

    for parent, grp in exon_df.groupby("Parent", sort=False):
        first = grp.iloc[0].copy()
        strand = first["Strand"]

        if strand == '+':
            # TSS is the left-most exon boundary on the positive strand.
            tss = int(grp["Start"].min())
            start = tss - promoter_length
            end = tss + tss_downstream_length
        elif strand == '-':
            # TSS is the right-most exon boundary on the negative strand.
            tss = int(grp["End"].max())
            start = tss - tss_downstream_length
            end = tss + promoter_length
        else:
            print(f"Warning: Unsupported strand '{strand}' for Parent '{parent}'. Skipping.")
            continue

        first["Start"] = int(start)
        first["End"] = int(end)
        first["Feature"] = "exon"
        promoter_records.append(first)

    if not promoter_records:
        raise ValueError("No promoter intervals could be built from exon features.")

    prom = pr.PyRanges(pd.DataFrame(promoter_records))

    print(f"Extending feautres by downstream length of {downstream_length}bp")
    if downstream_length > 0:
        d = cds.extend({'3': downstream_length}, group_by='Parent')
        d.Feature = 'exon'
        down = d.spliced_subsequence(-(downstream_length+1), -1, "Parent")
    else:
        down = pr.PyRanges()
    # Correct bounds to chromosome size
    print(f"Correcting bounds to chromosome size")
    pyf = pyfaidx.Fasta(fasta_file)
    prom = pr.genomicfeatures.genome_bounds(prom, chromsizes=pyf, clip=True)
    if len(down) > 0:
        down = pr.genomicfeatures.genome_bounds(down, chromsizes=pyf, clip=True)
        down_df = down.df.copy()
    else:
        down_df = pd.DataFrame()

    if remove_overlaps or discard_short_promoters:
        print(f"Removing overlaps with coding sequences")
        # Check the overlap of promoters and CDS
        prom_in_cds = prom.intersect(cds, strandedness=False)

        # If there are overlaps, adjust the promoters
        if not prom_in_cds.df.empty:
            prom_in_cds_merged = prom_in_cds.merge(by='Parent', slack=10000)
            prom_in_cds_merged_df = prom_in_cds_merged.df

            # For positive strand promoters, set the start to be the end of the overlap
            mask_positive = prom_in_cds_merged_df['Strand'] == '+'
            prom_in_cds_merged_df.loc[mask_positive, 'Start'] = prom_in_cds_merged_df.loc[mask_positive, 'End']

            # For negative strand promoters, set the end to be the start of the overlap
            mask_negative = prom_in_cds_merged_df['Strand'] == '-'
            prom_in_cds_merged_df.loc[mask_negative, 'End'] = prom_in_cds_merged_df.loc[mask_negative, 'Start']

            # Get the positive strand features
            prom_in_cds_merged_pos = prom_in_cds_merged_df[prom_in_cds_merged_df.Strand == '+']
            prom_overlap_start_pos = prom_in_cds_merged_pos[["Start", "Parent"]]

            # Get the negative strand features
            prom_in_cds_merged_neg = prom_in_cds_merged_df[prom_in_cds_merged_df.Strand == '-']
            prom_overlap_stop_neg = prom_in_cds_merged_neg[["End", "Parent"]]

            # Update the promoter dataframe
            prom_df = prom.df
            prom_df.set_index('Parent', inplace=True)
            prom_overlap_start_pos.set_index('Parent', inplace=True)
            prom_overlap_stop_neg.set_index('Parent', inplace=True)

            prom_df.update(prom_overlap_start_pos, overwrite=True)
            prom_df.update(prom_overlap_stop_neg, overwrite=True)

            # Reset index and prepare for output
            prom_df["Parent"] = prom_df.index
            prom_df.reset_index(drop=True, inplace=True)
        else:
            prom_df = prom.df
    else:
        prom_df = prom.df

    if discard_short_promoters:
        print(f"Discarding promoters shorter than {promoter_total_length}bp")
        # Filter out promoters shorter than the specified length
        prom_df = prom_df[(prom_df['End'] - prom_df['Start']) >= promoter_total_length]


    # Add Name and ID fields for the promoters
    
    print(prom_df.columns)
    #prom_df.to_csv('promoter_df.csv', index=False)  # Debugging line to check the DataFrame
    if 'gene_id' in prom_df.columns:
        prom_df["ID"] = "promoter_" + prom_df["gene_id"]
        prom_df["Parent"] = prom_df["gene_id"]
        prom_df["Name"] = "promoter_" + prom_df["gene_id"]
        
        if not down_df.empty:
            down_df["ID"] = "downstream_" + down_df["gene_id"]
            down_df["Parent"] = down_df["gene_id"]
            down_df["Name"] = "downstream_" + down_df["gene_id"]
        else:
            print("Warning: Downstream dataframe is empty.")
        # Drop Parent column to avoid problems with relationships
        #prom_df.drop(columns=['Parent'], inplace=True)

    else:
        print("Warning: No 'Parent' column found in the dataframe.")

     # concat prom and down
    prom_df = pd.concat([prom_df, down_df], ignore_index=True)
    # Sort and prepare for output
    prom_df.sort_values(by=['Chromosome', 'Start'], inplace=True)
    # Create PyRanges object and set feature type
    prom_updated = pr.PyRanges(prom_df)
    prom_updated.Feature = 'exon'

    # Write to output file
    print(f"Writing output to: {output_file}")
    prom_updated.to_gff3(output_file)

    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Total promoters extracted: {len(prom_df)}")
    if remove_overlaps:
        print(f"  Promoters adjusted for CDS overlaps: {len(prom_in_cds_merged_df) if 'prom_in_cds_merged_df' in locals() else 0}")

    return prom_updated


def main():
    parser = argparse.ArgumentParser(description='Extract promoter regions from genomic data.')
    parser.add_argument('--gff', required=True, help='Input GFF file')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output GFF file')
    parser.add_argument('--length', type=int, default=3000, help='Promoter and downstreamlength (default: 3000bp)')
    parser.add_argument('--tss-downstream', type=int, default=0, help='Extend the promoter region downstream of the TSS into the gene body (default: 0bp)')
    parser.add_argument('--remove-overlaps', action='store_true', help='Remove overlaps with coding sequences')
    parser.add_argument('--discard_short_promoters', action='store_true', help='Discard promoters shorter than the specified length due to overlaps')

    args = parser.parse_args()

    # Check if input files exist
    if not Path(args.gff).is_file():
        sys.exit(f"Error: GFF file '{args.gff}' not found.")
    if not Path(args.fasta).is_file():
        sys.exit(f"Error: FASTA file '{args.fasta}' not found.")

    # Check if output directory exists
    # output_dir = Path(args.output)
    # print(output_dir)
    # if not output_dir.exists():
    #     sys.exit(f"Error: Output directory '{output_dir}' does not exist.")

    # Extract promoters
    try:
        extract_promoters(
            args.gff,
            args.fasta,
            args.output,
            args.length,
            args.tss_downstream,
            0,  
            args.remove_overlaps,
            args.discard_short_promoters
        )
        print("Promoter extraction completed successfully.")
    except Exception as e:
        sys.exit(f"Error: {str(e)}")


if __name__ == "__main__":
    main()
