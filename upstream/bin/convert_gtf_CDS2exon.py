#!/usr/bin/env python3
"""
GTF/GFF CDS to Exon Converter

This script converts CDS features in a GTF/GFF file to exon features with optional
extension on both 5' and 3' ends. This is useful for normalizing transcript lengths
between different haplotypes by focusing on coding regions rather than variable UTRs.

Usage:
    python gtf_converter.py -i input.gtf -f genome.fa -o output.gtf [options]
"""

import argparse
import sys
import logging
from pathlib import Path
import pyranges as pr
import pandas as pd
import pyfaidx


def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def validate_files(input_file, fasta_file):
    """Validate that input files exist and are readable."""
    if not Path(input_file).exists():
        raise FileNotFoundError(f"Input GTF/GFF file not found: {input_file}")
    
    if not Path(fasta_file).exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    
    logging.info(f"Input files validated: {input_file}, {fasta_file}")


def load_annotation(input_file):
    """Load GTF/GFF file and return pyranges object."""
    logging.info(f"Loading annotation from: {input_file}")
    try:
        ann = pr.read_gtf(input_file)
        logging.info(f"Loaded {len(ann)} features from annotation file")
        return ann
    except Exception as e:
        logging.error(f"Failed to load annotation file: {e}")
        raise


def convert_cds_to_exons(ann, extension_5=0, extension_3=0, group_by='gene_id'):
    logging.info(f"Converting CDS to exons with {extension_5}bp 5' and {extension_3}bp 3' extensions")

    cds_selector = (ann.Feature == 'CDS')
    cds = ann[cds_selector]

    if len(cds) == 0:
        logging.warning("No CDS features found in the annotation file")
        return ann

    logging.info(f"Found {len(cds)} CDS features")

    # Extend CDS and convert to exons
    extended_cds = cds.extend({'5': extension_5, '3': extension_3}, group_by=group_by)
    extended_cds.Feature = 'exon'
    exon_df = extended_cds.df.copy()

    # Recalculate transcript boundaries from exons
    transcript_selector = ann.Feature.isin(["mRNA", "transcript"])
    transcripts = ann[transcript_selector].df.copy()

    if len(transcripts) > 0:
        # Get new min/max coordinates per transcript from the exons
        exon_bounds = (
            exon_df.groupby("transcript_id")
            .agg(new_Start=("Start", "min"), new_End=("End", "max"))
            .reset_index()
        )
        transcripts = transcripts.merge(exon_bounds, on="transcript_id", how="left")
        # Only update transcripts that have corresponding exons
        has_exons = transcripts["new_Start"].notna()
        transcripts.loc[has_exons, "Start"] = transcripts.loc[has_exons, "new_Start"].astype(int)
        transcripts.loc[has_exons, "End"]   = transcripts.loc[has_exons, "new_End"].astype(int)
        transcripts.drop(columns=["new_Start", "new_End"], inplace=True)

    # Recalculate gene boundaries from updated transcript boundaries
    gene_selector = ann.Feature == 'gene'
    genes = ann[gene_selector].df.copy()

    if len(genes) > 0 and len(transcripts) > 0:
        gene_bounds = (
            transcripts.groupby(group_by)
            .agg(new_Start=("Start", "min"), new_End=("End", "max"))
            .reset_index()
        )
        genes = genes.merge(gene_bounds, on=group_by, how="left")
        has_transcripts = genes["new_Start"].notna()
        genes.loc[has_transcripts, "Start"] = genes.loc[has_transcripts, "new_Start"].astype(int)
        genes.loc[has_transcripts, "End"]   = genes.loc[has_transcripts, "new_End"].astype(int)
        genes.drop(columns=["new_Start", "new_End"], inplace=True)

    # Combine all features
    parts = [df for df in [genes, transcripts, exon_df] if len(df) > 0]
    combined_df = pd.concat(parts)

    ann_updated = pr.PyRanges(combined_df)
    logging.info(f"Created {len(ann_updated)} features after conversion")

    return ann_updated

def filter_scaffolds(ann, exclude_pattern='scf_'):
    """Remove features from scaffolds matching the exclusion pattern."""
    if exclude_pattern:
        initial_count = len(ann)
        ann_filtered = ann[~ann.Chromosome.str.startswith(exclude_pattern)]
        removed_count = initial_count - len(ann_filtered)
        if removed_count > 0:
            logging.info(f"Removed {removed_count} features from scaffolds starting with '{exclude_pattern}'")
        return ann_filtered
    return ann


def clip_to_genome_bounds(ann, fasta_file):
    """Clip features to genome boundaries using the reference FASTA."""
    logging.info("Clipping features to genome boundaries")
    try:
        pyf = pyfaidx.Fasta(fasta_file)
        ann_clipped = pr.genomicfeatures.genome_bounds(ann, chromsizes=pyf, clip=True)
        logging.info("Successfully clipped features to genome bounds")
        return ann_clipped
    except Exception as e:
        logging.error(f"Failed to clip to genome bounds: {e}")
        return ann


def fix_annotation_attributes(ann_df, transcript_feature="mRNA"):
    """
    Fix various annotation attribute issues for different genome sources.
    
    Args:
        ann_df: DataFrame from PyRanges object
        transcript_feature: The feature type used for transcripts in input
    """
    logging.info("Fixing annotation attributes")
    

    
    # Fix gene_id for gene features if null
    gene_mask = ann_df["Feature"] == 'gene'
    gene_id_null_mask = ann_df["gene_id"].isnull()
    
    if gene_mask.any() and gene_id_null_mask.any():
        logging.info("Fixing null gene_id values using Name attribute")
        mask = gene_mask & gene_id_null_mask
        logging.info(ann_df.loc[mask,].head())
        try:
            ann_df.loc[mask, "gene_id"] = ann_df.loc[mask, "Name"].str.rsplit(".", n=1, expand=True)[0]
        except Exception as e:
            logging.error(f"Failed to fix gene_id from Name attribute: {e}")

    
    # Fix negative start positions
    negative_starts = ann_df["Start"] < 0
    if negative_starts.any():
        logging.info(f"Fixed {negative_starts.sum()} negative start positions")
        ann_df.loc[negative_starts, "Start"] = 0
    
    # Standardize feature names
    ann_df.loc[ann_df["Feature"] == transcript_feature, "Feature"] = "transcript"
    
    # Fix Parent and ID attributes
    transcript_mask = ann_df["Feature"] == "transcript"
    ann_df.loc[transcript_mask, "Parent"] = ann_df.loc[transcript_mask, "gene_id"]
    print(ann_df.head())
    gene_mask = ann_df["Feature"] == "gene"
    ann_df.loc[gene_mask, "ID"] = ann_df.loc[gene_mask, "gene_id"]
    
    return ann_df


def process_gtf(input_file, fasta_file, output_file, 
                extension_5=0, extension_3=0, 
                exclude_scaffolds='scf_', transcript_feature="mRNA",
                group_by='gene_id'):
    """
    Main processing function that converts CDS to exons with extensions.
    
    Args:
        input_file: Path to input GTF/GFF file
        fasta_file: Path to reference genome FASTA file
        output_file: Path for output GTF file
        extension_5: Extension on 5' end (bp)
        extension_3: Extension on 3' end (bp)
        exclude_scaffolds: Pattern to exclude scaffolds (None to disable)
        transcript_feature: Feature type for transcripts in input file
        group_by: Column to group by for extensions
    """
    
    # Validate input files
    validate_files(input_file, fasta_file)
    
    # Load annotation
    ann = load_annotation(input_file)
    
    # Convert CDS to exons with extensions
    ann_converted = convert_cds_to_exons(ann, extension_5, extension_3, group_by)
    
    # Filter scaffolds if requested
    ann_filtered = filter_scaffolds(ann_converted, exclude_scaffolds)
    
    # Clip to genome bounds
    ann_clipped = clip_to_genome_bounds(ann_filtered, fasta_file)
    
    # Fix annotation attributes
    ann_df = ann_clipped.df.copy()
    ann_df = fix_annotation_attributes(ann_df, transcript_feature)
    
    # Create final PyRanges object and sort
    ann_final = pr.PyRanges(ann_df)
    ann_final = ann_final.sort(by=['Chromosome', 'Start'])
    
    # Write output
    logging.info(f"Writing output to: {output_file}")
    ann_final.to_gtf(output_file)
    logging.info(f"Successfully wrote {len(ann_final)} features to output file")
    
    return ann_final


def main():
    """Main function with command-line interface."""
    parser = argparse.ArgumentParser(
        description='Convert CDS features to exon features with extensions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='Input GTF/GFF file path')
    parser.add_argument('-f', '--fasta', required=True,
                       help='Reference genome FASTA file path')
    parser.add_argument('-o', '--output', required=True,
                       help='Output GTF file path')
    
    # Optional arguments
    parser.add_argument('--ext5', type=int, default=0,
                       help='Extension on 5\' end (bp)')
    parser.add_argument('--ext3', type=int, default=0,
                       help='Extension on 3\' end (bp)')
    parser.add_argument('--exclude-scaffolds', default='scf_',
                       help='Exclude scaffolds starting with this pattern (use "none" to disable)')
    parser.add_argument('--transcript-feature', default='mRNA',
                       help='Feature type for transcripts in input file')
    parser.add_argument('--group-by', default='gene_id',
                       help='Column to group by for extensions')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(args.verbose)
    
    # Handle special case for disabling scaffold exclusion
    exclude_scaffolds = None if args.exclude_scaffolds.lower() == 'none' else args.exclude_scaffolds
    
    try:
        # Process the GTF file
        process_gtf(
            input_file=args.input,
            fasta_file=args.fasta,
            output_file=args.output,
            extension_5=args.ext5,
            extension_3=args.ext3,
            exclude_scaffolds=exclude_scaffolds,
            transcript_feature=args.transcript_feature,
            group_by=args.group_by
        )
        
        logging.info("Processing completed successfully!")
        
    except Exception as e:
        logging.error(f"Processing failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
