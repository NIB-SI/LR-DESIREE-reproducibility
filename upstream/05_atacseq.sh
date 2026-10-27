#!/usr/bin/env bash
# ATAC-seq (Désirée) with nf-core/atacseq 2.1.2. Multi-mapping reads are kept so that
# near-identical haplotype regions are not lost.
#
# Env: nextflow + singularity
set -euo pipefail

CONFIG=$(cd "$(dirname "$0")" && pwd)/config/atacseq.config

nextflow run nf-core/atacseq -r 2.1.2 -resume \
    -profile singularity \
    -c "$CONFIG" \
    --input /scratch/nadjafn/Desiree_ASE/input/NGS_stu_dna-seq_2024_timg_desiree-ladyrosetta_ATACseq_ADAPT/samplesheet_ATAC.csv \
    --outdir /scratch/nadjafn/Desiree_ASE/output/nf_core_atacseq \
    --gff /scratch/nadjafn/reference/Desiree_v1/Desiree.longest_isoforms_nodigits.gff \
    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
    --read_length 150 \
    --keep_multi_map \
    --mito_name "NC_008096.2|OR906543.1|OR906544.1|OR906545.1|OR906546.1|OR906547.1" \
    --ataqv_mito_reference NC_008096.2
