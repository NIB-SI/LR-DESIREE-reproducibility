#!/usr/bin/env bash
# Short-read RNA-seq (Désirée, control + heat) with nf-core/rnaseq 3.17.0, Salmon
# quantification against the haplotype-tagged liftoff annotation.
# config/rnaseq.config keeps duplicate transcripts in the Salmon index (--keepDuplicates),
# so identical alleles are not collapsed.
#
# Env: nextflow + apptainer
set -euo pipefail

CONFIG=$(cd "$(dirname "$0")" && pwd)/config/rnaseq.config

nextflow run nf-core/rnaseq -r 3.17.0 -resume \
    -profile apptainer \
    -c "$CONFIG" \
    --input /scratch/nadjafn/Desiree_ASE/input/NGS_stu_rna-seq_2024_timg_desiree-ladyrosetta_ADAPT/samplesheet_heat.csv \
    --outdir /scratch/nadjafn/Desiree_ASE/output/nf_core_rnaseq_short_read_heat_phased_desiree_liftoff \
    --gtf /scratch/nadjafn/reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.gtf \
    --fasta /scratch/nadjafn/reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
    --pseudo_aligner salmon \
    --skip_biotype_qc
