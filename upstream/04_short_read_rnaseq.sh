#!/usr/bin/env bash
# Short-read RNA-seq (Désirée, control + heat) with nf-core/rnaseq 3.17.0, Salmon
# quantification against the haplotype-tagged liftoff annotation.
# config/rnaseq.config keeps duplicate transcripts in the Salmon index (--keepDuplicates),
# so identical alleles are not collapsed.
#
# Env: nextflow + apptainer
set -euo pipefail

CONFIG=$(cd "$(dirname "$0")" && pwd)/config/rnaseq.config
SHEETS=$(cd "$(dirname "$0")" && pwd)/samplesheets

nextflow run nf-core/rnaseq -r 3.17.0 -resume \
    -profile apptainer \
    -c "$CONFIG" \
    --input "$SHEETS/samplesheet_shortread.csv" \
    --outdir Desiree_ASE/output/nf_core_rnaseq_short_read_heat_phased_desiree_liftoff \
    --gtf reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.gtf \
    --fasta reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa \
    --pseudo_aligner salmon \
    --skip_biotype_qc
