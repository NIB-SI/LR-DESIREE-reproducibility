#!/usr/bin/env bash
# Long-read RNA-seq (Désirée: PacBio Kinnex, Atlantic: ONT) with nf-core/plantlongrnaseq:
# https://github.com/nadjano/nf-core-plantlongrnaseq
#
# CDS features are first converted to exons (+ flanking sequence) so that reads are
# assigned to the CDS-based liftoff models.
#
# Run from a clone of the nf-core-plantlongrnaseq repository.
# Env: pyranges (CDS2exon conversion), nextflow + singularity
set -euo pipefail

BIN=$(cd "$(dirname "$0")" && pwd)/bin
eval "$(conda shell.bash hook)"

# ---- Désirée ----
DES=reference/Desiree_v1

conda activate pyranges
python "$BIN/convert_gtf_CDS2exon.py" \
    -i $DES/De_v1.unitato_liftoff_haplotap_gffread.with_chloroplast_and_mito.gtf \
    -o $DES/De_v1.unitato_liftoff_haplotap_gffread.with_chloroplast_and_mito.CDS2exon.gtf \
    -f $DES/De_v1_no_scaffold_chloroplast_mt.fa

conda activate sqanti3
nextflow run main.nf -resume -profile singularity \
    --input assets/samplesheet_desiree.csv \
    --outdir output_desiree_liftoff_phased_CDS2exon \
    --fasta $DES/De_v1_no_scaffold_chloroplast_mt.fa \
    --gtf $DES/De_v1.unitato_liftoff_haplotap_gffread.with_chloroplast_and_mito.CDS2exon.gtf \
    --centrifuge_db centrifuge/dbs_v2018/ \
    --sqanti_dir sqanti3/SQANTI3 \
    --technology PacBio \
    --skip_sqanti

# ---- Atlantic ----
ATL=reference/Atlantic

conda activate pyranges
python "$BIN/convert_gtf_CDS2exon.py" \
    -i $ATL/unitato2Atl.with_chloroplast_and_mito.no_scaffold.agat.gtf \
    -o $ATL/unitato2Atl.with_chloroplast_and_mito.no_scaffold.agat.CDS2exon.gtf \
    -f $ATL/ATL_v3.asm.with_chloroplast_and_mito.fa \
    --transcript-feature transcript

conda activate sqanti3
nextflow run main.nf -resume -profile singularity \
    --input assets/samplesheet_atlantic.csv \
    --outdir output_atlantic_liftoff_CDS2exon \
    --fasta $ATL/ATL_v3.asm.with_chloroplast_and_mito.fa \
    --gtf $ATL/unitato2Atl.with_chloroplast_and_mito.no_scaffold.agat.CDS2exon.gtf \
    --centrifuge_db centrifuge/dbs_v2018/ \
    --sqanti_dir sqanti3/SQANTI3 \
    --technology ONT \
    --downsample_rate 0.99 \
    --skip_sqanti
