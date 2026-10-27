#!/usr/bin/env bash
# Transposable-element class of promoter SVs: extract the sequence of each selected
# SyRI SV (> 50 bp, the only large SV in its promoter) and BLAST it against the
# EDTA TE libraries of the four Désirée haplotypes.
#
# Input: SV_IDS, list of selected SV ids (written by the promoter SV analysis notebook)
# Output: TE_sv50_blast_results_with_header.txt
# Env: blast (makeblastdb, blastn); GNU awk
set -euo pipefail

WD=/scratch/nadjafn/SV_whole_genome_desiree
SV_IDS=/scratch/nadjafn/LR_DESIREE_PAPER/ANALYSIS/ASE_results/SV_analysis/one_sv_50_id_chrom.csv
SYRI_DIR=$WD/output/syri_output
OUT_SEQ_DIR=$WD/output/SV_sequences
TE_DIR=$WD/output/TE_analysis
EDTA_LIBS=(/scratch/timg/desiree_annotation/edta/De_v1_hap*/*TElib.fa)
mkdir -p "$OUT_SEQ_DIR" "$TE_DIR"

> "$OUT_SEQ_DIR/sv_sequences_all.fa"

for HAP_DIR in "$SYRI_DIR"/hap*_asm20_optimized; do
    PREFIX=$(basename "$HAP_DIR")
    VCF_FILE="$HAP_DIR/${PREFIX}syri.vcf"
    [ -f "$VCF_FILE" ] || continue

    # unique SV id: CHROM_ID_ChrB
    IDS_VCF="$HAP_DIR/${PREFIX}syri_ids.vcf"
    awk -F'\t' 'BEGIN{OFS="\t"} {
        chrB = ""
        if (match($8, /ChrB=([^;]+)/, arr)) chrB = "_" arr[1]
        print $0, $1"_"$3 chrB
    }' "$VCF_FILE" > "$IDS_VCF"

    FILTERED_VCF="$HAP_DIR/${PREFIX}syri_ids_one_sv_50.vcf"
    grep -f "$SV_IDS" "$IDS_VCF" > "$FILTERED_VCF" || true

    # SV sequence = the longer of REF/ALT without the anchor base; skip symbolic alleles
    awk -F'\t' '
      $5 !~ /^</ {
        seq = (length($5) > length($4)) ? substr($5, 2) : substr($4, 2)
        print ">" $NF "\n" seq
      }
    ' "$FILTERED_VCF" >> "$OUT_SEQ_DIR/sv_sequences_all.fa"
done

cat "${EDTA_LIBS[@]}" > "$TE_DIR/De_v1_all_haps_chrs.fa.mod.EDTA.TElib.fa"
makeblastdb -in "$TE_DIR/De_v1_all_haps_chrs.fa.mod.EDTA.TElib.fa" -dbtype nucl -out "$TE_DIR/TE_lib"

blastn -query "$OUT_SEQ_DIR/sv_sequences_all.fa" -db "$TE_DIR/TE_lib" \
    -outfmt 6 -evalue 1e-5 -num_threads 20 > "$TE_DIR/TE_sv50_blast_results.txt"

echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tsgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
    | cat - "$TE_DIR/TE_sv50_blast_results.txt" > "$TE_DIR/TE_sv50_blast_results_with_header.txt"
