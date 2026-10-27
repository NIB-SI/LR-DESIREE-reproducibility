#!/usr/bin/env bash
# Lift the UniTato (DM v6.1) annotation onto each of the four haplotypes of a
# haplotype-resolved genome with Liftoff, then merge and tag gene IDs with the
# haplotype (Soltu.DM.01G000010 -> Soltu.DM.01G000010.Hap1).
#
# Usage: bash 01_liftoff_haplotypes.sh desiree|atlantic
# Env:   liftoff (v1.6.3), seqkit, python3
set -euo pipefail

GENOME=${1:?usage: $0 desiree|atlantic}
BIN=$(dirname "$0")/bin

# reference annotation (UniTato on DM v6.1)
REF_FA=/scratch/nadjafn/reference/ADAPT_liftoff/UniTato_nuc_mt_plast.manualfix.nfn.fa
REF_GFF=/scratch/nadjafn/reference/ADAPT_liftoff/UniTato_nuc_mt_plast.manualfix.nfn_chr.gff3

case $GENOME in
    desiree)
        TARGET_FA=/scratch/nadjafn/reference/Desiree_v1/De_v1.fa
        OUT=/scratch/nadjafn/reference/Desiree_v1/liftoff_unitato
        PREFIX=unitato2De
        chr_name() { echo "chr_${1}_${2}"; }   # chr_01_1
        ;;
    atlantic)
        TARGET_FA=/scratch/nadjafn/reference/Atlantic/ATL_v3.asm.fa
        OUT=/scratch/nadjafn/Atlantic_reference_liftoff/output/liftoff_unitato
        PREFIX=unitato2Atl
        chr_name() { echo "chr${1}_${2}"; }    # chr01_1
        ;;
    *) echo "unknown genome: $GENOME" >&2; exit 1 ;;
esac

mkdir -p "$OUT"
printf "mRNA\nexon\nCDS\nfive_prime_UTR\nthree_prime_UTR\n" > "$OUT/feature_types.txt"

chromosomes=(01 02 03 04 05 06 07 08 09 10 11 12)

for haplotype in 1 2 3 4; do
    # one FASTA per haplotype
    hap_fa=$OUT/${GENOME}_hap${haplotype}.fa
    seqkit grep -r -p "^$(chr_name '\d+' "$haplotype")$" "$TARGET_FA" > "$hap_fa"

    # reference chromosome -> haplotype chromosome
    chr_map_file=$OUT/chr_map_${haplotype}.csv
    > "$chr_map_file"
    for chr in "${chromosomes[@]}"; do
        echo "chr${chr},$(chr_name "$chr" "$haplotype")" >> "$chr_map_file"
    done

    liftoff -g "$REF_GFF" \
            -o "$OUT/${PREFIX}_${haplotype}_a0.9s0.9.gff" \
            -f "$OUT/feature_types.txt" \
            -chroms "$chr_map_file" \
            -u "$OUT/unmapped_features_${haplotype}.txt" \
            -dir "$OUT/intermediate_files_${haplotype}" \
            -p 24 \
            -a 0.9 \
            -s 0.9 \
            -copies \
            "$hap_fa" \
            "$REF_FA"
done

# merge haplotypes and add .HapN to gene/transcript IDs (haplotype taken from the chromosome suffix)
cat "$OUT"/${PREFIX}_[1-4]_a0.9s0.9.gff > "$OUT/${PREFIX}_a0.9s0.9.gff"
python3 "$BIN/update_gff_haplotypes.py" "$OUT/${PREFIX}_a0.9s0.9.gff" "$OUT/${PREFIX}_a0.9s0.9_haptag.gff"
