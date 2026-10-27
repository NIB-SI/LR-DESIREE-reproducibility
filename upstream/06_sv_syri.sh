#!/usr/bin/env bash
# Structural variants between the four Désirée haplotypes (all six pairs):
# minimap2 whole-haplotype alignment -> SyRI SV calling -> SVs overlapping promoters
# (2,000 bp upstream / 500 bp downstream of the TSS) of syntenic gene pairs.
#
# Env: seqkit, minimap2, samtools, syri, python3
set -euo pipefail

BIN=$(cd "$(dirname "$0")" && pwd)/bin
WD=/scratch/nadjafn/SV_whole_genome_desiree

ref=/scratch/nadjafn/reference/Desiree_v1/De_v1.fa
gff=/scratch/nadjafn/reference/Desiree_v1/De_v1.unitato_liftoff_haplotap.gff3
fasta_dir=$WD/individual_seqs
alignment_dir=$WD/alignments
syri_out=$WD/output/syri_output
parse_out=$WD/output/syri_parse
mkdir -p "$fasta_dir" "$alignment_dir" "$syri_out" "$parse_out"

# one FASTA per haplotype
for hap in 1 2 3 4; do
    seqkit grep -r -p "chr_\d+_${hap}" "$ref" > "$fasta_dir/hap${hap}.fasta"
done

pairs=(
    "hap1 hap2"
    "hap1 hap3"
    "hap1 hap4"
    "hap2 hap3"
    "hap2 hap4"
    "hap3 hap4"
)

for pair in "${pairs[@]}"; do
    read -r ref_hap qry_hap <<< "$pair"
    prefix=${ref_hap}_${qry_hap}_asm20_optimized
    echo "=== ${qry_hap} vs ${ref_hap} ==="

    minimap2 -ax asm20 -k15 -s 50 -z 1000 -t 16 --eqx -r3k,100k -O3,20 \
        "$fasta_dir/${ref_hap}.fasta" "$fasta_dir/${qry_hap}.fasta" \
        | samtools sort -o "$alignment_dir/${prefix}.bam"

    mkdir -p "$syri_out/$prefix"
    syri -c "$alignment_dir/${prefix}.bam" \
         -r "$fasta_dir/${ref_hap}.fasta" \
         -q "$fasta_dir/${qry_hap}.fasta" \
         -F B -f --prefix "$prefix" --dir "$syri_out/$prefix"

    python3 "$BIN/get_promoter_svs_syri.py" \
        -o "$parse_out/${prefix}_promoter_sv_pairs_f_2000up_500down.tsv" \
        --upstream 2000 \
        --ref-downstream 500 \
        --query-upstream 2000 \
        --query-downstream 500 \
        "$gff" \
        "$syri_out/$prefix/${prefix}syri.vcf"
done
