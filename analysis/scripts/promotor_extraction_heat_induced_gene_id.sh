
# ### Desiree
# # get the promotor sequences from the gff file

# # get the longest isofrom 
# conda activate agat
# agat_sp_keep_longest_isoform.pl -gff nf-core-plantlongrnaseq/output_desiree_liftoff_phased/gff/combined.fixed.gtf  -o nf-core-plantlongrnaseq/output_desiree_liftoff_phased/gff/combined.fixed.longest_isoforms.gtf

# conda activate gffread
# gffread -E nf-core-plantlongrnaseq/output_desiree_liftoff_phased/gff/combined.fixed.longest_isoforms.gtf -T -o nf-core-plantlongrnaseq/output_desiree_liftoff_phased/gff/combined.fixed.longest_isoforms.gff3



# DIR=Promotor_analysis_desiree
# conda activate pyranges
# python Promotor_analysis_desiree/scripts/promotor_extractor_exon.py  --gff nf-core-plantlongrnaseq/output_desiree_liftoff_phased/gff/combined.fixed.longest_isoforms.gtf  --fasta reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa  --output reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_2000_promotor.gff3 --length 2000



# conda activate gffread
# # # to gffread gff
# gffread reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_2000_promotor.gff3 > reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_2000_promotor_gffread.gff3


# # # extract the fasta of promotors with the gff file
# # conda activate gffread
# gffread -g reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa  -w  $DIR/output/desiree_2000_up_down_stream_gene_ids.fa reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_2000_promotor_gffread.gff3



# Get the promotor sequences of heat upregulated allles
conda activate minimap2

# we use the transcript ids for extracting the promtor sequences as this contains the UTRs also, but hen we extract also empty sequences so we need to filter them then

DIR=LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences

HEAT_UP_GENES="LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/upregulated_genes_voom_lr_combat_seq_corrected_lFC2_adjP0.01._tx_ids.tsv"

BACKGROUND_GENES="LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/background_genes_tx_ids.tsv"

# filter
# filter the list of genes for the ones that are present in the promtor fasta file 
grep -f  $HEAT_UP_GENES Promotor_analysis_desiree/output/desiree_2000_up_down_stream_transcript_ids.fa | sed 's/>//g'  > LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/upregulated_genes_voom_lr_combat_seq_corrected_lFC2_adjP0.01._tx_ids.filtered.tsv


grep -f  $BACKGROUND_GENES Promotor_analysis_desiree/output/desiree_2000_up_down_stream_transcript_ids.fa | sed 's/>//g' > LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/background_genes_tx_ids.filtered.tsv


HEAT_UP_GENES="LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/upregulated_genes_voom_lr_combat_seq_corrected_lFC2_adjP0.01._tx_ids.filtered.tsv"



BACKGROUND_GENES="LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/background_genes_tx_ids.filtered.tsv"

## now get the promotors
seqkit grep -i -f $HEAT_UP_GENES Promotor_analysis_desiree/output/desiree_2000_up_down_stream_transcript_ids.fa  > $DIR/desiree_heat_log2FC2_upregulated_2000_promotors.fa

# 

seqkit grep -i -f $BACKGROUND_GENES  Promotor_analysis_desiree/output/desiree_2000_up_down_stream_transcript_ids.fa > $DIR/desiree_heat_downregulated_not_changing_2000_promotors.fa



conda activate minimap2

DIR=LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences
# also get the kunitz genes
KUNITZ_GENES="LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/promotor_enrichment_input_sequences/kunitz_genes_tx_ids.tsv"

grep --no-group-separator -A29 -f $KUNITZ_GENES  Promotor_analysis_desiree/output/desiree_2000_up_down_stream_gene_ids.fa > $DIR/desiree_kunitz_2000_promotors.fa


# other NAOD genes
grep -E --no-group-separator -A29 "Soltu.DM.08G0228[0,2,3]" Promotor_analysis_desiree/output/desiree_2000_up_down_stream_gene_ids.fa 