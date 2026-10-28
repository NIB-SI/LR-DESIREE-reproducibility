
##### ATlantic
# get promotor coordinated form gff
conda activate pyranges
python Promotor_analysis_desiree/scripts/promotor_extractor.py  --gff reference/Atlantic/unitato2Atl.with_chloroplast_and_mito.no_scaffold.gff3    --fasta reference/Atlantic/ATL_v3.asm.with_chloroplast_and_mito.fa     --output reference/Atlantic/unitato2Atl.promoter1500.agat.gtf --length 1500

# extract the promotors form fasta
conda activate gffread
gffread  -g reference/Atlantic/ATL_v3.asm.with_chloroplast_and_mito.fa   -w  reference/Atlantic/unitato2Atl.promoter1500.fa reference/Atlantic/unitato2Atl.promoter1500.agat.gtf


## desiree 
conda activate pyranges
python Promotor_analysis_desiree/scripts/promotor_extractor.py  --gff  reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.gff     --fasta reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa     --output reference/Desiree_v1/De_v1.promoter1500.agat.gtf --length 1500

conda activate gffread
gffread  -g reference/Desiree_v1/De_v1_no_scaffold_chloroplast_mt.fa   -w  reference/Desiree_v1/De_v1.promoter1500.fa reference/Desiree_v1/De_v1.promoter1500.agat.gtf

# add cultivar prefix --> >Soltu... --> >ATL.Soltu
sed -i 's/^>/>ATL./' reference/Atlantic/unitato2Atl.promoter1500.fa 
sed -i 's/^>/>DES./' reference/Desiree_v1/De_v1.promoter1500.fa


SHARED_ASE_DES_ATL=LR_DESIREE_PAPER/ANALYSIS/ASE_results/comparison_DES_ATL/AT_DES_merged_cis_results_control_grouped_all_gene_ids.csv


grep --no-group-separator -A22  -f $SHARED_ASE_DES_ATL reference/Desiree_v1/De_v1.promoter1500.fa > LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/De_v1.promoter1500_filtered.fa

grep --no-group-separator -A22  -f $SHARED_ASE_DES_ATL reference/Atlantic/unitato2Atl.promoter1500.fa > LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/unitato2Atl.promoter1500_filtered.fa



conda activate blast

# align the ATL vs DES promoters with blast
makeblastdb   -dbtype nucl -parse_seqids	-in LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/unitato2Atl.promoter1500_filtered.fa

blastn -query LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/De_v1.promoter1500_filtered.fa -outfmt 6 -db LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/unitato2Atl.promoter1500_filtered.fa -out LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/BLAST_promotor_ATL_vs_DES.tsv



awk -F '\t' '$3 > 95.00 && $4 > 800' LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/BLAST_promotor_ATL_vs_DES.tsv  > LR_DESIREE_PAPER/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/BLAST_promotor_ATL_vs_DES_filtered_all.tsv



