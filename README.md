# LR-DESIREE-reproducibility

Reproducibility notebooks for:

**Long-read RNA-seq and haplotype-resolved reference genome resolve haplotype- and isoform-level heat responses in tetraploid potato**

## Contents

Figure notebooks:

- `Fig1_short_vs_long_read_comparison.ipynb`
- `Fig2_ASE_SV_promoter_analysis_main.ipynb`
- `FigS_ASE_SV_promoter_analysis_supplement.ipynb` (supplementary panels S1–S4)
- `Fig3_heat_responsive_alleles_and_motifs.ipynb`
- `Fig4_HSP_Kunitz_supplement.ipynb`
- `Fig5_isoform_expression.ipynb`

Code that produces the figure inputs:

- `upstream/`: annotation liftoff, syntelogs, RNA-seq, ATAC-seq, structural variants and TE
  analysis, with configs and sample sheets; see `upstream/README.md`
- `analysis/notebooks/`: analyses run on the pipeline outputs (outputs cleared)
  - `0_DEGene_*` / `0_DEAllele_*`: tximport, ComBat-seq batch correction, CPM filtering,
    voom/limma (voomLmFit, eBayes) and clusterProfiler enrichment, gene and allele level
  - `1_gene_count_comparision_*`: long- vs. short-read count comparison
  - `2_Motif_enrichment_*`: allele-level DE and promoter sets for the motif enrichment
  - `4_gene_set_enrichment.ipynb`: gene set enrichment of DEGs
  - `potato_polyase_desiree_batch_corr_CDS.ipynb`: allele-specific expression (PolyASE),
    control-condition and pairwise tests
  - `short_potato_polyase.ipynb`: the same with short-read EM counts; writes the multimapping ratios
  - `DIU_potato_polyase_desiree_batch_corr.ipynb`: differential isoform usage and isoform structure
  - `6_Atlantic_ASE.ipynb`, `7_Desiree_Atlantic_ASE_comparison.ipynb`: Atlantic ASE and the
    Désirée vs. Atlantic comparison
  - `peptidase_alignment.ipynb`: promoter alignment for Fig. 3F
- `analysis/scripts/`: promoter extraction (`promotor_extractor*.py`, 2 kb promoters of the
  heat-induced and non-changing alleles, 1.5 kb promoters and BLASTn for Désirée vs. Atlantic)
  and haplotype tagging of reads for the JBrowse views (`primary_alignment_haplotag.py`)
- `data_preparation/`: scripts that build the input tables of the figure notebooks from the
  outputs above

## Quality control

- `multiqc_report.html`: MultiQC report of the long-read RNA-seq (PacBio Kinnex) processing with
  nf-core/plantlongrnaseq against the haplotype-resolved Désirée genome (FastQC, Kraken, samtools
  alignment statistics, software versions). Download and open it in a web browser.
