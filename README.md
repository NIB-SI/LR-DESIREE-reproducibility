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
  - `01_DE_gene_level_*` / `02_DE_allele_level_*`: tximport, ComBat-seq batch correction, CPM
    filtering, voom/limma (voomLmFit, eBayes) and clusterProfiler enrichment
  - `03_gene_count_comparison_*`: long- vs. short-read count comparison
  - `04_motif_enrichment_gene_sets.ipynb`: allele-level DE and the gene sets for the motif enrichment
  - `05_gene_set_enrichment.ipynb`: gene set enrichment of DEGs
  - `06_ASE_long_read.ipynb`: allele-specific expression (PolyASE), control-condition and pairwise tests
  - `07_ASE_short_read.ipynb`: the same with short-read EM counts; writes the multimapping ratios
  - `08_isoform_usage_and_structure.ipynb`: differential isoform usage and isoform structure
  - `09_ASE_atlantic.ipynb`, `10_ASE_desiree_vs_atlantic.ipynb`: Atlantic ASE and the
    Désirée vs. Atlantic comparison
  - `11_peptidase_promoter_alignment.ipynb`: promoter alignment for Fig. 3F
- `analysis/scripts/`: promoter extraction (`promoter_extractor*.py`), 2 kb promoters of the
  heat-induced and non-changing alleles (`extract_promoters_heat_induced.sh`), and 1.5 kb
  promoters with BLASTn for Désirée vs. Atlantic (`extract_promoters_desiree_atlantic_blast.sh`)
- `data_preparation/`: scripts that build the input tables of the figure notebooks from the
  outputs above

Motif enrichment of the two promoter sets was run on the [MEME Suite](https://meme-suite.org)
web server. Haplotype tagging of reads for the JBrowse views is part of
[nf-core/plantlongrnaseq](https://github.com/nadjano/nf-core-plantlongrnaseq).

## Quality control

- `multiqc_report.html`: MultiQC report of the long-read RNA-seq (PacBio Kinnex) processing with
  nf-core/plantlongrnaseq against the haplotype-resolved Désirée genome (FastQC, Kraken, samtools
  alignment statistics, software versions). Download and open it in a web browser.
