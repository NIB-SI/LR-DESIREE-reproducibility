# LR-DESIREE-reproducibility

Reproducibility notebooks for:

**Long-read RNA-seq and haplotype-resolved reference genome resolve haplotype- and isoform-level heat responses in tetraploid potato**

## Contents

- `Fig1_short_vs_long_read_comparison.ipynb`
- `Fig2_ASE_SV_promoter_analysis_main.ipynb`
- `FigX_ASE_SV_promoter_analysis_supplement.ipynb` (supplementary panels for Fig. 2)
- `Fig3_heat_responsive_alleles_and_motifs.ipynb`
- `Fig4_HSP_Kunitz_supplement.ipynb`
- `Fig5_isoform_expression.ipynb`
- `upstream/`: scripts for the upstream processing (annotation liftoff, syntelogs, RNA-seq,
  ATAC-seq, structural variants, TE analysis); see `upstream/README.md`

## Quality control

- `multiqc_report.html`: MultiQC report of the long-read RNA-seq (PacBio Kinnex) processing with
  nf-core/plantlongrnaseq against the haplotype-resolved Désirée genome (FastQC, Kraken, samtools
  alignment statistics, software versions). Download and open it in a web browser.
