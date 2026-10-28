# Upstream processing

Scripts that produce the inputs of the figure notebooks, in the order they were run.
Paths point to the original locations on our cluster; adjust them before running.

| Script | Step | Tool |
|---|---|---|
| `01_liftoff_haplotypes.sh` | Lift the UniTato annotation onto the four haplotypes of Désirée / Atlantic and tag gene IDs with the haplotype | Liftoff 1.6.3 |
| `02_syntelogfinder.sh` | Group haplotype genes into syntelogs (`Synt_id`, synteny category) | [syntelogfinder](https://github.com/NIB-SI/syntelogfinder) |
| `03_long_read_rnaseq.sh` | Long-read RNA-seq: Désirée (PacBio Kinnex), Atlantic (ONT) | [nf-core/plantlongrnaseq](https://github.com/nadjano/nf-core-plantlongrnaseq) |
| `04_short_read_rnaseq.sh` | Short-read RNA-seq, Salmon quantification | nf-core/rnaseq 3.17.0 |
| `05_atacseq.sh` | ATAC-seq | nf-core/atacseq 2.1.2 |
| `06_sv_syri.sh` | SVs between all haplotype pairs and their overlap with promoters | minimap2, SyRI |
| `07_te_analysis.sh` | TE class of promoter SVs | BLAST against EDTA TE libraries |

`bin/` holds the helper scripts called above, `config/` the Nextflow configs, `params/` the
syntelogfinder parameter files and `samplesheets/` the pipeline sample sheets (file names only;
point them to your copies of the raw data).

Allele-specific expression, differential expression and isoform analyses were run on these
outputs with [polyase](https://github.com/nadjano/polyase); see `../analysis/`.
