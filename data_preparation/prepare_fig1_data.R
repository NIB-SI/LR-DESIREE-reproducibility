# Prepare the clean inputs for Fig1_short_vs_long_read_comparison.ipynb
#
# Re-runs the DE steps from
#   analysis/notebooks/01_DE_gene_level_long_vs_short_read.ipynb  (DM genome)
#   analysis/notebooks/02_DE_allele_level_long_vs_short_read.ipynb (haplotype-resolved)
# ComBat-seq on long-read counts -> voomLmFit (corrected + uncorrected), tximport + voomLmFit for short reads.
#
# Usage: Rscript prepare_fig1_data.R   (run from this directory, deseq2 conda env)

suppressPackageStartupMessages({
  library(sva); library(edgeR); library(limma); library(tximport)
  library(readr); library(dplyr); library(readxl)
})

OUT_DIR <- "."
LR_DIR  <- "nf-core-plantlongrnaseq"
SR_DIR  <- "Desiree_ASE/output"
ANNOTATION_XLSX <- "reference/ADAPT_liftoff/Phureja_v4-v6.1_translations.xlsx"

sample_names    <- c("DE_C_H_13", "DE_C_H_14", "DE_C_H_15", "DE_C_H_16",
                     "DE_H_H_17", "DE_H_H_18", "DE_H_H_19", "DE_H_H_20")
library_methods <- c(1, 1, 1, 1, 2, 2, 1, 2)
conditions      <- rep(c("control", "heat"), c(4, 4))
group           <- factor(conditions)

# --- functional annotation (DMv6.1), first entry per gene as in the Fig3 table ---
annotation <- read_excel(ANNOTATION_XLSX)[, 2:3]
colnames(annotation) <- c("gene_id_unitato", "functional_annotation")
annotation <- annotation[!duplicated(annotation$gene_id_unitato), ]

# Soltu.DM.01G000020.Hap2 -> Soltu.DM.01G000020 ; Soltu.DM.03G023500_1.Hap2 -> Soltu.DM.03G023500
to_unitato_id <- function(ids, haplotype_resolved) {
  if (haplotype_resolved) ids <- sub("\\.Hap[0-9]+$", "", ids)
  sub("_[0-9]+$", "", ids)
}

run_voom <- function(counts, min_cpm, min_samples) {
  dge <- DGEList(counts = counts, group = conditions)
  dge <- dge[rowSums(cpm(dge) > min_cpm) >= min_samples, , keep.lib.sizes = FALSE]
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- voomLmFit(dge, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(heat - control, levels = design)))
  res <- topTable(fit, sort.by = "P", n = Inf)
  res$gene_id <- rownames(res)
  cpms <- cpm(dge)
  res$avg_control <- rowMeans(cpms[res$gene_id, conditions == "control"])
  res$avg_heat    <- rowMeans(cpms[res$gene_id, conditions == "heat"])
  res
}

annotate_and_write <- function(res, haplotype_resolved, file) {
  res$gene_id_unitato <- to_unitato_id(res$gene_id, haplotype_resolved)
  res <- left_join(res, annotation, by = "gene_id_unitato")
  cols <- c("gene_id", if (haplotype_resolved) "gene_id_unitato",
            "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
            "avg_control", "avg_heat", "functional_annotation")
  write_csv(res[, cols], file.path(OUT_DIR, file))
  cat(file, nrow(res), "genes,", sum(!is.na(res$functional_annotation)), "annotated\n")
}

prepare_level <- function(label, lr_run, sr_run, min_cpm, min_samples, haplotype_resolved) {
  # long-read: ComBat-seq corrected and uncorrected
  uncorrected <- read.table(file.path(LR_DIR, lr_run, "tximeta/all_samples.gene_counts.tsv"),
                            header = TRUE, sep = "\t", as.is = c(1, 2))
  counts <- as.matrix(uncorrected[, sample_names])
  rownames(counts) <- uncorrected$gene_id
  corrected <- ComBat_seq(counts = counts, batch = library_methods,
                          group = as.integer(group))
  rownames(corrected) <- uncorrected$gene_id

  annotate_and_write(run_voom(corrected, min_cpm, min_samples), haplotype_resolved,
                     sprintf("DE_%s_long_corrected.csv", label))
  annotate_and_write(run_voom(counts, min_cpm, min_samples), haplotype_resolved,
                     sprintf("DE_%s_long_uncorrected.csv", label))

  # short-read: salmon -> tximport
  files <- file.path(SR_DIR, sr_run, "salmon", sample_names, "quant.sf")
  tx2gene <- read_tsv(file.path(SR_DIR, sr_run, "salmon/tx2gene.tsv"), show_col_types = FALSE)
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  colnames(txi$counts) <- sample_names
  annotate_and_write(run_voom(txi$counts, min_cpm, min_samples), haplotype_resolved,
                     sprintf("DE_%s_short.csv", label))
}

# DM genome: CPM > 3 in >= 5 samples ; haplotype-resolved: CPM > 4 in >= 6 samples (as in the analysis notebooks)
prepare_level("DM",  "output_desiree_unitato",        "nf_core_rnaseq_short_read_heat_unitato",                 3, 5, FALSE)
prepare_level("hap", "output_desiree_liftoff_phased", "nf_core_rnaseq_short_read_heat_phased_desiree_liftoff", 4, 6, TRUE)

# short-read mapping rate for the panel F bins (single-isoform genes, all 8 samples): see prepare_fig1_panelB.py

file.copy("LR_DESIREE_PAPER/DEG_results/avgTxlengths_desiree_unitato.tsv", OUT_DIR, overwrite = TRUE)
