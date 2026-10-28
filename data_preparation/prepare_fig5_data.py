# Prepare the inputs for clean_notebooks/Fig5_isoform_expression.ipynb
#
# Re-runs the isoform-level analysis of
#   ANALYSIS/notebooks/Polyase_analysis/DIU_potato_polyase_desiree_batch_corr.ipynb
# and keeps only what the figure needs:
#   DIU_results.csv                               differential isoform usage heat vs. control, all isoforms (panel C)
#   isoform_structure_results_control.csv        differential isoform structure between haplotypes, control (panel A)
#   DIU_plotting_selected_genes.pkl               DIU plotting table for the panel D-I genes (isoform structure plots)
#   isoform_structure_plotting_heat_ABH.pkl       allele-specific isoform structure plotting table, Soltu.DM.02G008950 (panel B)
#   annotation_selected_genes.parquet             transcript/exon annotation for the genes above (structure plots)
#   isoform_expression_selected_genes.h5ad        batch-corrected transcript counts/CPM for the panel D-I genes
#
# Usage: python prepare_fig5_data.py   (run from this directory, polyase conda env; takes a while)
import numpy as np
import pandas as pd
import polars as pl
import polyase as poly
import RNApysoforms as RNApy

LR_DIR = "nf-core-plantlongrnaseq/output_desiree_liftoff_phased"

DIU_GENES = ["Soltu.DM.09G020350",  # HSFA6b
             "Soltu.DM.01G000760",  # JAZ3
             "Soltu.DM.12G029670",  # SPA1-related 3
             "Soltu.DM.12G023350",  # COP1
             "Soltu.DM.05G008530",  # carbonic anhydrase
             "Soltu.DM.03G037350"]  # HSP40
STRUCTURE_GENE = "Soltu.DM.02G008950"  # alpha/beta hydrolase

SYNTELOG_CATEGORIES = ("1hap1_1hap2_1hap3_1hap4_s|0hap1_1hap2_1hap3_1hap4_s|1hap1_0hap2_1hap3_1hap4_s|"
                       "1hap1_1hap2_0hap3_1hap4_s|1hap1_1hap2_1hap3_0hap4_s|0hap1_0hap2_1hap3_1hap4_s|"
                       "1hap1_1hap2_0hap3_0hap4_s|0hap1_1hap2_1hap3_0hap4_s|1hap1_0hap2_0hap3_1hap4_s")

# --- load transcript counts (source cells 5-7) ---
sample_info = {s: "control" for s in ["DE_C_H_13", "DE_C_H_14", "DE_C_H_15", "DE_C_H_16"]} | \
              {s: "heat" for s in ["DE_H_H_17", "DE_H_H_18", "DE_H_H_19", "DE_H_H_20"]}
adata_isoform = poly.load_ase_data(
    var_obs_file="potato-allelic-orthogroups/output_desiree_liftoff/genespace/De_v1.unitato_liftoff_haplotap_genespace_categories.tsv",
    sample_info=sample_info,
    tx_to_gene_file=f"{LR_DIR}/custom/combined.tx2gene.tsv",
    isoform_counts_dir=f"{LR_DIR}/oarfish/",
    quant_dir=f"{LR_DIR}/oarfish/",
)

# --- replace counts by ComBat-seq batch-corrected counts (cell 9) ---
# columns are matched by sample name: the sample order differs between these files and adata.obs
# (the source notebook assigned them by position, which swapped replicates within each condition)
em = pd.read_csv(f"{LR_DIR}/tximeta/all_samples.transcript_counts_combat_seq_corrected.tsv", sep="\t", index_col=0)
em = em.drop(columns=["gene_id"])
unique = pd.read_csv(f"{LR_DIR}/oarfish_batch_corrected/allelic_unique_transcript_counts.tsv", sep="\t", index_col=0)
ambiguous = pd.read_csv(f"{LR_DIR}/oarfish_batch_corrected/allelic_ambiguous_transcript_counts.tsv", sep="\t", index_col=0)
for layer, counts in [("em_counts", em), ("unique_counts", unique), ("ambiguous_counts", ambiguous)]:
    adata_isoform.layers[layer] = counts.reindex(adata_isoform.var_names).fillna(0)[adata_isoform.obs_names].values.T

lib_sizes = adata_isoform.layers["em_counts"].sum(axis=1)[:, np.newaxis]
lib_sizes[lib_sizes == 0] = 1
for layer in ["em", "unique", "ambiguous"]:
    adata_isoform.layers[f"{layer}_cpm"] = adata_isoform.layers[f"{layer}_counts"] / lib_sizes * 1e6

# --- functional annotation (cell 11) ---
annotation = pd.read_excel("reference/ADAPT_liftoff/Phureja_v4-v6.1_translations.xlsx").iloc[:, 1:3]
annotation.columns = ["gene_id_unitato", "functional_annotation"]
var = adata_isoform.var.copy()
var["gene_id_unitato"] = var["gene_id"].str.rsplit(".", n=1).str[0]
var["original_index"] = var.index
var = var.merge(annotation, on="gene_id_unitato", how="left")
var = var.drop_duplicates(subset=["original_index", "gene_id"], keep="first").set_index("original_index")
var.index.name = None
adata_isoform.var = var

# --- differential isoform usage, heat vs. control (cells 36-39) ---
adata_filtered = poly.filter_low_expressed_genes(adata_isoform, min_expression=100.0, mode="all",
                                                 group_col="gene_id", verbose=True, layer="em_counts")
adata_filtered = poly.calculate_allelic_ratios(adata_filtered, "unique_counts")
adata_filtered = poly.calculate_multi_ratios(adata_filtered, "unique_counts", "ambiguous_counts")
adata_filtered = poly.calculate_per_allele_ratios(adata_filtered, unique_layer="unique_counts", multi_layer="ambiguous_counts")
mask = adata_filtered.var["tx_multimapping_ratio_per_allele_weighted_average"] <= 1
DIU_results, DIU_plotting = poly.test_isoform_DIU_between_conditions(adata_filtered[:, mask].copy(),
                                                                     layer="em_counts", inplace=True)
DIU_results.to_csv("DIU_results.csv", index=False)
DIU_plotting[DIU_plotting["gene_id"].str.contains("|".join(DIU_GENES), na=False)] \
    .to_pickle("DIU_plotting_selected_genes.pkl")
print("DIU results:", DIU_results.shape)


# --- differential isoform structure between haplotypes (cells 80/82) ---
def isoform_structure_test(condition, layer, min_similarity):
    adata_cond = poly.filter_low_expressed_genes(adata_isoform[adata_isoform.obs["condition"] == condition, :],
                                                 min_expression=100.0, mode="all", group_col="gene_id", verbose=True)
    adata_cond = poly.calculate_allelic_ratios(adata_cond, "unique_counts")
    adata_cond = poly.calculate_multi_ratios(adata_cond, "unique_counts", "ambiguous_counts")
    adata_cond = poly.calculate_per_allele_ratios(adata_cond, unique_layer="unique_counts", multi_layer="ambiguous_counts")
    poly.add_structure_from_gtf(adata_cond, f"{LR_DIR}/gff/combined.fixed.gtf", inplace=True, verbose=True)
    mask = (adata_cond.var["tx_multimapping_ratio_per_allele_weighted_average"] < 1) & \
           adata_cond.var["synteny_category"].str.contains(SYNTELOG_CATEGORIES, na=False)
    return poly.test_differential_isoform_structure(adata_cond[:, mask].copy(), layer=layer, test_condition="all",
                                                    min_similarity_for_matching=min_similarity, inplace=True)


# parameters as in the source notebook (heat: unique counts, similarity 0.90; control: EM counts, 0.95)
_, structure_plotting_heat = isoform_structure_test("heat", "unique_counts", 0.90)
structure_plotting_heat[structure_plotting_heat["gene_id"].str.startswith(STRUCTURE_GENE)] \
    .to_pickle("isoform_structure_plotting_heat_ABH.pkl")
structure_results_control, _ = isoform_structure_test("control", "em_counts", 0.95)
structure_results_control.to_csv("isoform_structure_results_control.csv", index=False)
print("isoform structure results (control):", structure_results_control.shape)

# --- transcript annotation for the structure plots (cell 43) ---
annotation_gtf = RNApy.read_ensembl_gtf("reference/Desiree_v1/De_v1.unitato_liftoff_haplotap_gffread.gtf")
bambu = RNApy.read_ensembl_gtf(f"{LR_DIR}/gff/combined.fixed.gtf").filter(pl.col("transcript_id").str.contains("Bambu"))
annotation_gtf.extend(bambu)
annotation_gtf = annotation_gtf.with_columns(pl.col("gene_id").str.replace(r"\.[^.]*$", "").alias("gene_id_unitato"))
annotation_gtf.filter(pl.col("gene_id_unitato").is_in(DIU_GENES + [STRUCTURE_GENE])) \
    .write_parquet("annotation_selected_genes.parquet")

# --- transcript expression for the panel D-I genes ---
selected = adata_isoform.var["gene_id_unitato"].isin(DIU_GENES)
subset = adata_isoform[:, selected].copy()
subset.uns = {}
subset.var = subset.var.astype({c: str for c in subset.var.columns if subset.var[c].dtype == object})
subset.write_h5ad("isoform_expression_selected_genes.h5ad")
print("expression subset:", subset.shape)
