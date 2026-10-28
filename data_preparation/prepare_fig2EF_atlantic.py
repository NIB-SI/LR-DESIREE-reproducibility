# Prepare inputs for Fig. 2E/F (Désirée vs. Atlantic allele-specific expression at shared promoters)
# in clean_notebooks/Fig2_ASE_SV_promoter_analysis_main.ipynb
#
# Follows ANALYSIS/notebooks/0_notebooks_for_paper/7_Desiree_Atlantic_ASE_comparison.ipynb:
#   - haplotype pairs tested for allelic imbalance (polyase pairwise cis test) in both cultivars
#   - whose promoters (1,500 bp) are identical between one Désirée and one Atlantic haplotype for both alleles (BLAST)
#   - pairs with infinite log2FC (allele not expressed) removed; chr08 genes removed
#
# Outputs:
#   atlantic_desiree_shared_promoter_pairs.csv  one row per gene: Désirée / Atlantic allele i and j,
#                                               significant_in_both (pair significant in both cultivars; panel F candidates)
#   atlantic_desiree_allele_cpm_control.csv     CPM per allele and control sample for those alleles
#
# Usage: python prepare_fig2EF_atlantic.py   (run from this directory, polyase conda env)
import re
import numpy as np
import pandas as pd

DIR = "LR_DESIREE_PAPER"
CMP_DIR = f"{DIR}/ANALYSIS/ASE_results/comparison_DES_ATL"
LR_DIR = "nf-core-plantlongrnaseq"

DES_SAMPLES = ["DE_C_H_13", "DE_C_H_14", "DE_C_H_15", "DE_C_H_16"]
ATL_SAMPLES = ["SRR14993894", "SRR14993895", "SRR14996168"]  # control samples used in the source notebook


def merged_cis_results(kind):
    des = pd.read_csv(f"{CMP_DIR}/desiree_pairwise_cis_results_{kind}.csv")
    atl = pd.read_csv(f"{CMP_DIR}/atlantic_pairwise_cis_results_{kind}.csv")
    for df in (des, atl):
        df["Unitato_gene_id"] = df["gene_id_i"].str.rsplit(".", n=1).str[0]
    return pd.merge(des, atl, on="Unitato_gene_id", suffixes=("_desiree", "_atlantic"))


def shared_promoter_genes(merged, id_filter, require_both=True):
    """Genes whose two tested alleles have identical 1,500 bp promoters in Désirée and Atlantic."""
    blast = pd.read_csv(f"{DIR}/DOWNSTREAM_ANALYSIS/shared_atlantic_desiree/BLAST_promotor_ATL_vs_DES_filtered_all.tsv",
                        sep="\t", header=None,
                        names=["Des_gene_id", "Atl_gene_id", "percent_identity", "alignment_length", "mismatches",
                               "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"])
    # DES.Soltu.DM.01G004650.1.Hap1 -> Soltu.DM.01G004650.Hap1
    for col, prefix in [("Des_gene_id", "DES."), ("Atl_gene_id", "ATL.")]:
        blast[col] = blast[col].str.removeprefix(prefix).map(lambda x: re.sub(r"\.[1-5]\.", ".", x, count=1))

    blast = blast[id_filter(blast)]
    des_in = blast["Des_gene_id"].isin(merged["gene_id_i_desiree"]) | blast["Des_gene_id"].isin(merged["gene_id_j_desiree"])
    atl_in = blast["Atl_gene_id"].isin(merged["gene_id_i_atlantic"]) | blast["Atl_gene_id"].isin(merged["gene_id_j_atlantic"])
    # source notebook: '&' for all pairs, '|' for significant pairs (kept as is)
    blast = blast[(des_in & atl_in) if require_both else (des_in | atl_in)]
    blast = blast.drop_duplicates(subset=["Des_gene_id"]).drop_duplicates(subset=["Atl_gene_id"])
    blast["Unitato_gene_id"] = blast["Des_gene_id"].str.rsplit(".", n=1).str[0]

    # exactly two matched haplotypes per gene, forming a tested pair in both cultivars
    n = blast["Unitato_gene_id"].value_counts()
    blast = blast[blast["Unitato_gene_id"].isin(n.index[n == 2])]
    des_pairs = set(zip(merged["gene_id_i_desiree"], merged["gene_id_j_desiree"])) | \
                set(zip(merged["gene_id_j_desiree"], merged["gene_id_i_desiree"]))
    atl_pairs = set(zip(merged["gene_id_i_atlantic"], merged["gene_id_j_atlantic"])) | \
                set(zip(merged["gene_id_j_atlantic"], merged["gene_id_i_atlantic"]))

    genes = {}
    for gene, grp in blast.groupby("Unitato_gene_id"):
        des, atl = grp["Des_gene_id"].tolist(), grp["Atl_gene_id"].tolist()
        if (des[0], des[1]) in des_pairs and (atl[0], atl[1]) in atl_pairs and "08G" not in gene:
            genes[gene] = {"DES_allele_i": des[0], "DES_allele_j": des[1],
                           "ATL_allele_i": atl[0], "ATL_allele_j": atl[1]}
    return genes


# all tested pairs (panel E): identical promoter = no mismatches over 1,500 bp, alleles shared in both cultivars
all_merged = merged_cis_results("all")
# pairs with an infinite/undefined allelic log2FC in either cultivar are dropped before matching promoters
all_merged = all_merged.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC_desiree", "log2FC_atlantic"])
all_genes = shared_promoter_genes(
    all_merged, lambda b: (b["mismatches"] == 0) & (b["alignment_length"] == 1500))

# significant pairs (panel F, orange points in E)
sig_merged = merged_cis_results("significant")
sig_genes = shared_promoter_genes(
    sig_merged, lambda b: (b["percent_identity"] == 100) & (b["alignment_length"] == 1500), require_both=False)
print("shared-promoter genes:", len(all_genes), "| significant in both:", len(sig_genes), sorted(sig_genes))

pairs = pd.DataFrame.from_dict(all_genes, orient="index").rename_axis("Unitato_gene_id").reset_index()
pairs["significant_in_both"] = pairs["Unitato_gene_id"].isin(sig_genes)
# significant genes not in the all-pairs table (different BLAST filter) keep their own allele order
missing = {g: a for g, a in sig_genes.items() if g not in all_genes}
if missing:
    extra = pd.DataFrame.from_dict(missing, orient="index").rename_axis("Unitato_gene_id").reset_index()
    extra["significant_in_both"] = True
    extra["in_panel_E"] = False
    pairs = pd.concat([pairs.assign(in_panel_E=True), extra], ignore_index=True)
else:
    pairs["in_panel_E"] = True
pairs.to_csv("atlantic_desiree_shared_promoter_pairs.csv", index=False)

# CPM per allele (computed over all samples of each run), control samples only
alleles = set(pairs[["DES_allele_i", "DES_allele_j"]].stack()) | set(pairs[["ATL_allele_i", "ATL_allele_j"]].stack())
cpm_tables = []
for cultivar, run, samples in [("Désirée", "output_desiree_liftoff_phased_CDS2exon", DES_SAMPLES),
                               ("Atlantic", "output_atlantic_liftoff_CDS2exon", ATL_SAMPLES)]:
    counts = pd.read_csv(f"{LR_DIR}/{run}/tximeta/all_samples.gene_counts.tsv", sep="\t", index_col="gene_id")
    counts = counts.drop(columns="gene_name")
    cpm = counts.div(counts.sum()) * 1e6
    cpm = cpm.loc[cpm.index.intersection(alleles), samples]
    cpm_tables.append(cpm.rename_axis("allele_id").reset_index()
                      .melt(id_vars="allele_id", var_name="sample", value_name="CPM").assign(cultivar=cultivar))
pd.concat(cpm_tables)[["cultivar", "allele_id", "sample", "CPM"]].to_csv("atlantic_desiree_allele_cpm_control.csv", index=False)
print("alleles with CPM:", pd.concat(cpm_tables)["allele_id"].nunique(), "of", len(alleles))
