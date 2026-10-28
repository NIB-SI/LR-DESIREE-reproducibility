# Prepare the per-locus unique mapping rate inputs for Fig1_short_vs_long_read_comparison.ipynb
#
# Unique mapping rate = unique / (unique + ambiguous) reads of single-isoform genes, summed over all 8 samples
# (control DE_C_H_13-16, heat DE_H_H_17-20).
#
# Outputs:
#   unique_mapping_rate_single_isoform_genes.csv     panel B: genes with > 20 reads, present in both technologies
#   hap_short_unique_mapping_rate_single_isoform.csv panel F bins: short-read, haplotype-resolved, >= 5 reads
#   single_isoform_genes_DM.csv / _hap.csv          genes with one transcript in the long-read tx2gene (panels C/E/F)
#
# Usage: python prepare_fig1_panelB.py   (run from this directory, polyase conda env)
import pandas as pd

LR_DIR = "nf-core-plantlongrnaseq"
SR_DIR = "Desiree_ASE/output"
SAMPLES = ["DE_C_H_13", "DE_C_H_14", "DE_C_H_15", "DE_C_H_16",
           "DE_H_H_17", "DE_H_H_18", "DE_H_H_19", "DE_H_H_20"]

tx2gene = {
    "haploid": pd.read_csv(f"{LR_DIR}/output_desiree_unitato/custom/tx2gene.tsv", sep="\t", index_col=0)["gene_id"],
    "phased": pd.read_csv(f"{LR_DIR}/output_desiree_liftoff_phased/tximeta/tx2gene.tsv", sep="\t", index_col=0)["gene_id"],
}


# single-isoform genes per reference (used to restrict the scatterplots, panels C/E/F)
for reference, t2g in tx2gene.items():
    n_isoforms = t2g.value_counts()
    single = n_isoforms.index[n_isoforms == 1].to_series(name="gene_id")
    single.to_csv(f"single_isoform_genes_{'DM' if reference == 'haploid' else 'hap'}.csv", index=False)
    print(reference, "single-isoform genes:", len(single), "of", len(n_isoforms))


def long_read_counts(run):
    # oarfish ambig_info: tname, unique_reads, ambig_reads
    tables = [pd.read_csv(f"{LR_DIR}/{run}/oarfish/{s}.ambig_info.tsv", sep="\t", index_col=0) for s in SAMPLES]
    return sum(t[["unique_reads", "ambig_reads"]] for t in tables).set_axis(["unique", "ambig"], axis=1)


def short_read_counts(run):
    # salmon aux_info/ambig_info.tsv has no transcript names; rows follow quant.sf
    tables = []
    for s in SAMPLES:
        t = pd.read_csv(f"{SR_DIR}/{run}/salmon/{s}/aux_info/ambig_info.tsv", sep="\t")
        t.index = pd.read_csv(f"{SR_DIR}/{run}/salmon/{s}/quant.sf", sep="\t")["Name"]
        tables.append(t)
    return sum(t[["UniqueCount", "AmbigCount"]] for t in tables).set_axis(["unique", "ambig"], axis=1)


def single_isoform_rate(counts, reference):
    # keep genes with exactly one isoform; gene counts are then the transcript counts
    t2g = tx2gene[reference]
    counts = counts.join(t2g, how="inner")
    counts = counts[counts.groupby("gene_id")["gene_id"].transform("size") == 1].set_index("gene_id")
    counts["value"] = counts["unique"] / (counts["unique"] + counts["ambig"])
    return counts


rates = {
    "long RNA-seq haploid":  single_isoform_rate(long_read_counts("output_desiree_unitato"), "haploid"),
    "long RNA-seq phased":   single_isoform_rate(long_read_counts("output_desiree_liftoff_phased"), "phased"),
    "short RNA-seq haploid": single_isoform_rate(short_read_counts("nf_core_rnaseq_short_read_heat_unitato"), "haploid"),
    "short RNA-seq phased":  single_isoform_rate(short_read_counts("nf_core_rnaseq_short_read_heat_phased_desiree_liftoff"), "phased"),
}

# --- panel B: > 20 reads, present in both technologies, DM-annotated genes only ---
expressed = {g: r[r["unique"] + r["ambig"] > 20] for g, r in rates.items()}
common = (set(expressed["long RNA-seq haploid"].index) & set(expressed["short RNA-seq haploid"].index)) | \
         (set(expressed["long RNA-seq phased"].index) & set(expressed["short RNA-seq phased"].index))
df_long = pd.concat([r[["value"]].assign(group=g, gene_id=r.index) for g, r in expressed.items()], ignore_index=True)
df_long = df_long[df_long["gene_id"].isin(common)]
df_long = df_long[~df_long["gene_id"].str.startswith(("OR", "gene"))].dropna()
df_long[["gene_id", "group", "value"]].to_csv("unique_mapping_rate_single_isoform_genes.csv", index=False)
print(df_long.groupby("group").size())

# --- panel F bins: short-read, haplotype-resolved, >= 5 reads ---
short_hap = rates["short RNA-seq phased"]
short_hap = short_hap[short_hap["unique"] + short_hap["ambig"] >= 5]
short_hap.rename(columns={"unique": "UniqueCount", "ambig": "AmbigCount", "value": "UniqueMappingRatio"}) \
    .rename_axis("gene_id").reset_index() \
    .to_csv("hap_short_unique_mapping_rate_single_isoform.csv", index=False)
print("panel F short-read single-isoform genes:", len(short_hap))
