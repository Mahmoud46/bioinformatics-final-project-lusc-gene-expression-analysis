print("Loading dependencies ...")
import pandas as pd
from analysis.fold_change import fold_change
from analysis.hypothesis_testing import hypothesis_testing
from config.paths import HYPOTHESIS_INDEPENDENT_DEGS_PATH, HYPOTHESIS_PAIRED_DEGS_PATH
from bioinfokit import visuz
import gseapy as gp
import numpy as np

print("Loading main data ...")

healthy_df = pd.read_csv('./project_main_data/lusc-rsem-fpkm-tcga_paired.txt', delimiter="\t")
cancerous_df = pd.read_csv("./project_main_data/lusc-rsem-fpkm-tcga-t_paired.txt", delimiter="\t")

healthy_samples = healthy_df[healthy_df.columns[0:]].values
cancerous_samples = cancerous_df[cancerous_df.columns[0:]].values

# ---------------------------------------------------------
# Hypothesis Testing (Paired vs. Independent)
# ---------------------------------------------------------
print("Hypothesis testing ...")

[independent_degs, paired_degs, gene_stats_summary] = hypothesis_testing(healthy_samples, cancerous_samples)

# Save 
pd.DataFrame(independent_degs).to_csv(HYPOTHESIS_INDEPENDENT_DEGS_PATH, index=False) # Save independent_degs into a csv file
pd.DataFrame(independent_degs).to_csv(HYPOTHESIS_PAIRED_DEGS_PATH, index=False) # Save paired_degs into a csv file


# ---------------------------------------------------------
# Fold Change
# ---------------------------------------------------------
print("Fold change ...")

[t, fc_degs, gene_fc_summary] = fold_change(healthy_samples, cancerous_samples)


# ---------------------------------------------------------
#  Gene Summary 
# ---------------------------------------------------------
print("Gene data summary ...")

gene_stats_summary_df = pd.DataFrame(gene_stats_summary)
gene_fc_summary_df = pd.DataFrame(gene_fc_summary)

gene_summary_df = pd.merge(gene_fc_summary_df, gene_stats_summary_df, on="Hugo_Symbol", how="inner")

# ---------------------------------------------------------
# Volcano Plot
# ---------------------------------------------------------
print("Volcano plot ...")

visuz.GeneExpression.volcano(
    df=gene_summary_df,
    lfc="Log2FC",
    pv="P_Paired",
    plotlegend=True,
    legendpos="upper right",
    lfc_thr=(t, t),
    pv_thr=(0.05, 0.05),
    show=False,
    figtype="png",
    figname="results/volcano_plot"
)

# ---------------------------------------------------------
# Gene Set Enrichment Analysis (GSEA)
# ---------------------------------------------------------
print("GSEA ...")

# # Indexing with Hugo_Symbol
indexed_gene_summary_df = gene_summary_df.set_index('Hugo_Symbol')

# # Rank all genes by Log2FC for GSEA Preranked
ranked_genes = indexed_gene_summary_df['Log2FC'].sort_values(ascending=False)

# print(ranked_genes.value_counts())

# FORCE STRICT UNIQUENESS (Fixes duplicate warning AND unlocks full PNG exports)
# Adding a micro-gradient (1e-12) guarantees zero identical float values
unique_offsets = np.linspace(1e-12, 1e-15, num=len(ranked_genes))
ranked_genes_clean = (ranked_genes + unique_offsets).sort_values(ascending=False)


# 5. Run GSEA Preranked with explicit plot export options
gsea_res = gp.prerank(
    rnk=ranked_genes_clean,
    gene_sets='MSigDB_Hallmark_2020',
    threads=4,                  # Avoids deprecation warning
    min_size=5,
    max_size=500,
    permutation_num=1000,
    outdir='results/gsea_results',
    format='png',
    graph_num=50,
    verbose=True                # Prints detailed file creation logs
)

print("Analysis Complete!")