print("Loading dependencies ...")
import pandas as pd
from analysis.fold_change import fold_change
from analysis.ht_fc_comparison import ht_fc_comparison
from analysis.hypothesis_testing import hypothesis_testing
from config.paths import HYPOTHESIS_INDEPENDENT_DEGS_PATH, HYPOTHESIS_PAIRED_DEGS_PATH, DEGS_COMPARISON_PATH, VOLCANO_PLOT_GRAPH_PATH
from bioinfokit import visuz
from config.settings import P_CUTOFF
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

[independent_ht_degs, paired_ht_degs, gene_ht_summary] = hypothesis_testing(healthy_samples, cancerous_samples)

# Save 
pd.DataFrame(independent_ht_degs).to_csv(HYPOTHESIS_INDEPENDENT_DEGS_PATH, index=False) # Save independent_degs as a csv file
pd.DataFrame(paired_ht_degs).to_csv(HYPOTHESIS_PAIRED_DEGS_PATH, index=False) # Save paired_degs as a csv file


# ---------------------------------------------------------
# Fold Change
# ---------------------------------------------------------
print("Fold change ...")

[t, fc_degs, gene_fc_summary] = fold_change(healthy_samples, cancerous_samples)


# ---------------------------------------------------------
#  Hypothesis testing vs. Fold change 
# ---------------------------------------------------------
print("Gene data summary ...")

gene_summary_df = ht_fc_comparison(pd.DataFrame(gene_ht_summary), pd.DataFrame(gene_fc_summary))
pd.DataFrame(gene_summary_df).to_csv(DEGS_COMPARISON_PATH, index=False) # Save gene comparison as a csv file

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
    pv_thr=(P_CUTOFF, P_CUTOFF),
    show=False,
    figtype="png",
    figname=VOLCANO_PLOT_GRAPH_PATH
)

# ---------------------------------------------------------
# Gene Set Enrichment Analysis (GSEA)
# ---------------------------------------------------------
print("GSEA ...")

# # Indexing with Hugo_Symbol
indexed_gene_summary_df = gene_summary_df.set_index('Hugo_Symbol')

indexed_paired_degs_df = indexed_gene_summary_df[indexed_gene_summary_df["Is_DEG_Paired"]]

indexed_paired_degs_df["score"] = np.sign(indexed_gene_summary_df["Log2FC"]) * - np.log10(indexed_gene_summary_df["P_Paired"])

# # Rank all genes by Log2FC for GSEA Preranked
ranked_genes = indexed_paired_degs_df['score'].sort_values(ascending=False)

# Run GSEA Preranked with explicit plot export options
gsea_res = gp.prerank(
    rnk=ranked_genes,
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