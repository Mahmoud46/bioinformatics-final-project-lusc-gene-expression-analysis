import pandas as pd
import numpy as np
from config.settings import FC_CUTOFF, P_CUTOFF
from utils.utils import categorize

def ht_fc_comparison(ht_df, fc_df):
    gene_comparison_df = pd.merge(fc_df, ht_df, on="Hugo_Symbol", how="inner")

    gene_comparison_df['Is_DEG_Paired'] = (gene_comparison_df['P_Paired'] < P_CUTOFF) & (gene_comparison_df['Log2FC'].abs() >= np.log2(FC_CUTOFF))
    gene_comparison_df['Is_DEG_Indep'] = (gene_comparison_df['P_Independent'] < P_CUTOFF) & (gene_comparison_df['Log2FC'].abs() >= np.log2(FC_CUTOFF))


    gene_comparison_df['Alignment_Category'] = gene_comparison_df.apply(categorize, axis=1)
    return gene_comparison_df