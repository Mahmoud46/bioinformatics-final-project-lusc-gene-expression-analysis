from scipy import stats
import numpy as np

def hypothesis_testing(healthy_samples, cancerous_samples):
    independent_degs = [] # case independent degs array
    paired_degs = [] # case paired degs array
    gene_ht_summary = []

    for i in range(len(healthy_samples)):
        healthy_sample = healthy_samples[i][2::] # Sample without Hugo_Symbol or Entrez_Gene_Id
        cancerous_sample = cancerous_samples[i][2::] # Sample without Hugo_Symbol or Entrez_Gene_Id

        t_independent, p_independent = stats.ttest_ind(np.array(healthy_sample, float), np.array(cancerous_sample, float)) # independent t test
        t_paired, p_paired = stats.ttest_rel(np.array(healthy_sample, float), np.array(cancerous_sample, float)) # paired t test

        # Null hypothesis => Population means are equal (μ1 = μ2)
        # p_value <= 0.05  => Statistically significant difference (Reject Null)

        gene_ht_summary.append({"Hugo_Symbol": cancerous_samples[i][0], "T_Paired": float(t_paired), "P_Paired": float(p_paired), "T_Independent": float(t_independent), "P_Independent": float(p_independent)})
        
        if p_independent <= 0.05:   # Reject Null Hypothesis for independent case
            independent_degs.append({"Hugo_Symbol": cancerous_samples[i][0], "T_Independent": float(t_independent), "P_Independent": p_independent})

        if p_paired <= 0.05:   # Reject Null Hypothesis for paired case
            paired_degs.append({"Hugo_Symbol": cancerous_samples[i][0], "T_Paired": float(t_paired), "P_Paired": p_paired})

    return [independent_degs, paired_degs, gene_ht_summary]