from config.settings import T_FOLD_CHANGE
import numpy as np

def fold_change(healthy_samples, cancerous_samples):
    t = np.log2(T_FOLD_CHANGE) # Calculate fc threshold
    fc_degs = [] 
    gene_fc_summary = []


    for i in range(len(healthy_samples)):
        healthy_sample = healthy_samples[i][2::] # Sample without Hugo_Symbol or Entrez_Gene_Id
        cancerous_sample = cancerous_samples[i][2::] # Sample without Hugo_Symbol or Entrez_Gene_Id

        healthy_sample_mean = np.mean(healthy_sample)
        cancerous_sample_mean = np.mean(cancerous_sample)

        if not (healthy_sample_mean == 0 or cancerous_sample_mean == 0):
            fc = cancerous_sample_mean / healthy_sample_mean
            log2fc = np.log2(fc)

            gene_fc_summary.append({"Hugo_Symbol":cancerous_samples[i][0], "Log2FC":log2fc})

            if log2fc > t or log2fc < -t: 
                fc_degs.append({"Hugo_Symbol":cancerous_samples[i][0], "Log2FC":log2fc})

    return [t, fc_degs, gene_fc_summary]
