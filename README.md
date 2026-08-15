# Differentially Expressed Genes (DEGs) Analysis in Lung Squamous Cell Carcinoma (LUSC)

A data-driven bioinformatics study analyzing paired RNA-seq Gene Expression (GE) data for **Lung Squamous Cell Carcinoma (LUSC)** from TCGA. This project evaluates differential expression under paired vs. independent statistical assumptions, identifies key up- and down-regulated genes, and identifies biological pathways using Gene Set Enrichment Analysis (GSEA).

## Overview

The primary objective of this project is to perform end-to-end differential gene expression (DEG) profiling on LUSC TCGA expression datasets:

1. **Hypothesis Testing:** Evaluate the impact of study design by comparing Paired $t$-tests against Independent Two-Sample $t$-tests.
2. **Fold Change Analysis:** Compute log2 fold changes ($\log_2 \text{FC}$) to measure expression magnitude changes.
3. **Volcano Plot Identification:** Isolate biologically and statistically significant DEGs.
4. **Pathway Enrichment (GSEA):** Run Functional Enrichment Analysis using MSigDB Hallmark Gene Sets via gseapy.

## Dataset Description

The dataset consists of two tab-separated `.txt` FPKM gene expression files located in `project_main_data/:File`
|Name|Description|
|---|---|
|`lusc-rsem-fpkm-tcga-t_paired.txt`|Gene expression levels in **Tumor / Cancer** tissues|
|`lusc-rsem-fpkm-tcga_paired.txt`|Gene expression levels in **Healthy / Control** tissues|

**Dataset Note:** _Samples are paired across both files in identical patient order. Genes are identified by standard_ `Hugo_Symbol` _identifiers._

---

© December 2022 Mahmoud Zakaria, All rights reserved.
