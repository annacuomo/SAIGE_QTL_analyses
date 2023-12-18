# SMR analysis for SAIGE-QTL paper

In order to compare SAIGE-QTL to Matrix eQTL (method used in the [original Science publication](https://www.science.org/doi/full/10.1126/science.abf3041)) further, here we want to compare how they perform in terms of identify Mendelian randomisation signals, using [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview).

## Overview

Like in the original paper, we want to test for whether our identified cell type-specific eQTLs are mediating SNP-disease associations with 7 autoimmune conditions, using GWAS summary statistics for: systemic lupus erythematosus (SLE), rheumatoid arthritis (RA), Crohn’s disease (CD), inflammatory bowel disease (IBD), multiple sclerosis (MS), ankylosing spondylitis (AS), and Type 1 diabetes mellitus (T1D), see Table 1 from Yazar et al.

Here, we run SMR again using summary stats from the same 7 GWAS results, for both Matrix eQTL and SAIGE-QTL, for comparison.

## Detailed workflow

### SAIGE-QTL

Preprocessing:

* [build chromosome-specific summary stats](chr_specific_saige_qtl_results.R)
* [make Matrix-eQTL like](make_matrix_eqtl_like_summary_stats.R)

SMR:

* eQTL p-value threshold: nominal p-value < ?
* no GWAS p-value threshold

### Matrix eQTL

* eQTL p-value threshold: nominal p-value < ?
* no GWAS p-value threshold
