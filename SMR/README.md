# SMR analysis for SAIGE-QTL paper

In order to compare SAIGE-QTL to Matrix eQTL (the method used in the [original Science publication](https://www.science.org/doi/full/10.1126/science.abf3041)) further, here we want to compare how they perform in terms of identify Mendelian randomisation signals, using [SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview).

## Overview

Like in the original paper, we want to test for whether our identified cell type-specific eQTLs are mediating SNP-disease associations with 7 autoimmune conditions, using GWAS summary statistics for: systemic lupus erythematosus (SLE), rheumatoid arthritis (RA), Crohn’s disease (CD), inflammatory bowel disease (IBD), multiple sclerosis (MS), ankylosing spondylitis (AS), and Type 1 diabetes mellitus (T1D), see Table 1 from Yazar et al.

Here, we rerun SMR using summary stats from the same 7 GWAS results, for both Matrix eQTL and SAIGE-QTL, for comparison.

## Detailed workflow

### SAIGE-QTL

Preprocessing:

* [Build chromosome-specific summary stats](chr_specific_saige_qtl_results.R)
* [Make Matrix-eQTL like](make_matrix_eqtl_like_summary_stats.R) 
* [Update the gene names in the eQTL summary](./SMR_using_SAIGE_QTL/format_saige_eQTL.R)
* [Create BESD format of eQTL summary](./SMR_using_SAIGE_QTL/convert_besd_format.qsub.sh)  (What is [BESD format](https://yanglab.westlake.edu.cn/software/smr/#BESDformat)?)
* [Update the epi and esi files for SMR](./SMR_using_SAIGE_QTL/update_epi_esi.qsub.sh)

SMR:

* eQTL p-value threshold: nominal p-value corresponding to FDR<5%
* no GWAS p-value threshold

[Run SMR using SAIGE-eQTL summary and further QCed GWAS summary](./SMR_using_SAIGE_QTL/main_ibd_smr.qsub.sh )

### Matrix eQTL

* eQTL p-value threshold: nominal p-value corresponding to FDR<5%
* no GWAS p-value threshold

[Re-run SMR using Matrix eQTL and further QCed GWAS summary](./SMR_using_MatrixeQTL/main_ibd_smr.qsub.sh )

