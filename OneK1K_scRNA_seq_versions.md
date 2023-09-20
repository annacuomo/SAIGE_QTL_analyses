# Data to use for SAIGE-QTL comparison

* Science paper original results, Spearman correlation + Matrix eQTL [1]
* PEER factor paper results, Matrix eQTL [2]
* vQTL paper results, TensorQTL [3]

## Science paper results

* [Yazar et al, Science 2022](https://www.science.org/doi/full/10.1126/science.abf3041)
* 982 individuals
* 14 cell types
* removed genes expressed in < 10% individuals
* pseudobulk expression (mean)
* sub-optimal PEER factors (n=2)
* Matrix eQTL as well as Spearman correlation
* For multiple testing correction (MTC), [manually implemented local FDR directly on nominal p-values](https://github.com/powellgenomicslab/onek1k_phase1/blob/main/single_cell_cis_eQTL_mapping/round1.run_spearman_rank_test.R#L218-L232), results reported for FDR<5%

## PEER factor results

* [Xue et al, Genome Biology 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02873-5)
* 980 individuals
* 14 cell types
* removed genes expressed in < 10% individuals. 
* pseudobulk expression (mean)
* new improved PEER factors. Number of PF for each cell type is different, which is decided by a local greedy algorithm based on the sensitivity analysis.
* Matrix eQTL results
* also local FDR for MTC (here implemented as part of the [qvalue](https://www.bioconductor.org/packages/release/bioc/html/qvalue.html) package)

## vQTL paper

* Xue et al, in preparation
* 980 individuals
* 14 cell types
* removed genes expressed in < 10% individuals
* psedobulk expression (mean), but only on individual-cell type combinations with >= 5 cells
* new PEER factors,  only on individual-cell type combinations with >= 5 cells (n=10)
* TensorQTL
* Beta approx for MTC (FastQTL, using 10,000 permutations)


Note: pseudobulk calculation (in all cases) should be: 1) SCT normalisation, 2) mean for each donor, gene, cell type, 3) ln(x+1) (base e log)
