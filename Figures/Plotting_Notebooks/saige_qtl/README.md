# Summary of notebooks here

The notebooks in this folder contain code to plot results of running SAIGE-QTL to map single-cell eQTLs across all 14 cell types from OneK1K, and comparison with performing the same task using [TensorQTL](https://github.com/broadinstitute/tensorqtl/).

* [ChiSquare SAIGE-QTL vs TensorQTL](ChiSquare_comparison_plot.ipynb): Comparison of Chi squared values between SAIGE-QTL and TensorQTL for each cell type.
* [CV conditional rounds](Conditional_results.ipynb): Number of eQTLs identified by each of 5 rounds of conditional analysis for each cell type.
* [ITGA4 forest plot](Forest_Plot.ipynb): Forest plot showing eQTL results (using SAIGE-QTL) for ITGA4 across cell types.
* [Scatter n eGenes vs n cells](Number_of_eGenes_by_number_of_cells.ipynb): Scatterplot of number of cells available for each cell type and number of eGenes identified by SAIGE-QTL for the same cell type.
* [SAIGE-QTL vs TensorQTL eGenes Venn](Venn_diagram.ipynb):
* [SAIGE-QTL vs Matrix eQTL pvs and betas scatter](pvals_and_betas_concordance_plots.ipynb): Scatterplot of negative log p-values and effect sizes (betas) between SAIGE-QTL results and results from the [original OneK1K paper](https://www.science.org/doi/full/10.1126/science.abf3041) using [Matrix eQTL](https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/).
