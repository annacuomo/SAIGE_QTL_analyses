# Summary of notebooks here

The notebooks in this folder contain code to plot and count the number of rare variant (RV) eGenes, across all 14 cell types from OneK1K.

* [RV results](RV_results_overview.ipynb): Barplot providing an overview of number of genes with at least one signal from rare variants (FDR<5%) across the 14 cell types.
* [RV results weight comparison](RV_eGenes_weights_comparison.ipynb): Three sets of weights were used for combining rare variant results: equal weights, Beta(1,25), distance from transcription start site. This notebook plots a comparison in performance between these across all cell types.
* [RV results conditioning on CV](RV_results_conditional.ipynb): Calculating and plotting how many of the RV signals are independent from common variant (CV) signals for the same gene, by re-running the analysis after conditioning on the top CV variant for each gene.
* [CV vs RV eGenes notebook](CV_vs_RV_eGenes.ipynb): Scatterplot comparing the number of RV vs CV eGenes across the 14 cell types.
