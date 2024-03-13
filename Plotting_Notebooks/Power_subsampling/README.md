# Power analysis by cell subsampling

In order to evaluate how the power to detect single-cell eQTLs is affected by the number of cells available, we considered the largest of the cell types (CD4 NC), and downsampled to different percentages, before re-running SAIGE-QTL.
We also built pseudobulk counts and ran TensorQTL for comparison on the same subsets.

* [preprocessing/subset_CD4_NC_cells.R](preprocessing/subset_CD4_NC_cells.R)
