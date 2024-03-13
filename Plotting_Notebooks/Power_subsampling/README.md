# Power analysis by cell subsampling

In order to evaluate how the power to detect single-cell eQTLs is affected by the number of cells available, we considered the largest of the cell types (CD4 NC), and downsampled to different percentages, before re-running SAIGE-QTL.
We also built pseudobulk counts and ran TensorQTL for comparison on the same subsets.

* [subset_CD4_NC_cells.R](../../preprocessing/subset_CD4_NC_cells.R): script to subset CD4 NC single-cell expression to 1, 5, 10, 20 and 50% of the cells (keeping the donor-to-cell ratio constant)
* [make_subset_pseudobulk.R](../../preprocessing/make_subset_pseudobulk.R): script to generated pseudobulk countd from the subsetted single-cell counts, in order to run TensorQTL
