# Brief description of scripts included here

* Power analysis subsampling scripts:
  * [R script](subset_CD4_NC_cells.R) for subsetting CD4 NC cells to 1, 5, 10, 20 and 50%.
  * [R script](make_subset_pseudobulk.R) to create equivalent pseudobulk counts.
* Add dTSS (distance from transcription start site) weights to rare variant gene group files:
  * [R script](add_weights_group_files.R), and
  * [qsub runner](add_weights_group_files_runner.qsub).
