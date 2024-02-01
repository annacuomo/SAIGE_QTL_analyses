# this script splits all results from SAIGE by chromosome
# they are also already split by cell type

library(data.table)

results_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"
out_dir_all = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/ct_chrom_summary_stats/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

for (celltype in celltypes){
  ct_results_dir = paste0(results_dir,"cis_",celltype,"/")
  ct_results_file = paste0(ct_results_dir, celltype, "_all_results_cis.txt.gz")
  df = read.csv(gzfile(ct_results_file))
  for (chrom in 1:22){
    df_chrom_ct = df[df$CHR==chrom,]
    out_dir = paste0(out_dir_all,celltype,"/")
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    out_file = paste0(out_dir,celltype,"_chr",chrom,".txt")
    data.table:::fwrite(df_chrom_ct, out_file)
  }
}
