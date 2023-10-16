# run as: Rscript step1_summarise_by_celltype.R {celltype}
# celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
#               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")


library(data.table)

mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"

# change, have this as an argument
args <- commandArgs(trailingOnly=TRUE)
celltype <- args[1]

ct_dir = paste0(mydir, "cis_", celltype,"/")
ct_files = list.files(ct_dir, pattern = "singleVar.txt")

# combine results across genes
file_df_list = list()
for (file in ct_files){
    filename = paste0(ct_dir,file)
    df_curr = as.data.frame(data.table:::fread(filename, header = T, stringsAsFactors = FALSE))
    if (nrow(df_curr) == 0){next}
    df_curr[["gene"]] = gsub("_.*","",file)
    file_df_list[[file]] = df_curr
}
df = rbindlist(file_df_list, fill=TRUE)

# save all results
new_file_all = paste0(ct_dir,celltype,"_all_results_cis.txt")
data.table:::fwrite(df,new_file_all)

# save top SNP per gene only
df = df[order(df$p.value),]
df0 = df[-which(duplicated(df$gene)),]
new_file_leads = paste0(ct_dir,celltype,"_top_snp_cis.txt")
data.table:::fwrite(df0,new_file_leads)
