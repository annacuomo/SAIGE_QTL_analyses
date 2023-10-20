results_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

for (celltype in celltypes){
    # open all results
    file = paste0(results_dir,"cis_",celltype,"/",celltype,"_all_results_cis.txt")
    df = data.table:::fread(file)
    # remove rare variants
    df = df[df$AF_Allele>0.05 & df$AF_Allele<0.95,]
    # take top (common) SNP only
    df = df[order(df$p.value),]
    df0 = df[-which(duplicated(df$gene)),]
    # write to file
    file = paste0(results_dir,"cis_",celltype,"/",celltype,"_top_common_snp_cis.txt")
    data.table:::fwrite(df0,file,sep="\t")
}
