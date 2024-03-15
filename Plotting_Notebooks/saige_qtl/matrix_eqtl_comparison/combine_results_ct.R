library(data.table)

celltypes = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R","Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC")

out_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/matrix_saige_comparison/"

#######################
###### SAIGE-QTL ######
#######################

# latest SAIGE-QTL results -- common variants, qv<5%
saige_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/from_wei/Feb24/"
all_results_dir = paste0(saige_dir,"cis_single_qvallt0.05/")

for (ct in celltypes){
    all_files = list.files(all_results_dir, pattern=paste0(ct,"_count"))
    df_list = list()
    for (file in all_files){
        df = as.data.frame(fread(paste0(all_results_dir,file)))
        gene = gsub(paste0(ct,"_count_saigeqtl_"),"",gsub("_cis_single.MAFgt0.05.txt","",file))
        df$gene = gsub("\\.","_",gsub("-","_",gene))
        df$celltype = ct
        df_list[[file]] = df
    }
    df_combine = rbindlist(df_list)
    out_file = paste0(out_dir,"saige_qtl_",ct,".csv")
    fwrite(df_combine, out_file)
}

#########################
###### Matrix eQTL ######
#########################

matrix_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/matrix_eqtl_results_snp_info/""

for (ct in celltypes){
    ct_me = gsub("_","",ct)
    df_list = list()
    for (chrom in 1:22){
        df = as.data.frame(fread(paste0(matrix_dir," ",ct_me,"chr",chrom,"_cis_eqtls_210211_snps.tsv")))
        df$gene = gsub("\\.","_",gsub("-","_",df$gene))
        df_list[[chrom]] = df
    }
    df_combine = rbindlist(df_list)
    out_file = paste0(out_dir,"matrix_eqtl_",ct,".csv")
    fwrite(df_combine, out_file)
}
