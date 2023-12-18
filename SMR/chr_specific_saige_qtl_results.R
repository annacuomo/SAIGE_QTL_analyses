library(data.table)

results_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"
unique_egenes_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/genes/egenes_saigeqtl_not_tensorqtl/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

for (celltype in celltypes){
    ct_results_dir = paste0(results_dir,"cis_",celltype,"/")
    ct_results_files = list.files(ct_results_dir, pattern="singleVar")
    select_genes_file = paste0(unique_egenes_dir,celltype,".csv")
    select_genes_df = as.data.frame(data.table:::fread(select_genes_file))
    for (chrom in 1:22){
        genes = select_genes_df[select_genes_df$seqid == chrom, "gene_name"]
        df_list = list()
        for (file in ct_results_files){
            gene = gsub("_.*","",file)
            if (!(gene %in% genes)){next}
            filename = paste0(ct_results_dir,file)
            df = as.data.frame(data.table:::fread(filename))
            df$gene = gene
            df_list[[file]] = df
        }
        df_chrom_ct = rbindlist(df_list)
        out_file = paste0(unique_egenes_dir,"ct_chrom_summary_stats/",celltype,"/",celltype,"_chr",chrom,".txt")
        data.table:::fwrite(df_chrom_ct, out_file)
    }
}
