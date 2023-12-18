library(qvalue)
library(data.table)

mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"
tensorqtl_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/output/tensorqtl/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

gene_loc_file = "/share/ScratchGeneral/anncuo/OneK1K/GeneLocations.tsv"
gene_loc_df = data.table:::fread(gene_loc_file)
gene_loc_df$gene_name = gsub("-","\\.",gene_loc_df$gene_name)

for (celltype in celltypes){
    # saige-qtl
    ct_dir = paste0(mydir, "cis_", celltype,"/")
    ct_acat_file = paste0(ct_dir,celltype,"_gene_acat_summary.csv")
    df = as.data.frame(data.table:::fread(ct_acat_file))
    df$qv = qvalue(df$p.value.cct, pi0=1)$qvalue
    genes = df[df$qv<0.05,"gene"]
    # tensorqtl
    celltype_dir = paste0(tensorqtl_dir, celltype,"/")
    df_list = list()
    for (chr in 1:22){
        myfile = paste0(celltype_dir,celltype,".sig_cis_qtl_pairs.chr",chr,".csv")
        df_list[[chr]] = read.csv(myfile, sep="\t")
    }
    summary_df = rbindlist(df_list,fill = TRUE)
    summary_df$qv = qvalue(summary_df$pval_beta, pi0=1)$qvalue
    genes1 = summary_df[summary_df$qv<0.05,"phenotype_id"]$phenotype_id
    print(c(celltype, length(genes[genes %in% genes1]), length(genes[!(genes %in% genes1)]), length(genes1[!(genes1 %in% genes)])))
    df_new = data.frame(gene_name = genes[!(genes %in% genes1)])
    add gene annos
    df_new = gene_loc_df[df_new, on="gene_name"]  # X[DT, on="x"]  left join  (DT is left, X is right)
    data.table:::fwrite(df_new, paste0("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/genes/egenes_saigeqtl_not_tensorqtl/", celltype,".csv"))

}
