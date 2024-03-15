library(data.table)

# Matrix eQTL results from original paper
matrix_dir = "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/shared_data/OneK1K_matrix_eQTL_results/"
matrix_files = list.files(matrix_dir)

# file mapping rsids to snp pos ids
snps = readRDS("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/hrc_ids_all.rds")
snps = as.data.frame(snps)
colnames(snps)[2] = "SNP"
# remove weird '.' SNPs
snps = snps[snps$SNP != ".",]

out_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/matrix_eqtl_results_snp_info/"
for (file in matrix_files){
    df = as.data.frame(fread(paste0(matrix_dir,file)))
    df_matrix_snps = as.data.frame(data.table(snps)[data.table(df), on="SNP", allow.cartesian=TRUE])
    df_matrix_snps$MarkerID = gsub("_.*","",df_matrix_snps$snpid)
    out_file = paste0(out_dir,gsub(".tsv","_snps.tsv",file))
    fwrite(df_matrix_snps, out_file)
}
