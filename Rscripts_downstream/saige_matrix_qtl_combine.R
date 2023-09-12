library(data.table)

# mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_highly_expressed_genes/"
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_genes_expressed_in_more_than_1pct_cells/"
myfiles = list.files(mydir)

for (chr in 1:22){
    print(chr)
    file = myfiles[grep(paste0("chr",chr,"_1M"), myfiles)]
    df0 = read.csv(paste0(mydir,file), sep="\t")
    filename = paste0("/share/ScratchGeneral/anncuo/OneK1K/Matrix_eQTL_results/BIN_chr",chr,"_cis_eqtls_210211_snp_ids.tsv")
    df1 = read.csv(filename, sep="\t")
    df1$MarkerID = gsub("_.*","", df1$snpid)
    df_inner_sc = as.data.frame(data.table(df0)[data.table(df1), on = c("MarkerID","gene"), nomatch=0])
    newfile = paste0(mydir,"joint_tables/new_chr",chr,".tsv")
    write.table(df_inner_sc, newfile, sep="\t")
}
