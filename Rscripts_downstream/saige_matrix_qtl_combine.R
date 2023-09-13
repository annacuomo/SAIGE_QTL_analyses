library(data.table)

# SAIGE-QTL results, naive B cells, cis results
# directories also already contain results aggregated by chromosomes, chr{1:22}_1Mb_all_results.txt

# can't remember what the criteria was, but this is 1,710 genes
# mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_highly_expressed_genes/"

# including all genes expressed in at least 1% of cells, 8,263 genes
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_genes_expressed_in_more_than_1pct_cells/"

# collect all files
myfiles = list.files(mydir)

# loop over (autosome) chromosomes
for (chr in 1:22){
    print(chr)

    # open chromosome specific SAIGE-QTL results
    file = myfiles[grep(paste0("chr",chr,"_1Mb"), myfiles)]
    df0 = read.csv(paste0(mydir,file), sep="\t")
    
    # open Matrix eQTL results (they are stored by chromosome)
    filename = paste0("/share/ScratchGeneral/anncuo/OneK1K/Matrix_eQTL_results/BIN_chr",chr,"_cis_eqtls_210211_snp_ids.tsv")
    df1 = read.csv(filename, sep="\t")

    # adjust SNP ID to match
    df1$MarkerID = gsub("_.*","", df1$snpid)

    # combine results using data.table library
    df_inner_sc = as.data.frame(data.table(df0)[data.table(df1), on = c("MarkerID","gene"), nomatch=0])

    # save combined file
    newfile = paste0(mydir,"joint_tables/new_chr",chr,".tsv")
    write.table(df_inner_sc, newfile, sep="\t")
}
