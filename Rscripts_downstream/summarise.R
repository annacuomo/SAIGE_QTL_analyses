library(data.table)

# mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_highly_expressed_genes/"
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_genes_expressed_in_more_than_1pct_cells/"
myfiles = list.files(mydir)

for (chrom in 1:22){
    chrom = paste0("chr",chrom,"_")
    df_list = list()
    for (file in myfiles){
        gene = gsub("_.*","",file)
        if (length(file[grep(chrom, file)]) == 0){next}
        if (file == "PRR4_B_IN_count_uncond_chr12_cis"){next}
        if (file == "ABCD4_B_IN_count_uncond_chr14_cis"){next}
        if (file == "C15orf40_B_IN_count_uncond_chr15_cis"){next}
        if (file == "HSPA13_B_IN_count_uncond_chr21_cis"){next}
        df_current = read.csv(paste0(mydir, file), sep="\t")
        if (nrow(df_current)==0){
            print(gene)
            next
        }
    #     add info
        df_current$gene = gene
        # combine results across all genes
        df_list[[gene]] = df_current
    }
    df0 = rbindlist(df_list)

    newfile = paste0(mydir,chrom,"1Mb_all_results.txt")
    write.table(df0, newfile, sep="\t", row.names=F, col.names=T, quote=F)
}
