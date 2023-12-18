library(qvalue)

dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/genes/egenes_saigeqtl_not_tensorqtl/ct_chrom_summary_stats/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

for (celltype in celltypes){
    for (chrom in 1:22){
        file = paste0(dir,celltype,"/",celltype,"_chr",chrom,".txt")
        df = data.table:::fread(file)
        if (nrow(df)==0){next}
        # filter to relevant SNPs
        df_bim <- data.table:::fread(sprintf("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/imputed_data/decompressed/filter_vcf_r08_maf005/plink_chr%s.bim", chrom))
        df = df[df$MarkerID %in% df_bim$V2,]
        # additionally, common SNPs only
        df = df[df$AF_Allele>0.05 & df$AF_Allele<0.95,]
        # filter to relevant genes
        input <- sprintf("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/smr/inputs/gene_files/%s/gene_chr%s.tsv", celltype, chrom)
        df_input <- data.table:::fread(input)
        df = df[df$gene %in% df_input$gene_name,]
        # add FDR
        df$qv = qvalue(df$p.value, pi0=1)$qvalue
        df$SNP = paste0(df$MarkerID,"_",df$Allele2)
        # select relevant columns only
        df1 = df[,c("SNP","gene","BETA","Tstat","p.value","qv")]
        # rename columns to fit Matrix eQTL format
        colnames(df1)[3:6] = c("beta", "t-stat", "p-value", "FDR")
        out_file = paste0(dir, "Matrix_eQTL_like/",celltype,"/",celltype,"_chr",chrom,".tsv")
        data.table:::fwrite(df1, out_file,sep="\t")
    }
}
