library(data.table)
library(ggplot2)
library(qvalue)

tensorqtl_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/output/tensorqtl/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

df_list = list()
for (celltype in celltypes){
    chrom_df = data.frame(celltype = celltype, chromosome = 1:22)
    celltype_dir = paste0(tensorqtl_dir, celltype,"/")
    for (chrom in 1:22){
        # significant results (main analysis)
        myfile = paste0(celltype_dir,celltype,".sig_cis_qtl_pairs.chr",chrom,".csv")
        df = read.csv(myfile, sep="\t")
        df$qv = qvalue(df$pval_beta, pi0=1)$qvalue
        chrom_df[chrom_df$celltype == celltype & chrom_df$chromosome == chrom,"n_sig_genes"] = nrow(df[df$qv<0.05,])
        # additional results (conditional analysis)
        myfile2 = paste0(celltype_dir,celltype,".independent_cis_qtl_pairs.chr",chrom,".csv")
        if (!file.exists(myfile2)){next}
        df = read.csv(myfile2, sep="\t")
        if (length(unique(df$rank)) == 1){
            if (nrow(df) == 1){df$qv = df$pval_beta}
            else{df$qv = qvalue(df$pval_beta, pi0=1)$qvalue}
            chrom_df[chrom_df$celltype == celltype & chrom_df$chromosome == chrom,"n_genes_rank1"] = nrow(df[df$qv<0.05,])
        }
        if (length(unique(df$rank)) > 1){
            for (rank in unique(df$rank)){
                df_curr = df[df$rank == rank,]
                if (nrow(df_curr) == 1){
                    df_curr$qv = df_curr$pval_beta
                }
                else{df_curr$qv = qvalue(df_curr$pval_beta, pi0=1)$qvalue}
                chrom_df[chrom_df$celltype == celltype & chrom_df$chromosome == chrom,paste0("n_genes_rank",rank)] = nrow(df_curr[df_curr$qv<0.05,])    
            }
        }
    }
    df_list[[celltype]] = chrom_df
}
summary_df = rbindlist(df_list,fill = TRUE)

for (celltype in celltypes){
    ct_df = data.frame()
    celltype_dir = paste0(tensorqtl_dir, celltype,"/")
    for (chrom in 1:22){
        # significant results (main analysis)
        myfile = paste0(celltype_dir,celltype,".sig_cis_qtl_pairs.chr",chrom,".csv")
        df = read.csv(myfile, sep="\t")
        df$qv = qvalue(df$pval_beta, pi0=1)$qvalue
        chrom_df = data.frame(egene = df[df$qv<0.1,"phenotype_id"], fdr = df[df$qv<0.1,"qv"])
        ct_df = rbind(ct_df, chrom_df)
    }
    write.csv(ct_df, paste0("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/tensorqtl_egenes/",celltype,".csv"))
}

# plot number of eGenes (barplot)
celltype_summary = summary_df[, sum(n_sig_genes),by=list(celltype)]

# reoder cell types
celltype_summary$celltype <- factor(celltype_summary$celltype, 
                                    levels = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK",
                                               "NK_R","Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))

# assign original colours to cell types
df_colours = data.frame(colours = c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265",
                                    "#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C"),
                        celltype = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
                                     "Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))

options(repr.plot.width = 12, repr.plot.height = 8) 
p = ggplot(as.data.frame(celltype_summary), aes(x=celltype, y=V1, fill=celltype)) + geom_bar(stat = "identity")
p = p + xlab("") + ylab("Number of genes at FDR<5%") + ylim(c(0,4200))
p = p + scale_fill_manual(values = df_colours$colours) + theme_classic()
# need to save this 
p = p + theme(text = element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

# plot number of eGenes by number of cells (scatter)

meta = read.csv("/share/ScratchGeneral/anncuo/OneK1K/metadata.csv")

ncells_df = data.frame(celltype_old = unique(meta$cell_type))
for (celltype in unique(meta$cell_type)){
    ncells_df[ncells_df$celltype_old == celltype, "n_cells"] = nrow(meta[meta$cell_type == celltype,])
}

# update cell names
celltype <- c("CD4_NC", "CD4_ET", "CD4_SOX4", 
                  "CD8_ET", "CD8_NC", "CD8_S100B", 
                  "NK", "NK_R", 
                  "Plasma", "B_Mem", "B_IN", 
                  "Mono_C", "Mono_NC", "DC")
celltype_old <- c("CD4+ KLRB1- T cell", "CD4+ KLRB1+ T cell", "CD4+ SOX4+ T cell", 
                  "CD8+ GNLY+ NKG7+ T cell",  "CD8+ LTB+ T cell", "CD8+ S100B+ T cell", 
                  "XCL1- NK","XCL1+ NK", 
                  "IgJ+ B cell", "TCL1A- FCER2- B cell", "TCL1A+ FCER2+ B cell", 
                  "Monocyte CD14+", "Monocyte FCGR3A+", "Dendritic cell")  
cells <- data.frame(celltype, celltype_old)

# combine info
ncells_df2 = as.data.frame(data.table(cells)[data.table(ncells_df), on="celltype_old"])
ncells_df3 = as.data.frame(data.table(celltype_summary)[data.table(ncells_df2), on="celltype"])

# exclude non-PBMCs
ncells_df3 = ncells_df3[!(ncells_df3$celltype_old %in% c("Erythrocytes","Platelets")),]
ncells_df3$n_sig_genes = ncells_df3$V1

ncells_df3$celltype <- factor(ncells_df3$celltype, 
                                    levels = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK",
                                               "NK_R","Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))

cor = round(cor(ncells_df3$n_cells, ncells_df3$n_sig_genes), digits=2)

options(repr.plot.width = 10, repr.plot.height = 8) 
p = ggplot(ncells_df3, aes(x=log10(n_cells), y=n_sig_genes, colour=celltype)) + geom_point(size=10) 
p = p + scale_colour_manual(values = df_colours$colours) + theme_classic()
p = p + xlab("Number of cells (log10)") + ylab("Number of significant eGenes (FDR<5%)")
p = p + theme(text = element_text(size=20)) + ggtitle(paste0("TensorQTL, cor=",cor))
print(p)
