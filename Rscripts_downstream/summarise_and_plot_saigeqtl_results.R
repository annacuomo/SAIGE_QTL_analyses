library(data.table)
library(ggplot2)
library(qvalue)

saigeqtl_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"

celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "NK_R", "NK", "Plasma", "Mono_C", "Mono_NC")

df_list = list()
for (celltype in celltypes){
    celltype_dir = paste0(saigeqtl_dir, "cis_", celltype,"/")
    ct_df = data.frame(celltype = celltype)
    myfile = paste0(celltype_dir,celltype,"_gene_acat_summary.csv")
    df = read.csv(myfile)
    df$qv = qvalue(df$p.value.cct, pi0=1)$qvalue
    ct_df[ct_df$celltype == celltype,"n_sig_genes"] = nrow(df[df$qv<0.05,])
    df_list[[celltype]] = ct_df
}
summary_df = rbindlist(df_list,fill = TRUE)


summary_df$celltype <- factor(summary_df$celltype, 
                                    levels = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK",
                                               "NK_R","Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))
df_colours = data.frame(colours = c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265",
                                    "#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C"),
                        celltype = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
                                     "Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))

options(repr.plot.width = 12, repr.plot.height = 8) 
p = ggplot(summary_df, aes(x=celltype, y=n_sig_genes, fill=celltype)) + geom_bar(stat = "identity")
p = p + xlab("") + ylab("Number of genes at FDR<5%")
p = p + scale_fill_manual(values = df_colours$colours) 
p = p + theme_classic() + ylim(c(0,4200))
p + theme(text = element_text(size=20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

for (celltype in celltypes){
    celltype_dir = paste0(saigeqtl_dir, "cis_", celltype,"/")
    myfile = paste0(celltype_dir,celltype,"_gene_acat_summary.csv")
    df = read.csv(myfile)
    df$qv = qvalue(df$p.value.cct, pi0=1)$qvalue
    ct_df = data.frame(egene = df[df$qv<0.1,"gene"], fdr = df[df$qv<0.1,"qv"])
    write.csv(ct_df, paste0("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/egenes/",celltype,".csv"))
}

meta = read.csv("/share/ScratchGeneral/anncuo/OneK1K/metadata.csv")
ncells_df = data.frame(celltype_old = unique(meta$cell_type))
for (celltype in unique(meta$cell_type)){
    ncells_df[ncells_df$celltype_old == celltype, "n_cells"] = nrow(meta[meta$cell_type == celltype,])
}
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

ncells_df2 = as.data.frame(data.table(cells)[data.table(ncells_df), on="celltype_old"])
ncells_df3 = as.data.frame(data.table(summary_df)[data.table(ncells_df2), on="celltype"])

# remove non-PBMC cell types
ncells_df3 = ncells_df3[!(ncells_df3$celltype_old %in% c("Erythrocytes","Platelets")),]

cor(ncells_df3$n_cells, ncells_df3$n_sig_genes)

# change cell type order to match colours
ncells_df3$celltype <- factor(ncells_df3$celltype, 
                                    levels = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK",
                                               "NK_R","Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))


options(repr.plot.width = 10, repr.plot.height = 8) 
p = ggplot(ncells_df3, aes(x=log10(n_cells), y=n_sig_genes, colour=celltype)) + geom_point(size=10) 
p = p + scale_colour_manual(values = df_colours$colours) + theme_classic()
p = p + xlab("Number of cells (log10)") + ylab("Number of significant eGenes (FDR<5%)")
p + theme(text = element_text(size=20)) + ggtitle("SAIGE-QTL")
