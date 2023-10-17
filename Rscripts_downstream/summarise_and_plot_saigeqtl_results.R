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
