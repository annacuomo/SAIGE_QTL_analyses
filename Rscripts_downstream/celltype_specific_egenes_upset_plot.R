library(qvalue)
library(UpSetR)

results_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"
celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
              "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

listInput <- list()
for (celltype in celltypes){
    file = paste0(results_dir, "cis_",celltype,"/",celltype,"_gene_acat_summary.csv")
    df = read.csv(file, row.names = 1)
    df$qv = qvalue(df$p.value.cct, pi0=1)$qvalue
    genes = df[df$qv<0.05,"gene"]
    listInput[[celltype]] <- genes
}

length(listInput)

# basic plot
options(repr.plot.width = 18, repr.plot.height = 10) 
upset(fromList(listInput), nsets = 14, order.by = "freq")

# add colours (a bit hacky, and some have no unique eGenes)

df_colours = data.frame(colours = c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265",
                                    "#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C"),
                        celltype = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
                                     "Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))

upset(fromList(listInput), nsets = 14, order.by = "freq", main.bar.color = "black", 
      queries = list(list(query = intersects, params = list(celltypes[1]), color=df_colours[df_colours$celltype==celltypes[1],"colours"], active = T),
                     list(query = intersects, params = list(celltypes[2]), color=df_colours[df_colours$celltype==celltypes[2],"colours"], active = T),
                     list(query = intersects, params = list(celltypes[3]), color=df_colours[df_colours$celltype==celltypes[3],"colours"], active = T),
                     list(query = intersects, params = list(celltypes[4]), color=df_colours[df_colours$celltype==celltypes[4],"colours"], active = T),
#                      list(query = intersects, params = list(celltypes[5]), color=df_colours[df_colours$celltype==celltypes[5],"colours"], active = T)#,
                     list(query = intersects, params = list(celltypes[6]), color=df_colours[df_colours$celltype==celltypes[6],"colours"], active = T),
                     list(query = intersects, params = list(celltypes[7]), color=df_colours[df_colours$celltype==celltypes[7],"colours"], active = T),
#                      list(query = intersects, params = list(celltypes[8]), color=df_colours[df_colours$celltype==celltypes[8],"colours"], active = T)#,
#                      list(query = intersects, params = list(celltypes[9]), color=df_colours[df_colours$celltype==celltypes[9],"colours"], active = T)#,
                     list(query = intersects, params = list(celltypes[10]), color=df_colours[df_colours$celltype==celltypes[10],"colours"], active = T),
                     list(query = intersects, params = list(celltypes[11]), color=df_colours[df_colours$celltype==celltypes[11],"colours"], active = T),
#                      list(query = intersects, params = list(celltypes[12]), color=df_colours[df_colours$celltype==celltypes[12],"colours"], active = T)#,
                     list(query = intersects, params = list(celltypes[13]), color=df_colours[df_colours$celltype==celltypes[13],"colours"], active = T)#,
#                      list(query = intersects, params = list(celltypes[14]), color=df_colours[df_colours$celltype==celltypes[14],"colours"], active = T)
                    ))
