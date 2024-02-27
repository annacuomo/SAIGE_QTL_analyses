library(data.table)

saige_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/from_wei/Feb24/"
results_dir = paste0(saige_dir,"cis_single_qvallt0.05/")

myfile = paste0(saige_dir, "allcelltype_CauchyPval_cis_MAFge0.05.txt.qvalue.lt0.05.txt.05")
df = as.data.frame(fread(myfile))
colnames(df) = c("celltype","gene","cauchy_pval","qval")

# unique cell types
celltypes = unique(df$celltype)

# loop over cell types
for (celltype in celltypes){
    df1 = as.data.frame(df[df$celltype == celltype,])
    genes = df1[df1$qval==max(df1[df1$qval<0.05,]$qval),]$gene
    max_pv = 0
    for (gene in genes){
        file = paste0(results_dir,celltype,"_count_saigeqtl_",gene,"_cis_single.MAFgt0.05.txt")
        df2 = fread(file)
        max_pv = max(max_pv, min(df2$p.value))
    }
    print(c(celltype, max_pv))
}
