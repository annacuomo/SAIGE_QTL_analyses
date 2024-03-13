library(data.table)

# pheno cov filenames (single-cell profiles by cell type, after SCT normalisation)
mydir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/"

# because it is such an abundant cell type, the CD4_NC file was split in two (with about half the genes in each)
# part1
file = paste0(mydir, "CD4_NC_sc_pheno_cov_part1.tsv")
df_pt1 = fread(file)
# part2
file = paste0(mydir, "CD4_NC_sc_pheno_cov_part2.tsv")
df_pt2 = fread(file)

# extract cell-to-individual map from one of the two files (the rows are identical, and so are these two columns)
df_ind_cells = as.data.frame(df_pt1[,c("individual","barcode")])

# for each individual, record the total number of cells (`n_cells_100`)
n_cells_df = data.frame(table(df0$individual))
colnames(n_cells_df) = c("individual","n_cells_100")

# for each individual, get the rounded up number of cells, 
# when subsetting to 1, 5, 10, 20 and 50% 
pct = c(50,20,10,5,1)
for (p in pct){
    n_cells_df[[paste0('n_cells_',p)]] = round(n_cells_df$n_cells_100 * p / 100, digits=0)
}

# after setting a seed, randomly subset cells given the numbers above
set.seed(101)

individuals = unique(df$individual)

# 50%
df_ind_cells_50 = data.frame()
for (ind in individuals){
    df_ind = df_ind_cells[df_ind_cells$individual == ind,]
    n = n_cells_df[n_cells_df$individual == ind,"n_cells_50"]
    df_new = df_ind[sample(nrow(df_ind), n),]
    if((n==nrow(df_new))==F){print(ind)}
    df_ind_cells_50 = rbind(df_ind_cells_50, df_new)
}

# 20%
df_ind_cells_20 = data.frame()
for (ind in individuals){
    df_ind = df_ind_cells[df_ind_cells$individual == ind,]
    n = n_cells_df[n_cells_df$individual == ind,"n_cells_20"]
    df_new = df_ind[sample(nrow(df_ind), n),]
    if((n==nrow(df_new))==F){print(ind)}
    df_ind_cells_20 = rbind(df_ind_cells_20, df_new)
}

# 10%
df_ind_cells_10 = data.frame()
for (ind in individuals){
    df_ind = df_ind_cells[df_ind_cells$individual == ind,]
    n = n_cells_df[n_cells_df$individual == ind,"n_cells_10"]
    df_new = df_ind[sample(nrow(df_ind), n),]
    if((n==nrow(df_new))==F){print(ind)}
    df_ind_cells_10 = rbind(df_ind_cells_10, df_new)
}

# 5%
df_ind_cells_5 = data.frame()
for (ind in individuals){
    df_ind = df_ind_cells[df_ind_cells$individual == ind,]
    n = n_cells_df[n_cells_df$individual == ind,"n_cells_5"]
    df_new = df_ind[sample(nrow(df_ind), n),]
    if((n==nrow(df_new))==F){print(ind)}
    df_ind_cells_5 = rbind(df_ind_cells_5, df_new)
}

# 1%
df_ind_cells_1 = data.frame()
for (ind in individuals){
    df_ind = df_ind_cells[df_ind_cells$individual == ind,]
    n = n_cells_df[n_cells_df$individual == ind,"n_cells_1"]
    df_new = df_ind[sample(nrow(df_ind), n),]
    if((n==nrow(df_new))==F){print(ind)}
    df_ind_cells_1 = rbind(df_ind_cells_1, df_new)
}

# extract the subsetted barcodes 
# and check that the overall numbers are as expected
all_barcodes = df_ind_cells$barcode
length(all_barcodes)
barcodes_50 = df_ind_cells_50$barcode
length(barcodes_50)*2
barcodes_20 = df_ind_cells_20$barcode
length(barcodes_20)*5
barcodes_10 = df_ind_cells_10$barcode
length(barcodes_10)*10
barcodes_5 = df_ind_cells_5$barcode
length(barcodes_5)*20
barcodes_1 = df_ind_cells_1$barcode
length(barcodes_1)*100

# for each of of the barcode lists just generated, 
# * extract single-cell expression 
# * get correct columns and save files

output_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K_from_ScratchGeneral/OneK1K/saige_eqtl/input/Sept23/SCT/"

# 50%
df_pt1_50 = as.data.frame(df_pt1[df_pt1$barcode %in% barcodes_50,])
df_pt2_50 = as.data.frame(df_pt2[df_pt2$barcode %in% barcodes_50,])
df_pt1_50 = df_pt1_50[,!(colnames(df_pt1_50) %in% c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'))]
df_pt2_50 = df_pt2_50[,!(colnames(df_pt2_50) %in% c('individual','barcode'))]
df_50 = cbind(df_pt1_50, df_pt2_50)
pheno_cov_filename = paste0(output_dir,"CD4_NC_sc_pheno_cov_50pct_subset.tsv")
write.table(df_50, pheno_cov_filename, sep="\t", row.names = F, col.names = T, quote = F)

# 20%
df_pt1_20 = as.data.frame(df_pt1[df_pt1$barcode %in% barcodes_20,])
df_pt2_20 = as.data.frame(df_pt2[df_pt2$barcode %in% barcodes_20,])
df_pt1_20 = df_pt1_20[,!(colnames(df_pt1_20) %in% c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'))]
df_pt2_20 = df_pt2_20[,!(colnames(df_pt2_20) %in% c('individual','barcode'))]
df_20 = cbind(df_pt1_20, df_pt2_20)
pheno_cov_filename = paste0(output_dir,"CD4_NC_sc_pheno_cov_20pct_subset.tsv")
write.table(df_20, pheno_cov_filename, sep="\t", row.names = F, col.names = T, quote = F)

# 10%
df_pt1_10 = as.data.frame(df_pt1[df_pt1$barcode %in% barcodes_10,])
df_pt2_10 = as.data.frame(df_pt2[df_pt2$barcode %in% barcodes_10,])
df_pt1_10 = df_pt1_10[,!(colnames(df_pt1_10) %in% c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'))]
df_pt2_10 = df_pt2_10[,!(colnames(df_pt2_10) %in% c('individual','barcode'))]
df_10 = cbind(df_pt1_10, df_pt2_10)
pheno_cov_filename = paste0(output_dir,"CD4_NC_sc_pheno_cov_10pct_subset.tsv")
write.table(df_10, pheno_cov_filename, sep="\t", row.names = F, col.names = T, quote = F)

# 5%
df_pt1_5 = as.data.frame(df_pt1[df_pt1$barcode %in% barcodes_5,])
df_pt2_5 = as.data.frame(df_pt2[df_pt2$barcode %in% barcodes_5,])
df_pt1_5 = df_pt1_5[,!(colnames(df_pt1_5) %in% c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'))]
df_pt2_5 = df_pt2_5[,!(colnames(df_pt2_5) %in% c('individual','barcode'))]
df_5 = cbind(df_pt1_5, df_pt2_5)
pheno_cov_filename = paste0(output_dir,"CD4_NC_sc_pheno_cov_5pct_subset.tsv")
write.table(df_5, pheno_cov_filename, sep="\t", row.names = F, col.names = T, quote = F)

# 1%
df_pt1_1 = as.data.frame(df_pt1[df_pt1$barcode %in% barcodes_1,])
df_pt2_1 = as.data.frame(df_pt2[df_pt2$barcode %in% barcodes_1,])
df_pt1_1 = df_pt1_1[,!(colnames(df_pt1_1) %in% c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'))]
df_pt2_1 = df_pt2_1[,!(colnames(df_pt2_1) %in% c('individual','barcode'))]
df_1 = cbind(df_pt1_1, df_pt2_1)
pheno_cov_filename = paste0(output_dir,"CD4_NC_sc_pheno_cov_1pct_subset.tsv")
write.table(df_1, pheno_cov_filename, sep="\t", row.names = F, col.names = T, quote = F)
