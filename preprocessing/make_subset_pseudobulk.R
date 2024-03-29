library(data.table)

# donor-level covariates for the CD4 NC cell type
peer_dir <- "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/flagship_paper/"
celltype = 'CD4all'
base.dir <- paste0(peer_dir, celltype, "/")
covs <- read.table(paste0(base.dir, celltype, "_peer_factors.tsv", collapse = ''), header = T)

# only the first two PEER factors were used in the original script and throughout our analyses too
cols = colnames(covs)[!(colnames(covs) %in% c(paste0("pf",3:10)))]
covs_sel = covs[,cols]

## for 1,5,10,20,50%
pcts = c(1,5,10,20,50)
remove_cols = c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2','barcode')
for (pct in pcts){
    pheno_cov_filename_sc = paste0(sct_dir,"CD4_NC_sc_pheno_cov_",pct,"pct_subset.tsv")
    pheno_cov_filename_pb = gsub("subset.tsv","subset_pseudobulk.tsv",pheno_cov_filename_sc)

    df = as.data.frame(fread(pheno_cov_filename_sc))
    rownames(df) = df$barcode
    
    df1 = df[,!(colnames(df) %in% remove_cols)]
    
    donors = unique(df1$individual)
    
    mat = matrix(0,nrow=length(donors),ncol=(ncol(df1)-1))
    rownames(mat) = donors
    for (donor in donors){
        df_curr = df1[df1$individual == donor,]
        df_curr$individual = c()
        mat[donor,] = colMeans(df_curr)
    }
    colnames(mat) = colnames(df1)[2:ncol(df1)]
    df_mean = as.data.frame(mat)
    
    df_mean$sampleid = donors
    pheno_cov_df = merge(df_mean, covs_sel, by.y = "sampleid")
    write.table(pheno_cov_df, pheno_cov_filename_pb, sep="\t", row.names = F, col.names = T, quote = F)
}

## do the same but for all 100% cells

file = paste0(sct_dir, "CD4_NC_sc_pheno_cov_part1.tsv")
df_pt1 = fread(file)

file = paste0(sct_dir, "CD4_NC_sc_pheno_cov_part2.tsv")
df_pt2 = fread(file)

pheno_cov_filename_pb = gsub("part2.tsv","pseudobulk.tsv",file)
donors = unique(df_pt2$individual)

df_pt1 = as.data.frame(df_pt1)
rownames(df_pt1) = df_pt1$barcode
remove_cols = c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2','barcode')
df1 = df_pt1[,!(colnames(df_pt1) %in% remove_cols)]

mat1 = matrix(0,nrow=length(donors),ncol=(ncol(df1)-1))
rownames(mat1) = donors
for (donor in donors){
    df_curr = df1[df1$individual == donor,]
    df_curr$individual = c()
    mat1[donor,] = colMeans(df_curr)
}
colnames(mat1) = colnames(df1)[2:ncol(df1)]
df_mean1 = as.data.frame(mat1)

df_pt2 = as.data.frame(df_pt2)
rownames(df_pt2) = df_pt2$barcode
remove_cols = c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2','barcode')
df2 = df_pt2[,!(colnames(df_pt2) %in% remove_cols)]

mat2 = matrix(0,nrow=length(donors),ncol=(ncol(df2)-1))
rownames(mat2) = donors
for (donor in donors){
    df_curr = df2[df2$individual == donor,]
    df_curr$individual = c()
    mat2[donor,] = colMeans(df_curr)
}

colnames(mat2) = colnames(df2)[2:ncol(df2)]
df_mean2 = as.data.frame(mat2)

df_mean = cbind(df_mean1, df_mean2)
df_mean$sampleid = donors
pheno_cov_df = merge(df_mean, covs_sel, by.y = "sampleid")
write.table(pheno_cov_df, pheno_cov_filename_pb, sep="\t", row.names = F, col.names = T, quote = F)
