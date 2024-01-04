library(data.table, quietly = T)
library(R.utils, quietly = T)

args = commandArgs(trailingOnly=TRUE)

pct <- args[1]
chrNumber <- args[2]

pheno_cov_file = paste0("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/CD4_NC_sc_pheno_cov_",pct,"pct_subset_pseudobulk.tsv")
pheno_cov = read.csv(pheno_cov_file, sep="\t")

# extract pseudobulk matrix only
rownames(pheno_cov) = pheno_cov$sampleid
remove_cols = c('sampleid','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2')
mx_t = pheno_cov[,!(colnames(pheno_cov) %in% remove_cols)]

# transpose to get genes by individuals
mx = t(mx_t)

# remove individuals with all NAs
if(sum(is.na(mx[1, ]))>0){
mx = mx[ ,-which(is.na(mx[1, ]))]
}

# remove genes expressed in < 10% individuals
pi0 = rowSums(mx==0)/ncol(mx)
# Remove genes with pi0 > 0.9
# mx = mx[which(pi0 <= 0.9), ]

# log(x+1) and standardization
mx = apply(mx, 2, function(x)log(x+1))
mx = as.data.frame(t(scale(t(mx))))


# get gene info
gene_loc_file = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/GeneLocations.tsv"
loc = read.csv(gene_loc_file, sep="\t")
loc = loc[loc$seqid == chrNumber,]

# extract overlapping genes
loc$gene_name = gsub("-","\\.",loc$gene_name)
if (nrow(loc)>length(unique(loc$gene_name))){
    loc = loc[-which(duplicated(loc$gene_name)),]
}
common_genes = unique(loc$gene_name[loc$gene_name %in% rownames(mx)]

# select relevant columns and adjust order
loc = loc[,c("seqid","start","end","gene_name")]
# e.g., chr22 instead of 22
loc$seqid = paste0("chr",loc$seqid)    

# subset and order both dataframes
loc = loc[loc$gene_name %in% common_genes,]
loc = loc[order(loc$gene_name),]
mx = mx[common_genes,]
mx = mx[order(rownames(mx)),]

# combine location and expression
data = cbind(loc, mx)
data = data[order(data$start, decreasing = F), ]
colnames(data)[1]="#chr"

# remove NAs
data = na.omit(data)

# create dir
output_dir = paste0("/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/input/Jan24/tensorqtl/CD4_NC_",pct,"pct_subset/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# save and zip
filename = paste0(output_dir,"chr",chrNumber,".bed")
write.table(data, filename, row.names = F, col.names = T, quote = F, sep = "\t")
gzip(filename, destname = paste0(filename,".gz"), overwrite = T)

# now onto the covariates
retain_cols = c('sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2')
covs_t = pheno_cov[,colnames(pheno_cov) %in% retain_cols]

# transpose to get covs by individuals
covs = t(covs_t) 

# write to file
cov_filename = paste0(output_dir,"chr",chrNumber,"_covs.txt")
write.table(covs, cov_filename, sep="\t", row.names = T, col.names = T, quote = F)
####                      
