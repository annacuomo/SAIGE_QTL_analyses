library(data.table)

args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

# function to obtain weights from distances
# using the dTSS approach described in the APEX paper
distance_to_weight <- function(d, gamma=1e-5){
    w = exp(-gamma*abs(d))
    return(w)
}

# file with gene locations
gene_file = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/GeneLocations.tsv"
gene_df = fread(gene_file)

# group file directory
gf_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/from_wei/gene_group_file_onek1k_forAnna/"
group_files = list.files(gf_dir)

# group file considered here
gf = group_files[i]

# new output directory (groups with weights)
out_dir = paste0(mydir,"with_dTSS_weights/")
out_file = paste0(out_dir, gsub(".grp","_dTSS_weights.grp",gf))
if (file.exists(out_file)){
  quit(save = "no", status = 1, runLast = FALSE)
}

# open file
df = as.data.frame(fread(paste0(mydir,gf), header=F))

# add one more line with weights
# first col is always the gene name
df[3,1] = df[1,1]
# second is description (weight def)
df[3,2] = "weight:dTSS"
gene = df$V1[1]
    if (nrow(gene_df[gene_df$gene_name == gene,])==0){
    print(gene)
    quit(save = "no", status = 1, runLast = FALSE)
}
TSS = gene_df[gene_df$gene_name == gene,"start"]
for (i in c(3:ncol(df))){
    snp = df[1,i]
    pos = as.numeric(unlist(strsplit(snp,split=":"))[2])
    distance = pos-TSS
    weight = distance_to_weight(distance)
    df[3,i] = weight
}
colnames(df) <- c()
fwrite(df, out_file)
