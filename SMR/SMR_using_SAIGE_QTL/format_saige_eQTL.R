####
args = commandArgs(trailingOnly=TRUE)
library(data.table)
library(dplyr)

ct = as.numeric(args[1])

new_names <- c("B_IN", "B_Mem", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")

new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff", "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma")

system(paste0("mkdir -p ", new_names2[ct]))

for(i in 1:22){

ref = fread(paste0("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/SMR_analysis/onek1kdata/plink_renames/names_chr",i,"_201007.tsv"),header=F)
ref = as.data.frame(ref)

res = fread(paste0("/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/ct_chrom_summary_stats/Matrix_eQTL_like/",new_names[ct],"/",new_names[ct],"_chr",i,".tsv"),header=T)
res=as.data.frame(res)

# How many unmatched SNPs
print(nrow(res) - sum(res$SNP %in% ref$V1))

# Remove dot in rsID in ref
ref = ref[!grepl("\\.",ref$V2),]

# Match SNPs
snp = Reduce(intersect,list(ref$V1,res$SNP))
ref = ref[ref$V1 %in% snp,]
res = res[res$SNP %in% snp,]

# Left join
colnames(ref)=c("SNP","rsID")
res = left_join(res,ref,by="SNP")
res$SNP = res$rsID
res = res[,-ncol(res)]

# Recover the dash from dot
replace_dots <- function(gene) {
  dots <- gregexpr("\\.", gene)[[1]]  # Find positions of dots

  if (length(dots) == 1) {
    gene <- gsub("\\.", "-", gene)  # Replace dot with dash
  } else if (length(dots) > 1) {
    gene <- sub("\\.", "-", gene)  # Replace the first dot with dash
  }

  return(gene)
}

res$gene = sapply(res$gene, replace_dots)

# Check
epi=read.table(paste0("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/smr_data/", new_names2[ct],"/", new_names2[ct],"_chr",i,"_210211.epi"),header=F)

res[!res$gene %in% epi$V2,"gene"] = gsub("-", "\\.", res[!res$gene %in% epi$V2,"gene"])

# Special case
# res$gene = gsub("hsa.mir.8072", "hsa-mir-8072", res$gene)
res$gene <- ifelse(grepl("hsa.mir", res$gene), gsub("\\.", "-", res$gene), res$gene)

write.table(res, paste0("./", new_names2[ct], "/", new_names2[ct], "_chr", i, ".tsv"),row.names = F, col.names = T, quote = F)
# write.table(res, paste0("./B_MEM/B_MEM_chr", i, ".tsv"),row.names = F, col.names = T, quote = F)

print(i)
}



####




