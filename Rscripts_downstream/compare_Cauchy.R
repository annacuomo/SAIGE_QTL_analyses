library(data.table)
library(ggplot2)

results_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"

celltype = "B_IN"

# all variants
file_all = paste0(results_dir, "cis_",celltype,"/",celltype,"_gene_acat_summary.csv")

# MAF>5% only
file_maf5 = paste0(results_dir, "cis_",celltype,"/",celltype,"_gene_acat_maf5_summary.csv")

df_all = read.csv(file_all, row.names = 1)
df_maf5 = read.csv(file_maf5)

# inner join
df_combined = as.data.frame(data.table(df_all)[data.table(df_maf5), on="gene", nomatch=0])

# remove NAs
df_combined <- df_combined[rowSums(is.na(df_combined)) == 0, ] 

# get correlation
cor = round(cor(-log10(df_combined$p.value.cct), -log10(df_combined$i.p.value.cct)), digits=2)

p = ggplot(df_combined, aes(x=-log10(p.value.cct), y=-log10(i.p.value.cct))) + geom_point() + theme_classic()
p = p + xlab("ALL") + ylab("MAF>5%") + theme(text = element_text(size=20)) 
p + ggtitle(paste0(celltype, ", cor=", cor))
