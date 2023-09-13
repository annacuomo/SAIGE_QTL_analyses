# this script combines joint results between Matrix eQTL and SAIGE-QTL
# to create one single comparison plot

# load useful R libraries
library(ggplot2)
library(data.table)

# this sub-folder contains joint tables constructed by the `saige_matrix_qtl_combine.R` script
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_genes_expressed_in_more_than_1pct_cells/joint_tables/"

# combine all results across all chromosomes
df_list = list()
for (chrom in 1:22){
    df0 = read.csv(paste0(mydir,"new_chr",chrom,".tsv"), sep="\t")
    df_list[[chrom]] = df0
}
df = rbindlist(df_list)

# plot p-values

n_tests = nrow(df)
n_genes = length(unique(df$gene))
pv = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['p.value']]
cor = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['estimate']]
text = paste0("all_chroms, ",n_tests," tests, ",n_genes," genes\ncor = ",round(cor,digits=2),", p-value = ",pv)
if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
p = ggplot(df, aes(x=-log10(p.value), y=-log10(i.p.value))) + geom_point() 
p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
p = p + theme(text = element_text(size=20))
myx = (-log10(min(df$p.value))*0.4)
myy = (-log10(min(df$i.p.value))*0.85)
p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
p

# plot betas (effect sizes, noting that the models are fundamentally different)

pv = cor.test(df$BETA, df$beta)[['p.value']]
cor = cor.test(df$BETA, df$beta)[['estimate']]
text = paste0("all_chroms, ",n_tests," tests, ",n_genes," genes\ncor = ",round(cor,digits=2),", p-value = ",pv)
if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
p = ggplot(df, aes(x=BETA, y=beta)) + geom_point() 
p = p + geom_vline(xintercept = 0, col = "firebrick") + geom_hline(yintercept = 0, col = "firebrick") 
p = p + theme_classic() + theme(text = element_text(size=20))
xrange = max(df$BETA)+abs(min(df$BETA))
myx = xrange*0.4 - abs(min(df$BETA))
yrange = max(df$beta)+abs(min(df$beta))
myy = yrange*0.5
p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
p
