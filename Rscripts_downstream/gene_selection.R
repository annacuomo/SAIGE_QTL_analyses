# load libraries
library(ggplot2)

# load relevant pheno cov filename (specific cell type)

# pheno_cov_filename = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/CD4_NC_sc_pheno_cov_part1.tsv"
# pheno_cov_filename = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/B_IN_sc_pheno_cov.tsv"
pheno_cov_filename = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/input/Sept23/SCT/Plasma_sc_pheno_cov.tsv"
df0 = data.table:::fread(pheno_cov_filename, header = T, stringsAsFactors = FALSE)

# make dataframe
df0 = data.frame(df0)

# get 100 values between 0 and 1 (0.01 step)
pct0s = seq(from = 0, to = 0.99, by = 0.01)

# determine number of genes with X percentage of 0s
n_genes = c()  # gene counter 
for (k in pct0s){
    n = 0
    # exclude non-gene columns (individual, cell, covs)
    for (i in 3:(ncol(df0)-10)){
        # determine number of cells with >0 expression for that gene
        x = sum(df0[,i]>0)
        # if number in relevant interval, e.g. (0,0.01] or (0.85,0.86], add 1
        if (x>k*nrow(df0) & x<=(k+0.01)*nrow(df0)){
            n = n+1
        }
    }
    n_genes = c(n_genes,n)
}

options(repr.plot.width = 12, repr.plot.height = 10) 
p = ggplot(data.frame(pct0s = pct0s, n_genes = n_genes), aes(x=pct0s, y=n_genes)) + geom_bar(stat = "identity")
p = p + xlab("Proportion of 0s") + ylab("Number of genes") + ylim(c(0,1500))
p + theme_classic() + theme(text = element_text(size=30))

# select ten genes for each bin
set.seed(100)
df_list = list()
for (k in df1$pct0s){
#     print(k)
    bin_genes = c()  # bin-specific gene collector
    for (i in 3:(ncol(df0)-10)){
        x = sum(df0[,i]>0)
        if (x>k*nrow(df0) & x<=(k+0.01)*nrow(df0)){
            # if genes in that bracket of expression, add to bin-specific gene collector
            bin_genes = c(bin_genes, colnames(df0)[i])
        }
    }
    # for each bin, select one gene at random
    sel_genes = bin_genes[sample(length(bin_genes), min(10,length(bin_genes)))]
    df_curr = data.frame(pct=k)
    n=length(sel_genes)
    if (n>=10){
        df_curr[paste0('gene',1:10)] = sel_genes
    }
    if (n<10){
        df_curr[paste0('gene',1:n)] = sel_genes
        df_curr[paste0('gene',(n+1):10)] = 'null'
    }
    # add selected genes to dataframe
    df_list[[paste0('pct',k)]] = df_curr
}
df_combine = rbindlist(df_list)

my_genes = sel_genes

for (gene in my_genes){
    df_to_plot = data.frame(gene = df0[,gene])
    p = ggplot(df_to_plot, aes(x=gene)) + geom_histogram(alpha = 0.8, bins=20) 
    prop = sum(df0[,gene]==0)/nrow(df0)
    p = p + xlab(paste0(gene, " pct of 0s: ", round(prop*100, digits=2))) + theme_classic()
    p = p + theme(text = element_text(size=20))
    print(p)
}
