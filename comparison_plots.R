library(cowplot)
library(ggplot2)


# Genes expressed in > 1% cells
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_genes_expressed_in_more_than_1pct_cells/joint_tables/"

# plot p-values (SAIGE-QTL vs Matrix eQTL)

for (chrom in 1:22){
    df = read.csv(paste0(mydir, "new_chr", chrom, ".tsv"), sep="\t")
    n_tests = nrow(df)
    n_genes = length(unique(df$gene))
    pv = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['p.value']]
    cor = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['estimate']]
    text = paste0("chrom", chrom, ", ", n_tests, " tests, ", n_genes, " genes\ncor = ", round(cor,digits=2), ", p-value = ", pv)
    if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
    p = ggplot(df, aes(x=-log10(p.value), y=-log10(i.p.value))) + geom_point() 
    p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
    p = p + theme(text = element_text(size=20))
    myx = (-log10(min(df$p.value, na.rm = TRUE))*0.4)
    myy = (-log10(min(df$i.p.value, na.rm = TRUE))*0.85)
    p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
    print(p)
}

# plot betas (SAIGE-QTL vs Matrix eQTL)

for (chrom in 1:22){
    df = read.csv(paste0(mydir,"new_chr",chrom,".tsv"), sep="\t")
    n_tests = nrow(df)
    n_genes = length(unique(df$gene))
    pv = cor.test(df$BETA, df$beta)[['p.value']]
    cor = cor.test(df$BETA, df$beta)[['estimate']]
    text = paste0("chrom", chrom,", ",n_tests," tests, ",n_genes," genes\ncor = ",round(cor,digits=2),", p-value = ",pv)
    if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
    p = ggplot(df, aes(x=BETA, y=beta)) + geom_point() 
    p = p + geom_vline(xintercept = 0, col = "firebrick") + geom_hline(yintercept = 0, col = "firebrick") 
    p = p + theme_classic() + theme(text = element_text(size=20))
    xrange = max(df$BETA)+abs(min(df$BETA))
    myx = xrange*0.4 - abs(min(df$BETA))
    yrange = max(df$beta)+abs(min(df$beta))
    myy = yrange*0.5
    p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
    print(p)
}


# Highly expressed genes only (def?)
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/B_IN_cis_results_highly_expressed_genes/joint_tables/"

# plot p-values (SAIGE-QTL vs Matrix eQTL)

for (chrom in 1:22){
    df = read.csv(paste0(mydir,"chr",chrom,".tsv"), sep="\t")
    n_tests = nrow(df)
    n_genes = length(unique(df$gene))
    pv = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['p.value']]
    cor = cor.test(-log10(df$p.value), -log10(df$i.p.value))[['estimate']]
    text = paste0("chrom", chrom,", ",n_tests," tests, ",n_genes," genes\ncor = ",round(cor,digits=2),", p-value = ",pv)
    if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
    p = ggplot(df, aes(x=-log10(p.value), y=-log10(i.p.value))) + geom_point() 
    p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
    p = p + theme(text = element_text(size=20))
    myx = (-log10(min(df$p.value))*0.4)
    myy = (-log10(min(df$i.p.value))*0.85)
    p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
    print(p)
}

# plot betas (SAIGE-QTL vs Matrix eQTL)

for (chrom in 1:22){
    df = read.csv(paste0(mydir,"chr",chrom,".tsv"), sep="\t")
    n_tests = nrow(df)
    n_genes = length(unique(df$gene))
    pv = cor.test(df$BETA, df$beta)[['p.value']]
    cor = cor.test(df$BETA, df$beta)[['estimate']]
    text = paste0("chrom", chrom,", ",n_tests," tests, ",n_genes," genes\ncor = ",round(cor,digits=2),", p-value = ",pv)
    if (pv==0){text = gsub("p-value = 0","p-value < 2e-16",text)}
    p = ggplot(df, aes(x=BETA, y=beta)) + geom_point() 
    p = p + geom_vline(xintercept = 0, col = "firebrick") + geom_hline(yintercept = 0, col = "firebrick") 
    p = p + theme_classic() + theme(text = element_text(size=20))
    xrange = max(df$BETA)+abs(min(df$BETA))
    myx = xrange*0.4 - abs(min(df$BETA))
    yrange = max(df$beta)+abs(min(df$beta))
    myy = yrange*0.5
    p = p + annotate("text", label = text, x = myx, y = myy, size = 7)
    print(p)
}


# Calibration plots

options(repr.plot.width = 14, repr.plot.height = 8) 

set.seed(12345) 
for (chrom in 1:22){
    df = read.csv(paste0(mydir,"chr",chrom,".tsv"), sep="\t")
    df <- df[rowSums(is.na(df)) == 0, ]
    df$pv_uniform <- runif(dim(df)[1], min = 0, max = 1)
    p1 = ggplot(df, aes(x=sort(-log10(pv_uniform)), y=sort(-log10(p.value)))) + geom_point() 
    p1 = p1 + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
    p1 = p1 + theme(text = element_text(size=20)) + ggtitle(paste0("chrom",chrom,", SAIGE-QTL"))
    p2 = ggplot(df, aes(x=sort(-log10(pv_uniform)), y=sort(-log10(i.p.value)))) + geom_point() 
    p2 = p2 + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
    p2 = p2 + theme(text = element_text(size=20)) + ggtitle(paste0("chrom",chrom,", Matrix eQTL"))
    print(plot_grid(p1, p2, ncol = 2))
}
