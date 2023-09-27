library(data.table)
library(ggplot2)
library(qvalue)

# load files
saigeqtl_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/cis_Plasma-B_IN-CD4_NC/"
tensorqtl_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/output/tensorqtl/"
fig_dir <- paste0(tensorqtl_dir,"figures/")

celltypes = c("Plasma", "B_IN", "CD4_NC")

for (celltype in celltypes){
  # open SAIGE_QTL file
  saigeqtl_file = paste0(mydir, celltype, "_gene_acat_summary.csv")
  df_saigeqtl = read.csv(saigeqtl_file, row.names=1)
  # re-define gene to match tensorQTL results
  df_saigeqtl$phenotype_id = df_saigeqtl$gene
  # calculate FDR
  df_saigeqtl$qv = qvalue(df_saigeqtl$p.value.cct)$qvalue
  # loop over chromosomes
  for (chrom in 1:22){
    # open tensorQTL file
    tensorqtl_file = paste0(tensorqtl__dir, celltype,"/",celltype,".sig_cis_qtl_pairs.chr",chrom,".csv")
    df_tensorqtl = read.csv(tensorqtl_file, sep="\t")
    # only consider significant results at FDR<10%
    df_saigeqtl = df_saigeqtl[df_saigeqtl$qv<0.1,]
    # combine results
    df_combined = as.data.frame(data.table(df_saigeqtl)[data.table(df_tensorqtl), on = c("phenotype_id"), nomatch=0])
    # extract smalles p-value for plotting later
    m = min(df1$pval_beta,df1$p.value.cct)
    # to avoid log(0), make m the smallest (R) value
    if (m==0){m=5e-324}
    # extract correlation between p-values to print out
    cor = round(cor(-log10(df_combined$pval_beta+5e-324), -log10(df_combined$p.value.cct+5e-324)), digits=2)
    # start building scatter plot
    p = ggplot(df_combined, aes(x=-log10(pval_beta), y=-log10(p.value.cct))) + geom_point() 
    # add diagonal
    p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
    # add title
    p = p + ggtitle(paste0(celltype,", chr",chrom,", cor=",cor))
    # add axes limits
    p = p + xlim(c(0,-log10(m))) + ylim(c(0,-log10(m)))
    p = p + theme(text = element_text(size=20))
    # save
    pdf(paste0(fig_dir,celltype,"_chr",chrom,"_ACAT_BetaPerm_scatter.pdf"), width=8, height=8)
    p
    dev.off()
  }
}  

