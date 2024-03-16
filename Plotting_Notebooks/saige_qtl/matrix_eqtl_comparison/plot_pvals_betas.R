library(data.table)
library(ggplot2)

# assign original colours to cell types
df_colours = data.frame(colours = c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265",
                                    "#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C"),
                        celltype = c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
                                     "Plasma","B_Mem","B_IN","Mono_C","Mono_NC","DC"))
celltypes = df_colours$celltype

files_dir = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/matrix_saige_comparison/"

# one SNP per gene (based on matrix results)
fig_dir <- "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/ms_figures/matrix_saige_comparison/top_snp_matrix/"
options(repr.plot.width = 8, repr.plot.height = 8)
cors_pvs = c()
cors_betas = c()
cors_abs_betas = c()
for (ct in celltypes){
    file_matrix = list.files(files_dir, pattern = paste0("matrix_eqtl_",ct))
    file_saige = list.files(files_dir, pattern = paste0("saige_qtl_",ct))
    # load Matrix eQTL results
    df_m = fread(paste0(files_dir,file_matrix))
    colnames(df_m) = gsub("-","_",colnames(df_m))
    # select top SNP only (Matrix eQTL)
    df_m = df_m[order(df_m$p_value),]
    df_m = df_m[-which(duplicated(df_m$gene)),]
    # load SAIGE-QTL results
    df_s = fread(paste0(files_dir,file_saige))
    # combine
    df_both = df_m[df_s, on=c("gene","MarkerID"), nomatch=0]
    df_both = as.data.frame(df_both)
    # change col names
    df_both$pval_matrix = df_both$p_value
    df_both$pval_saige = df_both$p.value
    df_both$beta_matrix = df_both$beta
    df_both$beta_saige = df_both$BETA
    # establish plot colour based on cell type
    col = df_colours[df_colours$celltype == ct,"colours"]
    # plot p-values
    m = -log10(min(min(df_both$pval_matrix),min(df_both$pval_saige)))
    p = ggplot(df_both, aes(x=-log10(pval_matrix), y=-log10(pval_saige))) + geom_point(col=col, size=2.5) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic()
    p = p + theme(text = element_text(size=20)) + xlim(c(0,m)) + ylim(c(0,m))
    p = p + xlab("-log10 p-value (Matrix eQTL)") + ylab("-log10 p-value (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_pvals_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    ggsave(paste0(fig_dir,ct,"_pvals_scatter.pdf"), p, width = 8, height = 8, dpi=1500)
    # plot betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=beta_matrix, y=beta_saige)) + geom_point(col=col, size=2) 
    p = p + geom_vline(xintercept = 0, col = "gray") + geom_hline(yintercept = 0, col = "gray") 
    p = p + theme_classic() + theme(text = element_text(size=20))
    p = p + xlim(c(-M_matrix,M_matrix)) + ylim(c(-M_saige,M_saige))
    p = p + xlab("effect size (Matrix eQTL)") + ylab("effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # plot abs betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=abs(beta_matrix), y=abs(beta_saige))) + geom_point(col=col, size=2) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic() 
    p = p + theme(text = element_text(size=20)) + xlim(c(0,M_matrix)) + ylim(c(0,M_saige))
    p = p + xlab("absolute effect size (Matrix eQTL)") + ylab("absolute effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_abs_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # print correlations
    cor_pvs = cor(-log10(df_both$pval_matrix), -log10(df_both$pval_saige))
    cor_betas = cor(df_both$beta_matrix, df_both$beta_saige)
    cor_abs_betas = cor(abs(df_both$beta_matrix), abs(df_both$beta_saige))
    print(paste0(ct,", cor pvs: ",round(cor_pvs,digits=2),", cor betas: ",round(cor_betas,digits=2)))
    print(paste0(ct,", cor abs betas: ",round(cor_betas,digits=2)))
    cors_pvs = c(cors_pvs, cor_pvs)
    cors_betas = c(cors_betas, cor_betas)
    cors_abs_betas = c(cors_abs_betas, cor_abs_betas)
}

# one SNP per gene (based on saige results)
fig_dir <- "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/ms_figures/matrix_saige_comparison/top_snp_saige/"
options(repr.plot.width = 8, repr.plot.height = 8)
cors_pvs = c()
cors_betas = c()
cors_abs_betas = c()
for (ct in celltypes){
    file_matrix = list.files(files_dir, pattern = paste0("matrix_eqtl_",ct))
    file_saige = list.files(files_dir, pattern = paste0("saige_qtl_",ct))
    # load Matrix eQTL results
    df_m = fread(paste0(files_dir,file_matrix))
    colnames(df_m) = gsub("-","_",colnames(df_m))
    # load SAIGE-QTL results
    df_s = fread(paste0(files_dir,file_saige))
    # select top SNP only (SAIGE-QTL)
    df_s = df_s[order(df_s$p.value),]
    df_s = df_s[-which(duplicated(df_s$gene)),]
    # combine
    df_both = df_m[df_s, on=c("gene","MarkerID"), nomatch=0]
    df_both = as.data.frame(df_both)
    # change col names
    df_both$pval_matrix = df_both$p_value
    df_both$pval_saige = df_both$p.value
    df_both$beta_matrix = df_both$beta
    df_both$beta_saige = df_both$BETA
    # establish plot colour based on cell type
    col = df_colours[df_colours$celltype == ct,"colours"]
    # plot p-values
    m = -log10(min(min(df_both$pval_matrix),min(df_both$pval_saige)))
    p = ggplot(df_both, aes(x=-log10(pval_matrix), y=-log10(pval_saige))) + geom_point(col=col, size=2.5) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic()
    p = p + theme(text = element_text(size=20)) + xlim(c(0,m)) + ylim(c(0,m))
    p = p + xlab("-log10 p-value (Matrix eQTL)") + ylab("-log10 p-value (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_pvals_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    ggsave(paste0(fig_dir,ct,"_pvals_scatter.pdf"), p, width = 8, height = 8, dpi=1500)
    # plot betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=beta_matrix, y=beta_saige)) + geom_point(col=col, size=2) 
    p = p + geom_vline(xintercept = 0, col = "gray") + geom_hline(yintercept = 0, col = "gray") 
    p = p + theme_classic() + theme(text = element_text(size=20))
    p = p + xlim(c(-M_matrix,M_matrix)) + ylim(c(-M_saige,M_saige))
    p = p + xlab("effect size (Matrix eQTL)") + ylab("effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # plot abs betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=abs(beta_matrix), y=abs(beta_saige))) + geom_point(col=col, size=2) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic() 
    p = p + theme(text = element_text(size=20)) + xlim(c(0,M_matrix)) + ylim(c(0,M_saige))
    p = p + xlab("absolute effect size (Matrix eQTL)") + ylab("absolute effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_abs_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # print correlations
    cor_pvs = cor(-log10(df_both$pval_matrix), -log10(df_both$pval_saige))
    cor_betas = cor(df_both$beta_matrix, df_both$beta_saige)
    cor_abs_betas = cor(abs(df_both$beta_matrix), abs(df_both$beta_saige))
    print(paste0(ct,", cor pvs: ",round(cor_pvs,digits=2),", cor betas: ",round(cor_betas,digits=2)))
    print(paste0(ct,", cor abs betas: ",round(cor_betas,digits=2)))
    cors_pvs = c(cors_pvs, cor_pvs)
    cors_betas = c(cors_betas, cor_betas)
    cors_abs_betas = c(cors_abs_betas, cor_abs_betas)
}

# all SNPs
fig_dir <- "/directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/ms_figures/matrix_saige_comparison/"
options(repr.plot.width = 8, repr.plot.height = 8)
cors_pvs = c()
cors_betas = c()
cors_abs_betas = c()
for (ct in celltypes){
    file_matrix = list.files(files_dir, pattern = paste0("matrix_eqtl_",ct))
    file_saige = list.files(files_dir, pattern = paste0("saige_qtl_",ct))
    # load Matrix eQTL results
    df_m = fread(paste0(files_dir,file_matrix))
    colnames(df_m) = gsub("-","_",colnames(df_m))
    # load SAIGE-QTL results
    df_s = fread(paste0(files_dir,file_saige))
    # combine
    df_both = df_m[df_s, on=c("gene","MarkerID"), nomatch=0]
    df_both = as.data.frame(df_both)
    # change col names
    df_both$pval_matrix = df_both$p_value
    df_both$pval_saige = df_both$p.value
    df_both$beta_matrix = df_both$beta
    df_both$beta_saige = df_both$BETA
    # establish plot colour based on cell type
    col = df_colours[df_colours$celltype == ct,"colours"]
    # plot p-values
    m = -log10(min(min(df_both$pval_matrix),min(df_both$pval_saige)))
    p = ggplot(df_both, aes(x=-log10(pval_matrix), y=-log10(pval_saige))) + geom_point(col=col, size=2.5) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic()
    p = p + theme(text = element_text(size=20)) + xlim(c(0,m)) + ylim(c(0,m))
    p = p + xlab("-log10 p-value (Matrix eQTL)") + ylab("-log10 p-value (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_pvals_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    ggsave(paste0(fig_dir,ct,"_pvals_scatter.pdf"), p, width = 8, height = 8, dpi=1500)
    # plot betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=beta_matrix, y=beta_saige)) + geom_point(col=col, size=2) 
    p = p + geom_vline(xintercept = 0, col = "gray") + geom_hline(yintercept = 0, col = "gray") 
    p = p + theme_classic() + theme(text = element_text(size=20))
    p = p + xlim(c(-M_matrix,M_matrix)) + ylim(c(-M_saige,M_saige))
    p = p + xlab("effect size (Matrix eQTL)") + ylab("effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # plot abs betas
    M_saige = max(abs(df_both$beta_saige))
    M_matrix = max(abs(df_both$beta_matrix))
    p = ggplot(df_both, aes(x=abs(beta_matrix), y=abs(beta_saige))) + geom_point(col=col, size=2) 
    p = p + geom_abline(slope = 1, intercept = 0, col = "gray") + theme_classic() 
    p = p + theme(text = element_text(size=20)) + xlim(c(0,M_matrix)) + ylim(c(0,M_saige))
    p = p + xlab("absolute effect size (Matrix eQTL)") + ylab("absolute effect size (SAIGE-QTL)")
    p = p + ggtitle(ct)
    print(p)
    # save
    pdf(paste0(fig_dir,ct,"_abs_betas_scatter.pdf"), width=8, height=8)
    print(p)
    dev.off()
    # print correlations
    cor_pvs = cor(-log10(df_both$pval_matrix), -log10(df_both$pval_saige))
    cor_betas = cor(df_both$beta_matrix, df_both$beta_saige)
    cor_abs_betas = cor(abs(df_both$beta_matrix), abs(df_both$beta_saige))
    print(paste0(ct,", cor pvs: ",round(cor_pvs,digits=2),", cor betas: ",round(cor_betas,digits=2)))
    print(paste0(ct,", cor abs betas: ",round(cor_betas,digits=2)))
    cors_pvs = c(cors_pvs, cor_pvs)
    cors_betas = c(cors_betas, cor_betas)
    cors_abs_betas = c(cors_abs_betas, cor_abs_betas)
}
