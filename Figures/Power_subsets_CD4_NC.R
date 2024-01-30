# Code below needs adjusting, but guidelines for power plot

pcts = c(1,5,10,20,50,100)
for (chrom in c(1:15)){
    df_summary = data.frame(pct_cells = pcts)
    for (pct in pcts){
        file = paste0(out_dir,'chr',chrom,'_',pct,'pct.csv')
        df = fread(file)
        df_summary[df_summary$pct_cells == pct,"n"] = nrow(df)
        df_summary[df_summary$pct_cells == pct,"n_saige_common"] = nrow(df[df$qv_saige_common<0.05,])
        df_summary[df_summary$pct_cells == pct,"n_tensor_common"] = nrow(df[df$qv_tensor_common<0.05,])
    }
    p = ggplot(df_summary, aes(x=pct_cells, y=n_saige_common)) 
    p = p + geom_line(col="forestgreen") + geom_point(size=3, col="forestgreen")
    p = p + geom_line(aes(x=pct_cells, y=n_tensor_common), col="cornflowerblue") 
    p = p + geom_point(aes(x=pct_cells, y=n_tensor_common), col="cornflowerblue", size=3) 
    p = p + theme_classic() + ggtitle(paste0("chrom",chrom)) + theme(text = element_text(size=20))
    p = p + xlab("% of cell subsetted") + ylab("# genes at FDR<5%")
    ypos = max(df_summary$n_saige_common)*0.3 
    p = p + annotate("text", label = "SAIGE-QTL", x=80, y=ypos, col="forestgreen", size=7)
    p = p + annotate("text", label = "TensorQTL", x=80, y=ypos*0.85, col="cornflowerblue", size=7)
    print(p)
}

fig_dir <- "/.../figures/"
pdf(paste0(fig_dir,"myplot.pdf"), width=8, height=6)
ggplot()
dev.off()
