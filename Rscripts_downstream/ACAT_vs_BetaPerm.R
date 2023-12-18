# this script compares performance of using the ACAT p-value 
# combination test as proposed in APEX vs the Beta approximation
# based on permutation proposed in the Fast QTL paper
# to obtain gene-level p-values in an eQTL study

# load R libraries
library(cowplot)
library(data.table)
library(ggplot2)

# ACAT implementation functions (from Wei)
# Code adapted from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value 
#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed.
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples CCT(pvals=pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association 115}(529), 393-402.
#' (\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
#' @export

CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  #if(is.zero && is.one){
  #  stop("Cannot have both 0 and 1 p-values!")
  #}
  if(is.zero){
    return(0)
  }
  if(is.one){
    #warning("There are p-values that are exactly 1!")
    return(min(1,(min(pvals))*(length(pvals))))
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

get_CCT_pvalue = function(pvalue, weights=NULL){
   pvals = pvalue
   notna = which(!is.na(pvals))
   if(length(notna) > 0){
     pvals = pvals[!is.na(pvals)]
     cctpval = CCT(pvals, weights=weights)
   }else{
     cctpval = NA
   }
   return(cctpval)
}

# Angli's TensorQTL results (implements a linear regression)
# using a slightly updated version of the OneK1K data
# * removing two low QC samples (tot n: 980)
# * removing from analysis samples (donor x cell type) with <5cells

# files copied from the following directory (in Brenner)
# Anglis_dir = "/share/ScratchGeneral/angxue/proj/vQTL/TensorQTL/scale_gene_wise/CD4_NC_cell_5_mean_mx/"
Anglis_dir = "/share/ScratchGeneral/anncuo/OneK1K/anglis_files/tensorqtl/"

# chromosome: 1
# cell type: CD4 NC 

# load file with significant results only
# top SNP per gene, contains beta-adjusted p-values
myfile = paste0(Anglis_dir,"OneK1K_CD4_NC.sig_cis_qtl_pairs.chr1.csv")
df_sig = read.csv(myfile, sep="\t")
nrow(df_sig)

# load file with all results 
# all SNPs, nominal p-values only
myfile = paste0(Anglis_dir,"OneK1K_CD4_NC.cis_qtl_pairs.chr1.csv")
df_all = read.csv(myfile, sep="\t", row.names=1)
nrow(df_all)

# create file for gene-specific p-values
# extract all unique genes
genes = unique(df_all$phenotype_id)
# make data frame
df_genes = data.frame(phenotype_id = genes)
nrow(df_genes)

# update data frame with CCT p-values (from nominal)
for (gene in genes){
    df_curr = df_all[df_all$phenotype_id == gene,]
    df_genes[df_genes$phenotype_id == gene,"pval_cct"] = get_CCT_pvalue(df_curr$pval_nominal)
}
head(df_genes)

# using data.table, combine with Beta approximated results
df_combined = as.data.frame(data.table(df_sig)[data.table(df_genes), on = c("phenotype_id"), nomatch=0])

# get p-value correlations
cor(-log10(df2$pval_beta+(10^(-16))), -log10(df2$pval_cct+(10^(-16))))

# plot p-values

p = ggplot(df2, aes(x=-log10(pval_beta), y=-log10(pval_cct))) + geom_point() 
p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
p = p + xlim(c(0,300)) + ylim(c(0,300)) 
p = p + theme(text = element_text(size=20))
p

# QQ plots

options(repr.plot.width = 12, repr.plot.height = 6) 
df2$pval_unif <- runif(dim(df2)[1], min = 0, max = 1)
p = ggplot(df2, aes(x=sort(-log10(pval_unif)), y=sort(-log10(pval_cct)))) + geom_point() 
p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
p = p + xlim(c(0,3)) + ylim(c(0,300)) 
p1 = p + theme(text = element_text(size=20))
p = ggplot(df2, aes(x=sort(-log10(pval_unif)), y=sort(-log10(pval_beta)))) + geom_point() 
p = p + geom_abline(slope = 1, intercept = 0, col = "firebrick") + theme_classic()
p = p + xlim(c(0,3)) + ylim(c(0,300)) 
p2 = p + theme(text = element_text(size=20))
plot_grid(p1, p2, ncol = 2)


