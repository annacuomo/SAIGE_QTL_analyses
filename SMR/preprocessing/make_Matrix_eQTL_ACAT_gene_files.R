library(data.table)
library(qvalue)

#Code adpated from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value 
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

# Matrix eQTL results (from original publication, Yazar et al Science 2022)
matrix_dir = "/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/smr_analysis/matrix_eQTL/"

# 14 cell types
celltypes = c('BimmNaive','Bmem','CD4all','CD4effCM','CD4TGFbStim','CD8all','CD8eff',
              'CD8unknown','DC','MonoC','MonoNC','NKact','NKmat','Plasma')

# given a celltype, chromosome, open file
# for each gene, calculate ACAT combined gene-level p-value from SNP-level p-values
get_acat_results <- function(celltype, chrom){
    ct_dir = paste0(matrix_dir, celltype, "/output_files/")
    matrix_file = paste0(ct_dir,celltype,"_chr",chrom,"_cis_eqtls_210211.tsv")
    matrix_df = read.csv(matrix_file, sep="\t")
    genes = unique(matrix_df$gene)
    df_gene = data.frame(gene = genes)
    for (gene in genes){
        df_curr = matrix_df[matrix_df$gene == gene,]
        df_gene[df_gene$gene == gene,"pval_cct"] = get_CCT_pvalue(df_curr$p.value)
    }
    return(df_gene)
}

# loop over celltypes and combine all results for each ct
# across chromosomes, then save
out_dir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/genes/matrixeqtl_egenes/"
for (celltype in celltypes){
    df_list = list()
    for (chrom in 1:22){
        df = get_acat_results(celltype = celltype, chrom = chrom)
        df$chrom = chrom
        df_list[[chrom]] = df 
    }
    df_all = rbindlist(df_list)
    df_all$qv = qvalue(df_all$pval_cct)$qvalues
    out_file = paste0(out_dir, celltype,".csv")
    fwrite(df_all, out_file)
}
