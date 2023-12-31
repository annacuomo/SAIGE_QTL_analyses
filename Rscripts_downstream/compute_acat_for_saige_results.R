# run as: Rscript compute_acat_for_saige_results.R {celltype}
# celltypes = c("B_IN","B_Mem","CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
#               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma")


# Code adpated from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value 
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

# new SAIGE-QTL results, cis results
mydir = "/share/ScratchGeneral/anncuo/OneK1K/saige_eqtl/from_wei/"

args <- commandArgs(trailingOnly=TRUE)
celltype <- args[1]
print(celltype)

ct_dir = paste0(mydir, "cis_", celltype, "/")
ct_files = list.files(ct_dir, pattern = "singleVar.txt")

gene_level_df = data.frame(gene = gsub("_.*", "", ct_files))
for (file in ct_files){
    filename = paste0(ct_dir,file)
    df_curr = as.data.frame(data.table:::fread(filename, header = T, stringsAsFactors = FALSE))
    df_curr = df_curr[df_curr$AF_Allele>0.05 & df_curr$AF_Allele<0.95,]
    gene = gsub("_.*","",file)
    gene_level_df[gene_level_df$gene == gene, "p.value.cct"] = get_CCT_pvalue(df_curr$p.value)
}

file = paste0(ct_dir, celltype, "_gene_acat_maf5_summary.csv")
data.table:::fwrite(gene_level_df, file)
