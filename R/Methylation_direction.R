#' Show the direction of methylation about which CpG site is hypermethylated or hypomethylated
#'
#' This function requires the beta matrix, list of differential CpG sites, and boolean vector for the contrast.
#' The boolean vector must have the same length as the number of columns of the beta matrix.
#'
#' @param beta_matrix Matrix of beta values.
#' @param gbdmr_list gbdmr's output or vector of differential CpGs.
#' @param contrast_boolean Boolean vector showing which columns to contrast.
#'
#' @return Character vector containing hypermethylated or hypomethylated.
#' @export

Methylation_direction = function(beta_matrix, gbdmr_list, contrast_boolean) {
  if (NCOL(beta_matrix) != length(contrast_boolean)) {
    stop("Beta matrix's column and the contrast input must have the same size.")
  }
  
  if (class(gbdmr_list) == "character") {
    differential.cpg = intersect(rownames(beta_matrix), gbdmr_list)
  }
  
  if (class(gbdmr_list) == "numeric" | class(gbdmr_list) == "boolean") {
    differential.cpg = rownames(beta_matrix)[gbdmr_list]
  }
  
  if (class(gbdmr_list) == "list") {
    if (names(gbdmr_list)[5] == "positive_cpg") {
      differential.cpg = gbdmr_list$positive_cpg
    }
  }
  
  averagebeta.T = rowMeans(subset(beta[differential.cpg, ], select=contrast_boolean))
  averagebeta.F = rowMeans(subset(beta[differential.cpg, ], select=!contrast_boolean))
  
  deltabeta = averagebeta.T - averagebeta.F
  
  return(ifelse(deltabeta >= 0, "Hypermethylated", "Hypomethylated"))
}
