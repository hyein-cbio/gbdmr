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
  
  vector_input = is.null(dim(gbdmr_list))
  if (vector_input && is.character(gbdmr_list)) {
    differential.cpg = intersect(rownames(beta_matrix), gbdmr_list)
  } else if (vector_input && (is.numeric(gbdmr_list) || is.logical(gbdmr_list))) {
    differential.cpg = rownames(beta_matrix)[gbdmr_list]
  } else if (is.list(gbdmr_list) && !is.data.frame(gbdmr_list) &&
             identical(names(gbdmr_list)[5], "positive_cpg")) {
    differential.cpg = gbdmr_list$positive_cpg
  } else {
    stop("gbdmr_list must be a character vector of CpG ids, a numeric or logical index vector, or a gbdmr output list with \"positive_cpg\" as its fifth element.")
  }
  
  averagebeta.T = rowMeans(beta_matrix[differential.cpg, contrast_boolean, drop = FALSE])
  averagebeta.F = rowMeans(beta_matrix[differential.cpg, !contrast_boolean, drop = FALSE])
  
  deltabeta = averagebeta.T - averagebeta.F
  
  return(ifelse(deltabeta >= 0, "Hypermethylated", "Hypomethylated"))
}
