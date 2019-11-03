#' Calculates the direct sum of matrices.
#' @param ... the matrices that should be summed.
#' @return The direct sum of the matrices in \code{...}.
#' @details The direct sum of matrices A, B, C, D ... is obtained by pasting
#'    them along the diagonal.
#' @keywords internal
#' @examples \dontrun{
#'    A = matrix(1:9, ncol = 3)
#'    B = matrix(1:20, ncol = 4)
#'    Sigma = direct_sum(A, B)
#' }

direct_sum = function(...) {

  matrices = list(...)

  are_matrices = lapply(matrices, function(mat) "matrix" %in% class(mat))
  if(!all(unlist(are_matrices))) stop("direct_sum can only take matrix inputs.")

  n_matrices = length(matrices)
  cum_ncol = c(0,cumsum(sapply(matrices, function(mat) ncol(mat))))
  cum_nrow = c(0,cumsum(sapply(matrices, function(mat) nrow(mat))))

  Sigma = matrix(0, ncol = cum_ncol[n_matrices+1],
                 nrow = cum_nrow[n_matrices+1])

  for (i in 1:length(matrices)) {
    row_indices = (cum_nrow[i]+1):cum_nrow[i+1]
    col_indices = (cum_ncol[i]+1):cum_ncol[i+1]
    Sigma[row_indices, col_indices] = matrices[[i]]
  }

  return(Sigma)

}
