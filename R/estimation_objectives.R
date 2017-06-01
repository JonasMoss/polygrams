#' Transforms a symmetric matrix into the form demanded by Rmosek.
#'
#' @param mat A matrix.
#' @return A matrix in triangular sparse triplet form.

as.mosek_mat = function(mat) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix needed for this function to work. Please install it.",
         call. = FALSE)
  }

  new = mat*0
  new[upper.tri(mat, diag = TRUE)] = mat[upper.tri(mat, diag = TRUE)]
  new = Matrix::as(new, "dtTMatrix")
  mat_ = list()
  mat_$i = attr(new, "j") + 1
  mat_$j = attr(new, "i") + 1
  mat_$v = attr(new, "x")
  mat_
}


#' Calculates the objective matrix for use in a polygram.
#'
#' @param ms Vector of Bernstein orders.
#' @param s The vector of splits.

polygram_objective_matrix = function(ms, s) {
  # The objective matrix is a weighted direct sum of |s| + 1
  # matrices on the same form: A_

  A_ = function(m) {

    FUN = Vectorize(function(eta,nu) {
      (m+1)/(choose(2*m+1,m+1))*choose(eta+nu,eta)*choose(2*m-eta-nu,m-eta)
    })

    outer(0:m, 0:m, FUN)
  }

  # The next step is to do the direct sum of the weighted matrices.
  s_aug = c(0,s,1)

  weights = 1/sapply(1:(length(s_aug) - 1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  K = length(s) + 1

  A = do.call(direct_sum, lapply(1:K, function(k) weights[k]*A_(ms[k])))

  return(A)
}

