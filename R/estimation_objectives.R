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

polygram_objective_vector_ = function(data, s, m, d = NULL, M = NULL) {

  matrix_ = outer(0:m, 1:(length(s) + 1))

  if (!is.null(data)) {
    # I begin by splitting the data set into |s|-1 groups
    # by using the splits vector s.

    data = sort(data)
    N = length(data)

    if (is.null(s)) {
      FUN = Vectorize(function(nu, index) {
        1/N*sum(dbernsteinbasis(data, nu, m))
      })

      matrix_ = outer(0:m, 1, FUN = FUN)
    } else {

      s_aug = c(0, s, 1)

      # Now we will group the data!

      len = length(s) + 1
      bin_index = 1
      groupings = replicate(len, rep(1,0))


      for (i in 1:length(data)) {
        while (TRUE) {
          if (data[i] > s_aug[bin_index+1]) {
            bin_index = bin_index + 1
          } else {
            break
          }
        }
        groupings[[bin_index]] = c(groupings[[bin_index]],
                                   (data[i] -s_aug[bin_index])/(s_aug[bin_index+1]-s_aug[bin_index]))
      }

      FUN = Vectorize(function(nu, index) {
        1/N*sum(dbernsteinbasis(groupings[[index]], nu, m))
      } )

      matrix_ = outer(0:m, 1:len, FUN = FUN)
    }
  }

  if (!is.null(d)) {
    N = length(data)
    len = length(s) + 1
    s_aug = c(0,s,1)
    integral = Vectorize(function(nu, i) {
      integrate(function(x) dbernsteinbasis((x-s_aug[i])/(s_aug[i+1]-s_aug[i]), nu, m)*d(x),
                lower = s_aug[i], upper = s_aug[i+1])$value
    })

    matrix_ = N/(N+M)*matrix_ + M/(N+M)*outer(0:m, 1:len, FUN = integral)
  }

  s_aug = c(0,s,1)
  weights = 1/sapply(1:(length(s_aug) - 1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  vector = t(t(matrix_) * weights)

  dim(vector) = c((m+1)*(length(s)+1))

  return(vector)
}

#' Calculates the vector B used in the optimisation procedure. It supports both
#' data and distribution functions as input.
#'
#' @param data The data.
#' @param s The vector of splits.
#' @param m The order of the Bernstein polynomial.
#' @param d An optional density.
#' @return An (|s|+1)*(m+1) vector.
#'
polygram_objective_vector = function(data, ms, s) {

  if(length(s) != (length(ms)-1)) {
    stop("length(s) + 1 must equal length(ms)!")
  }

  data = sort(data)
  N = length(data)

  if (is.null(s)) {
    FUN = Vectorize(function(nu, index) {
      1/N*sum(dbernsteinbasis(data, nu, ms))
    })

    matrix_ = outer(0:ms, 1, FUN = FUN)
  } else {
    s_aug = c(0, s, 1)
    K = length(s) + 1

    weights = 1/sapply(1:(length(s_aug) - 1), function(i) {
      s_aug[i+1] - s_aug[i]
    })

    len = length(s) + 1
    bin_index = 1
    groupings = replicate(len, rep(1,0))

    for (i in 1:length(data)) {
      while (TRUE) {
        if (data[i] > s_aug[bin_index+1]) {
          bin_index = bin_index + 1
        } else {
          break
        }
      }
      groupings[[bin_index]] = c(groupings[[bin_index]],
              (data[i] -s_aug[bin_index])/(s_aug[bin_index+1]-s_aug[bin_index]))
    }

    vec_list = lapply(1:K, function(k) {
      vec = sapply(0:ms[k], function(nu) {
        1/N*sum(dbernsteinbasis(groupings[[k]], nu, ms[k]))
      })
      vec*weights[k]
    })
    vector = unlist(vec_list)
  }

  return(vector)

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

