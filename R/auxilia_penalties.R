#' Calculates the penalty matrix for s, m and p.
#'
#' @param s A vector of splits. Sould not include end points.
#' @param m The order of the Bernstein polynomial.
#' @param p Derivative order.
#' @return A matrix

polygram_penalty_matrix = function(s, m, p) {

  K = length(s) + 1

  # The penalty matrix is a direct sum of K square matrices
  # of size m[k] x m[k]. Each of these matrices multplied with
  # a scalar after production.

  if(length(m) == 1) {
    m = rep(m, K)
  }

  s_aug = c(0, s, 1)

  Sigma = Vectorize(function(nu, eta, mk, p) {
    lower_bound_nu  = max(0, nu + p - mk)
    lower_bound_eta = max(0, eta + p - mk)
    upper_bound_nu  = min(nu, p)
    upper_bound_eta = min(eta, p)

    if ((lower_bound_nu  > upper_bound_nu) |
        (lower_bound_eta > upper_bound_eta)) {
      return(0)
    }

    is = lower_bound_nu:upper_bound_nu
    js = lower_bound_eta:upper_bound_eta

    summands = function(i, j) {
      first  = (-1)^(j+i)*choose(mk - p, nu - i)*choose(mk - p, eta - j)
      second = choose(p, i)*choose(p, j)
      third  = choose(2*(mk-p), nu - i + eta -j)
      first*second/third
    }

    sum(outer(is, js, FUN = summands))
  })

  matrix_list = lapply(1:K, function(k) {
    mk = m[k]
    wk = s_aug[k+1] - s_aug[k]
    mk_weight_a = (choose(mk, p)*factorial(p))^2*(mk - p + 1)^2
    mk_weight_b = (2*(mk - p) + 1)*wk^(2*p+1)
    mk_matrix = outer(0:mk, 0:mk, Sigma, mk = mk, p = p)
    mk_matrix*mk_weight_a/mk_weight_b
  })


  A = do.call(directSum, matrix_list)

  return(A)
}

