#' Evaluates the k-th Bernstein basis density of degree m at x.
#'
#' @export
#' @param x A double where the polynomial is evaluated.
#' @param k Specifies which among the m+1 Bernstein polynomials to use.
#' Runs from 0,1 ... m.
#' @param m The degree of the Bernstein polynomial.
#' @param log = FALSE Returns the logarithm of the density if true.
#' @details Normalized Bernstein polynomials are positive and orthogonal on
#' [0,1], and integrates to 1.
#'
#' @return A double
dbernsteinbasis = function(x, k, m, log = FALSE) {
  if (!log) {
    choose(m, k)*x^k*(1-x)^(m-k)*(m+1)
  } else {
    lchoose(m,k) + k*log(x) + (m-k)*log(1-x) + log(m+1)
  }
}

#' Evaluates a m-ary Bernstein density with specified weights.
#'
#' @export
#' @param x A double where the density is evaluated.
#' @param lambda An m-ary vector of weights.
#' @param log Returns the logarithm of the density if true.
#' @details The vector lambda gives weights to the m Bernstein
#' polynomials of degree m.
#'
#' @return A double
dbernstein = function(x, lambda, log = FALSE) {
  m = length(lambda) - 1
  bernsteins = sapply(0:m, function(k) dbernsteinbasis(x, k = k, m = m))
  bernsteins = matrix(bernsteins, ncol = (m+1))
  density = bernsteins%*%lambda

  if(log) return(log(density))
  else return(density)
}

#' @export
#' @describeIn base_functions Calculates the density of a Bernstein polygram.
dpolygram = function(x, polygram_object, log = FALSE) {

  support = attr(polygram_object, "support")
  s       = c(attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])

  # If s_aug or weights are NULL, the polygram equals a
  # Bernstein density.

  if (length(s) == 1 & !is.null(s) & any(s == 0)) {
    if ( (ceiling(s) == s)) {
      s = as.integer(s)
      s = (1:s)/(s+1)
    }
  } else if (is.null(s) | all(s == 0)) {
    return(dbernstein(x, polygram_object[[1]], log = log))
  } else if (!("numeric" %in% class(s))){
    stop("s must either be numeric, null or an integer.")
  }

  x = sort(x)
  weights = sapply(polygram_object, sum)
  n = length(x)
  return_vector = rep(NA,n)
  index = 1
  i = 1

  while (i <= n) {
    while (TRUE) {
      if (x[i] > s_aug[index+1]) {
        index = index + 1
      } else {
        break
      }
    }

    y = (x[i] - s_aug[index])/(s_aug[index+1] - s_aug[index])
    current_weight = weights[index]/(s_aug[index+1] - s_aug[index])
    return_vector[i] = current_weight*dbernstein(y, polygram_object[[index]]/weights[index])
    i = i + 1

  }

  return(return_vector)

  }
