#' Evaluates the p-th derivative of the k-th Bernstein basis density of
#' degree m at x.
#'
#' @export
#' @param x A double where the polynomial is evaluated.
#' @param k Specifies which among the m+1 Bernstein polynomials to use.
#' Runs from 0,1 ... m.
#' @param m The degree of the Bernstein polynomial.
#' @param p Degree of derivative.
#' @param log = FALSE Returns the logarithm of the density if true.
#'
#' @return A double
dxbernsteinbasis = Vectorize(function(x, k, m, p, log = FALSE) {
  y = NA
  if(p >= m) {
    y = x*0
  } else {
    lower_bound = max(0, k + p - m)
    upper_bound = min(k, p)
    if(lower_bound > upper_bound) {
      y = rep(0, length())
    } else {
      is = lower_bound:upper_bound
      pochhammer = choose(m, p)*factorial(p)
      ys = choose(p,is)*(-1)^(is+p)*dbernsteinbasis(x, k-is,m-p)
      y = pochhammer*sum(ys)
    }
  }

  if(!log) {
    return(y)
  } else {
    return(log(y))
  }
})

#' Evaluates the pth derivative of an m-ary Bernstein density with specified weights.
#'
#' @export
#' @param x A double where the density is evaluated.
#' @param lambda An m-ary vector of weights.
#' @param p Desired derivative. Defaults to 1.
#' @param log Returns the logarithm of the density if true.
#' @details The vector lambda gives weights to the m Bernstein
#' polynomials of degree m.
#'
#' @return A double
dxbernstein = function(x, lambda, p = 1, log = FALSE) {
  m = length(lambda) - 1
  bernsteins = sapply(0:m, function(k) dxbernsteinbasis(x, k = k, m = m, p = p))
  bernsteins = matrix(bernsteins, ncol = (m+1))
  density = bernsteins%*%lambda

  if(log) return(log(density))
  else return(density)
}


#' @export
#' @describeIn base_functions Calculates the pth derivative of a Bernstein polygram.
dxpolygram = function(x, polygram_object, p, log = FALSE) {

  if(!is.integer(p) & !is.numeric(p)) {
    stop("Option 'p' must be an integer greater than -2.")
  } else if (p < -1) {
    stop("Option 'p' must be an integer greater than -2.")
  }

  if(p == 0) {
    return(dpolygram(x, polygram_object, log = log))
  } else if (p == -1) {
    return(ppolygram(x, polygram_object, log.p = log))
  }

  support = attr(polygram_object, "support")
  m       = c(attr(polygram_object, "m"))
  s       = c(attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])

  if(all(p >= m)) {
    return(rep(0,length(x)))
  }

  # If s_aug or weights are NULL, the polygram equals a
  # Bernstein density.

  if (length(s) == 1 & !is.null(s) & any(s == 0)) {
    if ( (ceiling(s) == s)) {
      s = as.integer(s)
      s = (1:s)/(s+1)
    }
  } else if (is.null(s) | all(s == 0)) {
    return(dxbernstein(x, polygram_object[[1]], p = p, log = log))
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

    # This loop iterates through x and updates the index counter whenever a new
    # bin is encountered.

    while (TRUE) {
      if (x[i] > s_aug[index+1]) {
        index = index + 1
      } else {
        break
      }
    }

    # This is the actual derivatives.
    y = (x[i] - s_aug[index])/(s_aug[index+1] - s_aug[index])
    current_weight = weights[index]/(s_aug[index+1] - s_aug[index])^(p+1)
    return_vector[i] = current_weight*dxbernstein(y, lambda = polygram_object[[index]]/weights[index], p = p)
    i = i + 1

  }

  return(return_vector)

}
