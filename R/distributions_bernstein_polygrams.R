#' Density, distribution function, quantile function, random generation, and
#' pth derivative of the density for a polygram object.
#'
#' @name base_functions
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If code{length(n) > 1}, the length is taken to be the number required.
#' @param polygram_object A polygram object.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if  \code{TRUE} (default), probabilities are \eqn{P[x\leX]} otherwise, \eqn{P[X<x]}.
#' @details A Bernstein polygram on [a, b] is disjoint mixture of rescaled Bernstein densities on subintervals
#' of [a, b]. Such polygrams can be used for non-parametric density estimation, as an alternative to logsplines
#' and kernel density estimation.
#' @return \code{dbernstein} gives the density, \code{pbernstein} the distribution function, \code{qbernstein}
#' the quantile function, \code{rbeta} generates random deviates, and \code{dxbernstein} gives the pth
#' derivative of the density function.
#' @examples
#' ## Do a parametric bootstrap.
#' set.seed(1337)
#' data = rbeta(200, 2, 7)
#' polygram_object = polygram(data, s = 3, m = 4)
#' current_median = qpolygram(0.5, polygram_object)
#' medians = replicate(100, {
#'   new_data = rpolygram(200, polygram_object)
#'   new_median = qpolygram(0.5, polygram(new_data, s = 3, m = 4))
#'   new_median
#' })
#' plot(polygram(sqrt(200)*(current_median - medians), s = 4, m = 4,
#'               support = c(-1,1)))
NULL

#' @describeIn distributions_bernstein_densities Calculates pth derivative of a
#'  Bernstein density.
dxpolygram = function(x, polygram_object, p = 1, support = c(0, 1), log = FALSE) {
  # This R function calls a C++ function which does not handle list objects
  # easily, hence the list in polygram_object will be changed into a vector,
  # and its split points stored in ms.

  w       = unlist(polygram_object,
                   use.names = FALSE)

  ms      = attr(polygram_object, "m")

  s_aug   = c(attr(polygram_object, "support")[1],
              attr(polygram_object, "s"),
              attr(polygram_object, "support")[2])

  deriv   = dxpolygram_cpp(x = x,
                           w = w,
                           ms = ms,
                           s_aug = s_aug,
                           p = p,
                           log = log)

  # Rcpp::Armadillo column vectors are interpreted as column vectors in
  # R as well, so they have dim(n, 1). The desired return type is vector,
  # not matrix, hence we return:

  as.vector(deriv)
}

#' @describeIn distributions_bernstein_densities Calculates a Bernstein
#' polygram density.
dpolygram = function(x, polygram_object, support = c(0, 1), log = FALSE) {
  dxpolygram(x, polygram_object = polygram_object, support = support, p = 0, log = log)
}

#' @describeIn distributions_bernstein_densities Calculates a Bernstein
#' polygram distribution function.
ppolygram = function(q, polygram_object, support = c(0, 1), lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail){
    dxpolygram(q, polygram_object = polygram_object, support = support, p = -1, log = log.p)
  } else {
    1 - dxpolygram(q, polygram_object = polygram_object, support = support, p = -1, log = log.p)
  }
}




