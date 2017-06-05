#' Density, distribution function, quantile function, random generation, and
#' pth derivative for a Bernstein density with parameter vector \code{lambda}.
#'
#' @name distributions_bernstein_densities
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param lambda A vector of positive weights summing to 1.
#' @param support The support of the
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#' @param lower.tail logical; if  \code{TRUE} (default), probabilities are \eqn{P[x\leX]} otherwise, \eqn{P[X<x]}.
#' @details A Bernstein density on [a, b] is a subset Such polygrams can be used for non-parametric density estimation, as an alternative to logsplines
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
dxbernstein = function(x, lambda, p = 1, support = c(0, 1), log = FALSE) {
  dxbernstein_cpp(x, lambda, support, p, log)
}

#' @describeIn distributions_bernstein_densities Calculates a Bernstein
#' density.
dbernstein = function(x, lambda, p = 1, support = c(0, 1), log = FALSE) {
  dbernstein_cpp(x, lambda, support, log)
}

#' @describeIn distributions_bernstein_densities Calculates the distribution
#' function of a Bernstein density.
pbernstein = function(x, lambda, p = 1, support = c(0, 1), log.p = FALSE) {
  pbernstein_cpp(x, lambda, support, log.p)
}




