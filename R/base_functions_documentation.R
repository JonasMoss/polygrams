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





