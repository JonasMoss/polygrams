#' The Polygram Class of Distributions
#'
#' Density, distribution function, random generation, and pth derivative of the
#' density for a polygram distribution object.
#'
#' @name Polygram
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If code{length(n) > 1}, the length is taken
#'    to be the number required.
#' @param object A \code{\link{polygram}} object.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as
#'    \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'    \eqn{P[x\leX]} otherwise, \eqn{P[X<x]}.
#' @details A Bernstein polygram on [a, b] is disjoint mixture of rescaled
#'   Bernstein densities on subintervals of [a, b]. Such polygrams can be used
#'   for non-parametric density estimation as an alternative to logsplines and
#'   kernel density estimation.
#' @return \code{dbernstein} gives the density, \code{pbernstein} the
#'   distribution function, \code{qbernstein} the quantile function,
#'   \code{rbeta} generates random deviates, and \code{dxbernstein} gives
#'   the pth derivative of the density function.
NULL

#' @rdname Polygram
#' @export
dxpolygram =function(x, object, p = 1, log = FALSE) {

  w = unlist(object,  use.names = FALSE)
  m = attr(object, "m")
  support = attr(object, "support")
  s = attr(object, "s")
  s_aug = c(support[1], s, support[2])
  dxpolygram_cpp(x = x, w = w, ms = m, s_aug = s_aug, p = p, log = log)

}

#' @rdname Polygram
#' @export
dpolygram = function(x, object, log = FALSE)
  dxpolygram(x, object = object, p = 0, log = log)


#' @rdname Polygram
#' @export
ppolygram = function(q, object, lower.tail = TRUE, log.p = FALSE) {

  out = dxpolygram(q, object = object, p = -1, log = log.p)

  if(!lower.tail) 1 - out else out

}

#' @rdname Polygram
#' @export
rpolygram = function(n, object) {

  support = attr(object, "support")
  s_aug   = c(support[1], attr(object, "s"), support[2])
  K       = length(s_aug) - 1
  lower_s = c(support[1], attr(object, "s"))
  weights = 1/sapply(1:(length(s_aug)-1), function(i) s_aug[i+1] - s_aug[i])
  m       = attr(object, "m")

  # The sampling procedure works by calling rbeta many times, every time with
  # argument (nu,m-nu) for each nu.

  N = sum(m + 1)

  probabilities = unlist(lapply(object, function(obj) {
    obj[obj < 0] = 0
    obj
  }))

  samples = tabulate(sample.int(n = N, size = n,
                                prob = probabilities, replace = TRUE))
  samples = c(samples, rep(0, N - length(samples)))
  samples = split(samples, rep(1:(length(m)), m + 1))
  names(samples) = NULL

  # The row sums of the samples matrix tells us how many times to call
  # rbeta(nu, m - nu).

  vector = unlist(sapply(1:K, function(index) {
    current_samples = unlist(sapply(0:m[index], function(nu) {
      rbeta(samples[[index]][nu+1], nu + 1, m[index] - nu + 1)
    }))
    current_samples*1/weights[index] + lower_s[index]
  }))

  return(vector)

}



