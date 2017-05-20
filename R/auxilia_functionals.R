#' Moments, central moments, mean, variance, standard deviation, skewness and
#' (excess) kurtosis for Bernstein polygrams.
#'
#' @name moment_functionals
#' @param polygram_object a polygram object.
#' @param moment positive integer, specifying which moment to calculate.
#' @param excess logical, only for kurtosis; if  \code{TRUE} (default), the excess kurtosis is calculated.
#' @details A Bernstein polygram on [a, b] is disjoint mixture of rescaled Bernstein densities on subintervals
#' of [a, b]. Such polygrams can be used for non-parametric density estimation, as an alternative to logsplines
#' and kernel density estimation. These functions help with calculating moments.
#' @return \code{mean} returns the expected value, \code{var} the variance, \code{sd} the standard deviation,
#' \code{moment} the pth moment, \code{central_moment} the pth central moment, \code{skewness} the skewness,
#' and \code{kurtosis} the kurtosis.
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

#' @export
#' @describeIn moment_functionals Calculates the mean of a Bernstein polygram.
mean_polygram = function(polygram_object) {

  support = attr(polygram_object, "support")
  s_lower = c(support[1], attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  K       = length(s_aug) - 1
  m       = attr(polygram_object, "m")
  weights = sapply(1:(length(s_aug)-1), function(i) s_aug[i+1] - s_aug[i])

  nus   = 0:m
  means = (nus+1)/(m+2)
  first = sum(apply(polygram_object, 2, function(x) sum(x*means))*weights)
  last  = sum(colSums(polygram_object)*s_lower)
  first + last

}

#' @export
#' @describeIn moment_functionals Calculates the pth moment of a Bernstein polygram.
moment_polygram = function(polygram_object, moment = 1) {
  support = attr(polygram_object, "support")
  s_lower = c(support[1], attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  K       = length(s_aug) - 1
  m       = attr(polygram_object, "m")
  weights = sapply(1:(length(s_aug)-1), function(i) s_aug[i+1] - s_aug[i])

  S_matrix = outer(0:moment, 1:K, FUN =
                     function(i,j) choose(moment,i)*weights[j]^i*s_lower[j]^(moment - i))
  PI_matrix = outer(0:m, 1:moment, FUN =
                     Vectorize(function(nu, r) prod((nu + 1:r)/(m + 1 + 1:r))))
  PI_matrix = cbind(rep(1,m+1), PI_matrix)
  Sigma_matrix = (PI_matrix %*% S_matrix)

  sum(Sigma_matrix * polygram_object)
}

#' @export
#' @describeIn moment_functionals Calculates the pth central moment of a Bernstein polygram.
central_moment_polygram = function(polygram_object, moment = 2) {
  if (moment == 1) return(0)
  if (moment == 2) {
    variance = moment.polygram(polygram_object, 2) - moment.polygram(polygram_object, 1)^2
    return(variance)
  }
}

#' @export
#' @describeIn moment_functionals Calculates the variance of a Bernstein polygram.
var_polygram = function(polygram_object) {
  central_moment.polygram(polygram_object, moment = 2)
}

#' @export
#' @describeIn moment_functionals Calculates the standard deviation of a Bernstein polygram.
sd_polygram = function(polygram_object) {
  sqrt(var.polygram(polygram_object))
}

#' @export
#' @describeIn moment_functionals Calculates the skewness of a Bernstein polygram.
skewness_polygram = function(polygram_object) {
  sigma = sd.polygram(polygram_object)
  third_moment = moment.polygram(polygram_object, moment = 3)
  mu = mean.polygram(polygram_object)
  skewness = 1/sigma^3*(third_moment - 3*mu*sigma^2 - mu^3)

  return(skewness)
}

#' @export
#' @describeIn moment_functionals Calculates the (excess) kurtosis of a Bernstein polygram.
kurtosis_polygram = function(polygram_object, excess = FALSE) {
  second_moment = moment.polygram(polygram_object, moment = 2)
  third_moment = moment.polygram(polygram_object, moment = 3)
  fourth_moment = moment.polygram(polygram_object, moment = 4)
  mu = mean.polygram(polygram_object)
  var = second_moment - mu^2
  kurtosis = 1/var^2*(fourth_moment - 4*mu*third_moment +
                        6*mu^2*second_moment - 3*mu^4)

  if(excess) {
    return(kurtosis - 3)
  } else {
    return(kurtosis)
  }
}
