#' @export
#' @describeIn base_functions Calculates the distribution function of a Bernstein polygram.

ppolygram = function(q, polygram_object, lower.tail = TRUE, log.p = FALSE) {

  # We will need to calculate the weights for each bin, which is done by using w.

  support = attr(polygram_object, "support")
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  m       = attr(polygram_object, "m")
  weights = sapply(polygram_object, sum)

  n = length(q)
  q = sort(q)
  return_vector = rep(NA,n)
  index = 1 # Keeps track of which bin we are in.
  i = 1     # Keeps track of which quantile we use.
  base = 0  # Sum of contributions from previous bins.
  polynomial_weights = polygram_object[[1]] # Vector of probabilities for the first bin.

  while (i <= n) {
    while (TRUE) {
      if (q[i] > s_aug[index+1]) {
        base = base + weights[index]
        index = index + 1
        polynomial_weights = polygram_object[[index]]
      } else {
        break
      }
    }

    y = (q[i] - s_aug[index])/(s_aug[index+1] - s_aug[index])
    betas = sapply(0:m[index], function(nu) stats::pbeta(y, nu + 1, m[index] + 1 - nu))
    return_vector[i] = base + sum(betas*polynomial_weights)
    i = i + 1

  }

  if(!lower.tail) {
    return_vector = 1 - return_vector
  }

  if(log.p) {
    return_vector = log(return_vector)
  }

  return(return_vector)
}
