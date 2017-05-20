#' @export
#' @describeIn base_functions Calculates quantiles of a Bernstein polygram. NOT IMPLEMENTED.

qpolygram = function(p, polygram_object, lower.tail = TRUE, log.p = FALSE) {
  stop("qpolygram has not been correctly implemented yet.")

  # We will need to calculate the weights for each bin, which is done by using w.

  support = attr(polygram_object, "support")
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  m       = attr(polygram_object, "m")
  weights = sapply(polygram_object, sum)

  n = length(p)
  p = sort(p)
  return_vector = rep(NA,n)
  index = 1 # Keeps track of which bin we are in.
  i = 1     # Keeps track of which quantile we use.
  base = 0  # Sum of contributions from previous bins.
  polynomial_weights = polygram_object[[1]] # Vector of probabilities for the first bin.

  objective = function(y, m = NULL, polynomial_weights = NULL, rho = 0) {
    betas = sapply(0:m, function(nu) stats::pbeta(y, nu + 1, m + 1 - nu))
    (sum(betas*polynomial_weights) - rho)^2
  }

  while (i <= n) {
    while (TRUE) {
      if (p[i] > sum(weights[1:(index)])) {
        base = base + weights[index]
        index = index + 1
        print(p[i])
        print(index)
        polynomial_weights = polygram_object[[index]]
      } else {
        break
      }
    }

    extra = stats::nlm(objective, p = 0.5, m = m[index],
                polynomial_weights = polynomial_weights, rho = p[i] - base)$estimate
    return_vector[i] = extra
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
