#' Samples from a polygram.
#'
#' @param n Number of desired samples.
#' @param polygram_object A polygram object fitted by polygram.
#' @return A vector of n random samples.

rpolygram = function(n, polygram_object) {
  support = attr(polygram_object, "support")
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  K       = length(s_aug) - 1
  lower_s = c(support[1], attr(polygram_object, "s"))
  weights = 1/sapply(1:(length(s_aug)-1), function(i) s_aug[i+1] - s_aug[i])
  m       = attr(polygram_object, "m")

  # The sampling procedure works by calling rbeta many times, every time with
  # argument (nu,m-nu) for each nu.

  N = sum(m + 1)

  probabilities = unlist(lapply(polygram_object, function(obj) {
    obj[obj < 0] = 0
    obj}))

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
