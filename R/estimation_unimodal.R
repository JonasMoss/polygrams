library("hisemi")
library("Rmosek")

# Iteratively calls polygram with different inputs.

polygram_unimodal = function(data, s = length(data)^(1/3), m = NULL, p = NULL, d = d, M = M, support = NULL,
                    symmetric = FALSE, monotone = NULL, shape = NULL, lower_boundary = NULL,
                    upper_boundary = NULL, strict = NULL) {

  len = length(s) + 1
  min_loss = Inf
  index = NA
  polygram_obj = NULL
  for (k in 1:(len)) {
    for (nu in 0:(m)) {
      a = polygram(data = data, s = s, m = m, p = p, d = d, M = M, symmetric = symmetric, monotone = monotone,
               shape = shape, unimodal = TRUE, lower_boundary = lower_boundary, support = support,
               upper_boundary = upper_boundary, nu = nu, k = k, strict = strict)
      loss = attr(a,"loss")
      if(loss < min_loss) {
        min_loss = loss
        index = c(k = k,nu = nu)
        polygram_obj = a
      }
    }
  }

  polygram_obj
}

polygram_biimodal = function(data, s = length(data)^(1/3), m = NULL, p = NULL,
                             symmetric = FALSE, monotone = NULL, shape = NULL) {

  len = length(s) + 1
  min_loss = Inf
  index = NA
  polygram_obj = NULL
  for (k in 1:(len)) {
    for (nu in 0:(m)) {
      a = polygram(data = data, s = s, m = m, p = p, symmetric = symmetric, monotone = monotone,
                   shape = shape, unimodal = TRUE, nu = nu, k = k)
      loss = attr(a,"loss")
      if(loss < min_loss) {
        min_loss = loss
        index = c(k = k,nu = nu)
        polygram_obj = a
      }
    }
  }

  polygram_obj
}

