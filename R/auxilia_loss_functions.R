#' Calculates the empirical or theoretical Bernstein contributions for the
#' data "data" (or density "d"), splits vector s, with respect to
#' Bernstein polynomials of order m. Used as a helper function in the
#' polygram_objective_vector function.
#'
#' @param data The data. Can be either a
#' @param s The vector of splits.
#' @param m The order of the Bernstein polynomial.
#' @param d Optional density, which will override the data if supplied.
#' @return An (|s|-1) x (m+1)-dimensional matrix
empirical_bernstein = function(data, s, m, d = NULL) {

  # If a d is null, the data is used. If not, the density is used.

  if (is.null(d)) {
    # I begin by splitting the data set into |s|-1 groups
    # by using the splits vector s.

    data = sort(data)
    N = length(data)

    if (is.null(s)) {
      FUN = Vectorize(function(nu, index) {
        1/N*sum(dbernsteinbasis(data, nu, m))
      })

      matrix_ = outer(0:m, 1, FUN = FUN)
    } else {

      s_aug = c(0, s, 1)

      # Now we will group the data!

      len = length(s) + 1
      bin_index = 1
      groupings = replicate(len, rep(1,0))


      for (i in 1:length(data)) {
        while (TRUE) {
          if (data[i] > s_aug[bin_index+1]) {
            bin_index = bin_index + 1
          } else {
            break
          }
        }
        groupings[[bin_index]] = c(groupings[[bin_index]],
                                   (data[i] -s_aug[bin_index])/(s_aug[bin_index+1]-s_aug[bin_index]))
      }

      FUN = Vectorize(function(nu, index) {
        1/N*sum(dbernsteinbasis(groupings[[index]], nu, m))
      } )

      matrix_ = outer(0:m, 1:len, FUN = FUN)
    }
  } else {
    # When no data is supplied, we integrate the bernstein polynomials w.r.t
    # the derivative d.
    len = length(s) + 1
    s_aug = c(0,s,1)
    integral = Vectorize(function(nu, i) {
      integrate(function(x) dbernsteinbasis((x-s_aug[i])/(s_aug[i+1]-s_aug[i]), nu, m)*d(x),
                lower = s_aug[i], upper = s_aug[i+1])$value
    })

    matrix_ = outer(0:m, 1:len, FUN = integral)

  }

  s_aug = c(0,s,1)
  weights = 1/sapply(1:(length(s_aug) - 1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  t(t(matrix_) * weights)
}

#' Returns the true L2-loss for a weight matrix V and underlying distribution
#' d. Useful for theoretical investigations or simulation studies.
#'
#' @param d The true density.
#' @param s The vector of splits. If a natural number K, s is chosen as
#' (1:K)/(K+1). Defaults to length(data)^(1/3)
#' @param V Weight matrix for the Bernstein polygram.
#' @param discrepancy If true, only returns the discrepancy.
#' @return A real number, the L2-loss.

polygram_L2_loss = function(d, polygram_object, discrepancy = TRUE) {

  s       = c(attr(polygram_object, "s"))
  m       = attr(polygram_object, "m")
  v       = unlist(polygram_object)

  first = integrate(function(x) d(x)^2, lower = 0, upper = 1)$value
  middle = -sum(2*empirical_bernstein(NULL, s, m, d)*polygram_object)
  last = t(v)%*%polygram_objective_matrix(m, s)%*%v

  if (!discrepancy) {
    first + middle + last
  } else {
    middle + last
  }
}


#' Returns the penalty of order p for a polygram object.
#'
#' @param polygram_object A polygram object.
#' @param p Derivative. Only supported for p = 0 right now.
#' @return The pth penalty.

penalty = function(polygram_object, p = 0) {
  obj_      = polygram_object
  s         = c(attr(polygram_object, "s"))
  m         = attr(polygram_object, "m")
  obj_      = unlist(obj_)
  mat       = polygram_objective_matrix(m, s)
  c(t(obj_)%*%mat%*%obj_)
}

intersections = function(s1, s2) {
  # This is a list of |s1| vectors which specify which elements of |s2| intersects
  # which elements of |s1|
  intersections = c()
  s1_index = 1
  s2_index = 1
  s1_upper = c(s1, 1)
  s2_aug   = c(0, s2, 1)
  s2_lower = c(0, s2)

  while(s1_index <= length(s1_upper) & s2_index <= length(s2_lower)) {
    # On a new iteration the current s1 bin always intersects the current
    # s2 bin.
    intersections = c(intersections, s1_index, s2_index)

    # If the upper limit of the current s2 bin is larger than the upper
    # limit of the current s1 bin, the s1 index should be increased
    if(s2_aug[s2_index+1] > s1_upper[s1_index]) {
      s1_index = s1_index + 1
    } else if (s2_lower[s2_index] <= s1_upper[s1_index]) {
      s2_index = s2_index + 1
    }
  }

  dim(intersections) = c(2, length(intersections)/2)
  intersections = t(intersections)
  colnames(intersections) = c("k","j")
  intersections
}

#' Calculates the L2-loss between two polygram objects. The first argument
#' is interpreted as the "true" polygram. Only supported when m are equal for the
#' moment.
#'
#' @param polygram_object_1 A polygram object, interpreted as the true density.
#' @param polygram_object_2 A polygram object.
#' @param discrepancy If true, only returns the discrepancy.
#' @return The pth penalty.
#'
polygram_vs_loss = function(polygram_object_1, polygram_object_2, discrepancy = FALSE) {

  s1       = c(attr(polygram_object_1, "s"))
  s2       = c(attr(polygram_object_2, "s"))
  s1_aug   = c(0, s1, 1)
  s2_aug   = c(0, s2, 1)
  m        = attr(polygram_object_1, "m")

  v1       = polygram_object_1
  dim(v1)  = length(polygram_object_1)
  v2       = polygram_object_2
  dim(v2)  = length(polygram_object_2)

  first = t(v1)%*%polygram_objective_matrix(s1, m)%*%v1
  last  = t(v2)%*%polygram_objective_matrix(s2, m)%*%v2


  integral_pair = function(k, j) {
    pairs = expand.grid(0:m, 0:m)

    mat = apply(pairs, 1, function(pair) {
      nu = pair[1]
      eta  = pair[2]
      lower = max(s1_aug[k], s2_aug[j])
      upper = min(s1_aug[k+1], s2_aug[j+1])
      integrate(function(x) dbeta((x-s1_aug[k])/(s1_aug[k+1]-s1_aug[k]), nu + 1, m + 1 - nu)*
                  dbeta((x-s2_aug[j])/(s2_aug[j+1]-s2_aug[j]), eta + 1, m + 1 - eta),
                lower = lower, upper = upper)$value*1/(s1_aug[k+1]-s1_aug[k])*1/(s2_aug[j+1]-s2_aug[j])
    })

    dim(mat) = c(m+1,m+1)
    mat
  }

  ints = intersections(s1, s2)
  mats = lapply(1:nrow(ints), function(i) integral_pair(k = ints[i, 1], j = ints[i, 2]))
  mat_prods = sapply(1:nrow(ints), function(i) polygram_object_1[,ints[i,1]]%*%mats[[i]]%*%polygram_object_2[,ints[i,2]])
  middle = -2*sum(mat_prods)


  if (!discrepancy) {
    first + middle + last
  } else {
    middle + last
  }
}


#' Obtains the empirical discrepancy of a polygram object. If d is non-null,
#' it is calculated with d as the true density.
#'
#' @export
#' @param polygram_object A polygram object.
#' @return The least squares discrepancy.
discrepancy = function(polygram_object, d = NULL) {
  if(is.null(d)) attr(polygram_object,"loss")
  else polygram_L2_loss(d = d, polygram_object, discrepancy = TRUE)
}

#' Makes an elbow plot for a polygram object.
#'
#' @export
#' @param data The underlying data
#' @param N The maximal number of bins.
#' @param type Can be "regular" or "quantiles".
#' @param loss The type of discrepancy used. Defaults to L2, and KL will be
#' supported.
#' @param ... Passed to polygram.
#' @return An elbow plot.

elbow_plot = function(data, lower = NULL, upper = NULL, type = "regular", loss = "L2", ...) {

  n = length(data)
  if(is.null(lower)) lower = 0
  if(is.null(upper)) upper = floor(2*sqrt(n))

  if(type == "quantiles") {
    losses = sapply(lower:upper, function(s) {
      if(s == 0) s = NULL
      else s = quantile(data, (1:s)/(s+1))
      obj = polygram(data, s = s, ...)
      discrepancy(obj)
    })
  } else {
    losses = sapply(lower:upper, function(s) {
      if(s == 0) s = NULL
      tryCatch({
        obj = polygram(data, s = s, ...)
        discrepancy(obj)
      }, error = function(e) NA)
    })
  }

  na_omit_losses = na.omit(losses)
  bounds = c(min(na_omit_losses), max(na_omit_losses))
  xx = lower:upper
  plot(xx, losses, ylim = bounds, main = "Elbow plot", xlab = "k",
       ylab = paste0("Empirical ", loss, "-discrepancy"), bty = "l")
  grid()
}
