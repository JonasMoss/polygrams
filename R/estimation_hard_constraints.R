## ============================================================================
## This file contains functions that enforce different sorts of contraints.
## There are three different kinds of functions, and each of them follow a
## uniform setup.
##
## 1.) Linear Matrix constraints.
##   These functions correspond to linear constraints that require matrices to
##   work, so they are not box constraints. Each function returns a list
##   containing three items. (Or NULL, corresponding to no constraints.)
##   - constraint: The constraint matrix.
##   - lower: The lower bound of constraint %*% w.
##   - upper: The upper bound of constraint %*% w.
##   It is important that these have the correct dimensions, and both lower and
##   must be supplied, even if only is needed. (The matrices can always be
##   constructed so that only one is needed, but is this not the approach taken
##   here.) The output of these functions can be sent into a penalizer together
##   with a lambda in order to make penalizers of them. (Sometimes referred to
##   as soft constraints.)
##
##   Currently this section holds:
##
##   _polygram_connectivity_constraint(s, m, p) :
##      Here p can be either a vector or a natural number. Returns constraints
##      so that the ith split is connected to order p[i]. The most important
##      function obviously!
##
##   _polygram_monotone_constraints(s, m, p, monotone):
##      Returns monotonicity constratins decided by "monotone", which can be
##      either increasing or decreasing. An additional parameter specifying
##      over what bins it should be monotonic in conjunction with a vector
##      for the monotone parameter might be implemented later, in addition
##      to a strictness parameter.
##
##   _polygram_shape_constraints(s, m, p, monotone):
##      Shape is either convex or concave.
##
##   _
## 2.) Box constraint. Involves modification of bounds for each separate
##   variable. Only used for forcing boundary constraints. (But these are
##   quite important!)
##
## 3.) Quadratic constraints. None of these are implemented, but can be used to
##   control the variance, for instance.
##
## ============================================================================

#' Returns the connectivity constraints for Cp.
#'
#' @param s The vector of splits.
#' @param m The order of the Bernstein polynomial.
#' @param p Desired order of connectivity. Must be strictly less than m.
#' @return A matrix of constraints.
#'

get_connectivity_constraints = function(ms, s, ps, directions = NULL) {
  # If m is equal to zero at all places, the return value is NULL.
  if(all(ms == 0) | all(ps == -1) | is.null(s) | length(s) == 0) {
    return(NULL)
  }

  K = length(s) + 1      # K is the number of bins in the polygram.
  s_aug = c(0, s ,1)     # The support is always assumed to be (0, 1).

  weights = sapply(1:(K+1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  lower_matrices = lapply(1:(K-1), function(split_index) {
    if(ps[split_index] == -1) {
      lower = rep(0,ms[split_index]+1)
      dim(lower) = c(1, length(lower))
      lower
    } else {t(sapply(0:ps[split_index], function(rho) {
      m = ms[split_index]

      # This fills up the matrix with zero rows when the
      # desired derivative has to be zero and the current density
      # is of lower order than the derivative itself. In addition, it
      # must be zero when rho is larger than the next density's m!

      if(rho >= (m+1)) {
        return(rep(0, m+1))
      }

      nus = seq(0,rho)
      signs = (-1)^(nus + rho)
      binoms = choose(rho, nus)

      lower = c(rep(0,ms[split_index]-rho),
                weights[split_index+1]^(rho+1)*rev(binoms)*signs)
      lower*prod((m-rho+1):(m+1))
    }))
   }
  })

  upper_matrices = lapply(2:K, function(split_index) {
    if(ps[split_index-1] == -1) {
      upper = rep(0,ms[split_index]+1)
      dim(upper) = c(1, length(upper))
      upper
    } else {
      t(sapply(0:ps[split_index-1], function(rho) {
      m = ms[split_index]

      if(rho >= (m+1)) {
        return(rep(0, m+1))
      }

      nus = seq(0,rho)
      signs = (-1)^(nus + rho)
      binoms = choose(rho, nus)

      upper = -c(weights[split_index-1]^(rho+1)*binoms*signs,
                 rep(0,ms[split_index]-rho))
      upper*prod((m-rho+1):(m+1))
      }))
    }
  })

  lower_matrix_ = do.call(hisemi::directSum, lower_matrices)
  lower_matrix  = cbind(lower_matrix_,
                       matrix(0, ncol = ms[K] + 1, nrow = nrow(lower_matrix_)))

  upper_matrix_  = do.call(hisemi::directSum, upper_matrices)
  upper_matrix  = cbind(matrix(0, ncol = ms[1] + 1, nrow = nrow(upper_matrix_)),
                       upper_matrix_)

  constraint_matrix = lower_matrix + upper_matrix

  lower = rep(0, nrow(constraint_matrix))
  upper = rep(0, nrow(constraint_matrix))

  if(is.null(directions)) {
    directions = rep("both", length(ps))
  }

  if(length(directions) == 1) {
    directions = rep(directions, length(ps))
  }

  current_index = 1
  for(split_index in 1:(K-1)) {
    current_length = nrow(lower_matrices[[split_index]])

    if(directions[split_index] == "left") {
      upper[current_index:(current_index + current_length - 1)] = Inf
    } else if(directions[split_index] == "right") {
      lower[current_index:(current_index + current_length - 1)] = -Inf
    }

    current_index = current_index + current_length
  }

  list(constraint = constraint_matrix,
       lower = lower,
       upper = upper)
}

#' Returns the edge constraints.
#'
#' @param s The vector of splits.
#' @param ms The order of the Bernstein polynomial.
#' @param ps A list of desired values for the derivatives.
#' @param directions Directions for the derivatives. Can be "equal", "below",
#' or "above" for each derivative in ps.
#' @return A matrix of constraints.

get_edge_constraints = function(ms, s, directions = NULL) {
  # If m is equal to zero at all places, the return value is NULL.
  if(all(ms == 0) | all(ps == -1) | is.null(s) | length(s) == 0) {
    return(NULL)
  }

  K = length(s) + 1      # K is the number of bins in the polygram.
  s_aug = c(0, s ,1)     # The support is always assumed to be (0, 1).

  weights = sapply(1:(K+1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  t(sapply(0:ps[split_index], function(rho) {
    m = ms[split_index]

    nus = seq(0,rho)
    signs = (-1)^(nus + rho)
    binoms = choose(rho, nus)

    lower = c(rep(0,ms[split_index]-rho),
                weights[split_index+1]^(rho+1)*rev(binoms)*signs)
    lower*prod((m-rho+1):(m+1))
    }))

  list(constraint = constraint_matrix,
       lower = lower,
       upper = upper)
}


#' Returns the constraints needed to enforce symmetry of a Bernstein polygram.
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points. Must be symmetric for the procedure to work.
#' @return A constraint matrix.

get_symmetry_constraints = function(ms, s, symmetric = FALSE) {

  if(!symmetric) return(NULL)

  # Tests if the supplied s is symmetric.
  k = length(s)
  s_aug = c(0, s, 1)
  lengths = sapply(1:(k+1), function(j) s_aug[j+1] - s_aug[j])
  diffs = sapply(1:ceiling(k/2), function(j) lengths[k+2-j] - lengths[j])

  tol = 10^-8

  if(abs(sum(diffs)) > tol) {
    warning("The supplied s does not appear symmetric, a requirement for
            the symmetry constraints to work.")
  }

  # Tests if the supplied ms are symmetric:
  diffs = sapply(1:ceiling(k/2), function(j) ms[k+2-j] - ms[j])

  tol = 10^-8

  if(abs(sum(diffs)) > tol) {
    warning("The supplied m does not appear symmetric, a requirement for
            the symmetry constraints to work.")
  }

  len = floor(sum(ms+1)/2)
  A = diag(len)
  B = -apply(A, 1, rev)

  if (sum(ms+1) %% 2 == 0) {
    constraint = cbind(A,B)
  } else {
    constraint = cbind(A, rep(0,len) ,B)
  }

  #
  # # Now can make the symmetries.
  # if (k %% 2 == 1 | m %% 2 == 1) {
  #   A = diag(sum(ms[1:floor(k/2)]+1))
  #   B = -apply(A, 1, rev)
  #   constraint = cbind(A,B)
  # } else {
  #   A = diag(((m+1)*(k+1)+1)/2-1)
  #   B = -t(apply(A, 1, rev))
  #   constraint = cbind(A, rep(0, nrow(B)), B)
  # }

  bound = rep(0, nrow(constraint))
  list(constraint = constraint,
       lower      = bound,
       upper      = bound)
}

#' Returns the constraints needed to enforce monotonicity of a Bernstein
#' polygram.
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points.
#' @param monotone Either "decreasing" or "increasing". Can also be NULL,
#' but that should have been caught be now.
#' @return A list containing a constraint matrix and upper and lower bounds
#' for the monotone

get_monotonicity_constraints = function(ms, s, monotone) {
  # If p > -1, we won't need the sum constraints.
  if (is.null(monotone)) return(NULL)
  if (monotone != "increasing" & monotone != "decreasing") {
    stop("Only 'increasing', 'decreasing' and NULL are allowed as arguments
          to 'monotone'.")
  }

  # First we find the local constraints. These are constraints for each
  # Bernstein density in the mixture.
  matrices = lapply(ms, function(m) {
    if (m == 0) {
      return(0)
    } else {
    A = diag((m+1))
    B = cbind(rep(0, nrow(A)), -A)[,1:(ncol(A))]
    block = -(A + B)[1:(nrow(A)-1),]
    dim(block) = c(m, m+1)
    return(block)
    }
  })

  local_constraint  = do.call(hisemi::directSum, matrices)
  local_lower       = rep(0, nrow(local_constraint))
  local_upper       = local_lower + Inf

  # The global constraints are connectivity constraints: When the density is
  # increasing, the left-most point of the next density must be larger than or
  # equal to the the right most point of the current one.

  global = get_connectivity_constraints(ms = ms, s = s, ps = rep(0,length(s)),
                                        direction = "left")

  global_constraint  = global$constraint
  global_lower       = global$lower
  global_upper       = global$upper

  constraint = rbind(global_constraint, local_constraint)
  lower      = c(global_lower, local_lower)
  upper      = c(global_upper, local_upper)

  if(monotone == "decreasing") {
    constraint = -constraint
  }

  list(constraint = constraint,
       lower      = lower,
       upper      = upper)
}

#' Returns the constraints needed to enforce convexity / concavity of
#' a Bernstein histogram.
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points.
#' @param shape Can be "convex" or "concave" or NULL.
#' @return A constraint matrix.
get_shape_constraints = function(ms, s, shape) {

  if (is.null(shape)) return(NULL)
  if (shape != "convex" & shape != "concave") {
    stop("Only 'convex', 'concave' and NULL are allowed as arguments
         to 'shape'.")
  }

  if (all(ms == 0)) {
    # When m == 0, the only solution is the uniform distribution. That happens
    # when each weight equals the length its bin.
    K = length(s)
    s_aug = c(0, s, 1)
    lengths = sapply(1:(K + 1), function(i) s_aug[i+1] - s_aug[i])
    constraint = diag(K+1)
    upper = lengths
    lower = lengths
    return (list(constraint = constraint, lower = lower, upper = upper))
  }

  # First we find the local constraints. These are constraints for each
  # Bernstein density in the mixture.
  matrices = lapply(ms, function(m) {
    if (m == 0 ){
      return(0)
    } else if (m == 1) {
      block = c(0,0)
      dim(block) = c(1,2)
      return(block)
    } else {
      A = diag((m+1))
      B = cbind(rep(0, nrow(A)), -2*A)[,1:(ncol(A))]
      C = cbind(rep(0, nrow(A)),rep(0, nrow(A)), A)[,1:(ncol(A))]
      block = (A + B + C)[1:(nrow(A)-2),]
      dim(block) = c(m-1,m+1)
      return(block)
    }
  })

  local_constraint  = do.call(hisemi::directSum, matrices)
  local_lower       = rep(0, nrow(local_constraint))
  local_upper       = local_lower + Inf

  # The global constraints are connectivity constraints: When the density is
  # convex, the left-most derivative of the next density must be larger than or
  # equal to the the right most derivative of the current one.

  global = get_connectivity_constraints(ms = ms, s = s, ps = rep(1,length(s)),
                                        direction = "right")

  global_constraint  = global$constraint
  global_lower       = global$lower
  global_upper       = global$upper

  constraint = rbind(global_constraint, local_constraint)
  lower      = c(global_lower, local_lower)
  upper      = c(global_upper, local_upper)

  if(shape == "concave") {
    constraint = -constraint
  }

  list(constraint = constraint,
       lower      = lower,
       upper      = upper)
}

#' Returns the constraints needed to enforce the pth moment to be in
#' the interval [a,b]
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points.
#' @param moment The moment desired.
#' @param a Lower bound for the mean
#' @param b Upper bound for the mean. Defaults to a.
#' @return A constraint matrix.

get_moment_constraints = function(m, s, moments, lower_bounds, upper_bounds) {

  # Add checks here.

  if(is.null(upper_bounds) & !is.null(lower_bounds)) {
    upper_bounds = lower_bounds
  }

  if(!is.null(upper_bounds) & is.null(lower_bounds)) {
    lower_bounds = upper_bounds
  }

  if(is.null(moments) | (is.null(upper_bounds) & is.null(lower_bounds))) {
    return(NULL)
  }

  if(!all(m==m[1])) {
    stop("Moment conditions only supported for equal ms.")
  }

  len        = length(moments)
  constraint = NULL
  lower      = NULL
  upper      = NULL

  s_aug   = c(0, s, 1)
  K       = length(s_aug) - 1
  s_lower = c(0, s)
  weights = sapply(1:(length(s_aug)-1), function(i) s_aug[i+1] - s_aug[i])

  for (index in 1:len) {

    lower_bound = lower_bounds[index]
    upper_bound = upper_bounds[index]
    moment = moments[index]

    S_matrix = outer(0:moment, 1:K, FUN =
                     function(i,j) choose(moment,i)*weights[j]^i*s_lower[j]^(moment - i))
    PI_matrix = outer(0:m, 1:moment, FUN =
                      Vectorize(function(nu, r) prod((nu + 1:r)/(m + 1 + 1:r))))
    PI_matrix = cbind(rep(1,m+1), PI_matrix)
    Sigma_matrix = (PI_matrix %*% S_matrix)
    dim(Sigma_matrix) = length(Sigma_matrix)

    constraint = rbind(constraint, Sigma_matrix)
    lower = c(lower, lower_bound)
    upper = c(upper, upper_bound)
  }

  return(list(constraint = constraint, lower = lower, upper = upper))
}

#' Returns the constraints needed to enforce the qth quantile to be between
#' the splits s_i and s_k
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points.
#' @param quantiles The quantiles desired. Should be an increasing vector in (0,1).
#' @param lower_indices Increasing vector of lower indices for the involved split
#' points.
#' @param upper_indices Increasing vector of lower indices for the involved split
#' points.
#' @return A list containting a constraint matrix and two vectors of lower and
#' upper bounds.

get_quantile_constraints = function(m, s, quantiles, lower_indices, upper_indices) {

  # Add checks here.

  if(is.null(upper_indices) & !is.null(lower_indices)) {
    upper_indices = lower_indices
  }

  if(!is.null(upper_indices) & is.null(lower_indices)) {
    lower_indices = upper_indices
  }

  if(is.null(quantiles) | (is.null(upper_indices) & is.null(lower_indices))) {
    return(NULL)
  }

  if(!all(m==m[1])) {
    stop("Quantile conditions only supported for equal ms.")
  }

  len = length(quantiles)
  constraint = NULL
  lower      = NULL
  upper      = NULL

  s_aug   = c(0, s, 1)
  K       = length(s_aug) - 1

  for(index in 1:len) {
    i = lower_indices[index]
    j = upper_indices[index]
    q = quantiles[index]


    is = rep(0,K)
    is[1:(i)] = 1

    js = rep(0,K)
    js[1:(j)] = 1

    row_one = sapply(is, function(x) rep(x,m+1))
    dim(row_one) = length(row_one)

    row_two = sapply(js, function(x) rep(x,m+1))
    dim(row_two) = length(row_two)

    constraint = rbind(constraint, row_one, row_two)
    lower = c(lower, 0, q)
    upper = c(upper, q, 1)
  }

  return(list(constraint = constraint, lower = lower, upper = upper))
}


#' Returns the constraints needed to enforce modality at (nu,k).
#'
#' @param m The order of the Bernstein densities.
#' @param s The split points.
#' @param p Order of connectedness.
#' @param nu nu-index of the mode.
#' @param k k-index of the mode
#' @return A constraint matrix.

get_unimodal_constraints = function(m, s, p, unimodal, nu, k) {
  if(!unimodal) return(NULL)
  if(is.null(nu) | is.null(k)) return(NULL)
  # If p > -1, we won't need the sum constraints.
  if (p > -1 | length(s) == 0) {
    if (nu == 0 & k == 1) {
      return(get_monotonicity_constraints(m, s, p, monotone = "decreasing"))
    }

    if (nu == m & k == length(s) + 1) {
      return(get_monotonicity_constraints(m, s, p, monotone = "increasing"))
    }

    len = length(s) + 1
    # The constraints building up to the mode. They start on the first item
    # and continue all the way to the ((m+1)*(k-1) + nu)th item.
    A = -diag((m+1)*(k-1)+nu+1)
    B = cbind(rep(0, nrow(A)), -A)[,1:(ncol(A))]
    block = (A + B )[1:(nrow(A)-1),]
    dim(block) = dim(A) - c(1,0)

    #dim(block) = c(m-1,m+1)

    # The constraints building down from the mode. They start on the
    # ((m+1)*(k-1) + nu + 1)th item and continue all the way to last one at
    # (m+1)*len

    A_ = diag((m+1)*len - dim(A)[1] + 1)
    B_ = cbind(rep(0, nrow(A_)), -A_)[,1:(ncol(A_))]
    block_ = (A_ + B_)[1:(nrow(A_)-1),]
    dim(block_) = dim(A_) - c(1,0)

    # The constraints are constructed, but have duplicate columns in
    # middle: The column with index (m+1)*(k-1) + nu + 1 and
    # (m+1)*(k-1) + nu + 2 must be summed.
    constraints_ = hisemi::directSum(block, block_)
    constraints_[,(m+1)*(k-1) + nu + 1] = constraints_[,(m+1)*(k-1) + nu + 1] +
      constraints_[,(m+1)*(k-1) + nu + 2]
    constraint = constraints_[,-((m+1)*(k-1) + nu + 2)]

    upper = rep(1, nrow(constraint))
    lower = 0*upper
    list(constraint = constraint, lower = lower, upper = upper)

  } else {
    stop("Modality constraints only supported for connected densities.")
  }
}



