#' Estimates a polygram with supplied data, m, p and split points. Several
#' hard constraints are supported. The interface is prone to changes.
#' @export
#' @param data The data. Can be \code{NULL} if a density \code{d} is supplied.
#' @param s The vector of split points. If a natural number \code{K}, \code{s} is chosen as
#' \code{(1:K)/(K+1)}. Defaults to \code{length(data)^(1/3)}. If \code{0} or \code{NULL}, Bernstein
#' density is computed.
#' @param support A vector \code{c(l, u)} specifying the support of the density. Defaults to \code{c(0, 1)}.
#' @param m The order of the Bernstein polynomial. Defaults to 3.
#' @param p Connectedness order. Defaults to \code{m-1}.
#' @param d Optional distribution function. Can be the base function in a
#' Dirichlet process if used in conjunction with M. Defaults to NULL. Currently not supported.
#' @param M The concentration parameter in a Dirichlet process. Defaults to 1. Ignored
#' if either \code{d} or \code{data} is \code{NULL}.
#' @param symmetric Hard symmetry constraint. Requires that s is symmetric as well.
#' @param monotone Hard monotonicity constraint. Can be \code{"increasing"},
#' \code{"decreasing"} or \code{NULL}.
#' @param shape Hard convexity / concavity constraint. Can be \code{"convex"},
#' \code{"concave"} or \code{NULL}.
#' @param moment_conditions A list of moment conditons. The argument "moment" tells
#' which moment is being constrained, while a and be give upper and lower bounds.
#' If only one a or b is given, it is assumed that a = b and the condition is
#' constant.
#' @param quantile_conditions A list of quantile conditions. List of "quantile", i and j specifing
#' the indices of the splits for where the quantile is.
#' @param unimodal Hard unimodality constraint. Not supported yet.
#' @param lower_boundary Hard lower boundary constraint. Defualts to NULL.
#' @param upper_boundary Hard upper boundary constraint. Defualts to NULL.
#' @param nu,k Specific unimodality constraints.
#' @param lambdas A vector of lambdas. The first element is lambda1, second lambda2, etc.
#' @param verbose An integer between 0 and 10. Modifies Rmoseks verbose command.
#' @param method Optimisation method. "quadprog" (default) and "mosek" supported. Note that mosek requires
#' a license.
#' @return A polygram object.

polygram = function(data, s = NULL, support = NULL, m = NULL, p = NULL, d = NULL,
                    M = 1, symmetric = FALSE, monotone = NULL, shape = NULL, moment_conditions = NULL,
                    unimodal = FALSE, lower_boundary = NULL, quantile_conditions = NULL, lambdas = NULL,
                    upper_boundary = NULL, nu = NULL, k = NULL, verbose = 0,
                    method = "QP") {

  if(unimodal) {
    stop("Unimodality is not implemented yet.")
  }
  # Fill in standard values
  if (is.null(m)) m = 3
  if (is.null(p)) p = m-1

  if (length(m) != 1 & is.null(s)) {
    len = length(m) - 1
    s = (1:len)/(len+1)
  } else if (is.null(s)) {
    s = floor(length(data)^(1/3))
  }
  #if (p >= m) stop("p must be smaller than m.")

  # Her we manipulate the support and the data. The data is molded into the
  # [0,1]-paradigm by dividing removing the lower support and dividing by the
  # length of the support. The same is done with the vector of splits.

  if(!is.null(data) & (NA %in% data)) {
      data = data[!is.na(data)]
  }

  if(is.null(support)) {
    if(is.null(data)) {
      support = c(0,1)
    } else {
      support = c(min(c(0, data)), max(1, data))
    }
  }

  if(support[1] > support[2]) {
    warning(paste0("The support vector appears to be flipped: c(",support[1],support[2],").
                   We will flip it back for you."))
    support = c(support[2], support[1])

  }

  if(!is.null(data)) {

    if (NA %in% data) {
      data = data[!is.na(data)]
    }

    if(max(data) > support[2]) {
      warning("Support does not contain the max data. Support adjust accordingly.")
      support[2] = max(data)
    } else if(min(data) < support[1]) {
      warning("Support does not contain the min data. Support adjust accordingly.")
      support[1] = min(data)
    }
    data  = (data - support[1])/(support[2] - support[1])
  }

  if (length(s) == 1 & !is.null(s)) {
    if ( (ceiling(s) == s) ) {
      if (s == 0) {
        s = NULL
        p = -1
      } else {
        s = as.integer(s)
        s = (1:s)/(s+1)
        s = s*(support[2] - support[1]) + support[1]
      }
    }
  } else if (is.null(s)) {
    p = -1
  } else if (!("numeric" %in% class(s) | "integer" %in% class(s))){
    stop("s must either be numeric, null or an integer.")
  }

  s_org = s
  s     = (s - support[1])/(support[2] - support[1])

  K = length(s) + 1


  # This is done in order to support variable ms and ps. Usually, polygram will be called with
  # only a single value for m and p, but other values should be supported as well.

  if(!is.numeric(m) & !is.integer(m)) {
    stop("m must be numeric or integer.")
  }

  if(length(m) == 1) {
    ms = rep(m, K)
  } else {
    ms = m
    if(length(ms) != (length(s) + 1)) {
      stop("When length(m)!=1, length(m) = (length(s) + 1) is
           required.")
    }
  }

  if(!is.numeric(p) & !is.integer(p)) {
    stop("p must be numeric or integer.")
  }

  if(length(p) == 1) {
    ps = rep(p, K-1)
  } else {
    ps = p
    if(length(ps) != length(s)) {
      stop("When length(p)!=1, length(p) = length(s) is
           required.")
    }
  }

  constraint_list = list(
    sum_constraint           = list(constraint = rep(1,sum(ms+1)),
                                                 lower      = 1,
                                                 upper      = 1),

    connectivity_constraints = get_connectivity_constraints(ms = ms, s = s, ps = ps,
                                    directions = NULL),

    symmetry_constraints     = get_symmetry_constraints(ms = ms, s = s,
                                    symmetric = symmetric),

    monotone_constraints     = get_monotonicity_constraints(ms = ms, s = s,
                                    monotone = monotone),

    shape_constraints        = get_shape_constraints(ms = ms, s = s,
                                    shape = shape),

    moment_constraints       = get_moment_constraints(m = m, s = s,
                                    moments       = moment_conditions$moments,
                                    lower_bounds  = moment_conditions$lower_bounds,
                                    upper_bounds  = moment_conditions$upper_bounds),

    quantile_constraints     = get_quantile_constraints(m = m, s = s,
                                    quantiles     = quantile_conditions$quantiles,
                                    lower_indices = quantile_conditions$lower_indices,
                                    upper_indices = quantile_conditions$upper_indices),

    unimodal_constraints     = get_unimodal_constraints(m = m, s = s, p = p,
                                    unimodal = unimodal,
                                    nu = nu,
                                    k = K)
  )

  constraint_matrix = do.call(rbind, sapply(constraint_list, function(object) object$constraint))
  lower_bounds      = do.call(c,     sapply(constraint_list, function(object) object$lower))
  upper_bounds      = do.call(c,     sapply(constraint_list, function(object) object$upper))

  # This is the basic objective matrix - without any additional penalties.
  qobj_ = polygram_objective_matrix(ms, s)

  # Here comes additional penalties from the lambda list.
  if (!is.null(lambdas)) {
    for(i in 1:(length(lambdas))) {
      qobj_ = qobj_ + lambdas[i]*polygram_penalty_matrix(s, m, p = i)
    }
  }

  dvec = -polygram_objective_vector(data, ms, s) # The objective vector.

  if (method == "mosek" | method == "rmosek" | method == "mosek") {
    qobj = as.mosek_mat(qobj_)

    N = length(dvec)

    # These are the box constraints. The basic box constraints for the coefficients to be positive,
    # but they can be augmented by lower and upper boundaries - forced values at the edges of the
    # support region.

    blx = c(rep(0,N))
    bux = c(rep(1,N))

    if(!is.null(lower_boundary)) {
      blx[1] = lower_boundary
      bux[1] = lower_boundary
    }

    if(!is.null(upper_boundary)) {
      blx[length(bux)] = upper_boundary
      bux[length(bux)] = upper_boundary
    }

    qo1 <- list()
    qo1$sense <- "min"
    qo1$c <- dvec
    qo1$A <- Matrix::Matrix(constraint_matrix, sparse=TRUE)
    qo1$bc <- rbind(blc = lower_bounds,
                    buc = upper_bounds)

    # Box constraints. blx are the lower bounds for the weights, bux upper bounds.
    # blx is only 0 if no lower_bound is supplied, and similarily for bux.

    qo1$bx <- rbind(blx = blx,
                  bux = bux)

    qo1$qobj <- qobj

    r = Rmosek::mosek(qo1, opts = list(soldetail = 1,
                             verbose = verbose))
    v = r$sol$itr$xx
    if(r$sol$itr$prosta == "ILL_POSED") {
      stop("The problem is ill-posed! Try a different s.")
    }

    if(r$sol$itr$solsta != "OPTIMAL" & r$sol$itr$solsta != "NEAR_OPTIMAL") {
      warning("The solution is not optimal or near optimal.")
    }

    v = split(v, rep(1:(length(ms)), ms + 1))
    names(v) = NULL
    attr(v, "loss") = 2*r$sol$itr$pobjval
    attr(v, "m") = ms
    attr(v, "K") = K

    if (is.null(s_org)) {
      attr(v, "s") = rep(1,0)
    } else {
      attr(v, "s") = s_org
    }

    attr(v, "p") = ps
    attr(v, "support") = support
    attr(v, "call") = sys.call()
    #attr(v, "sol") = r$sol
    class(v) = c("polygram")
    v
  } else if (method == "quadprog" | method == "QP" | method == "qp") {
    qp_list = mosek2qp(constraint_matrix, lower_bounds, upper_bounds)
    Dmat = qobj_
    Amat = qp_list$Amat
    bvec = qp_list$bvec
    meq  = qp_list$meq

    # In addition, quadprog needs positivity constraints:
    bvec = c(bvec, positive = rep(0,sum(ms+1)))
    Amat = rbind(Amat, positive = diag(sum(ms+1)))

    r = quadprog::solve.QP(Dmat = Dmat, dvec = -dvec, Amat = t(Amat), bvec = bvec, meq = meq)
    v = r$sol

    v = split(v, rep(1:(length(ms)), ms + 1))
    names(v) = NULL
    attr(v, "loss") = 2*r$value
    attr(v, "m") = ms
    attr(v, "K") = K

    if (is.null(s_org)) {
      attr(v, "s") = rep(1,0)
    } else {
      attr(v, "s") = s_org
    }

    attr(v, "p") = ps
    attr(v, "support") = support
    attr(v, "method") = method
    attr(v, "call") = sys.call()
    class(v) = c("polygram")
    v

    # list(v = v,
    #      Amat = Amat,
    #      bvec = bvec,
    #      meq  = meq)
  }
}

#' Makes quadprog-compliant matrices from rmosek-compliant ones.
#' @param constraint_matrix matrix of constraints.
#' @param lower_bounds lower bounds for the constraints.
#' @param upper_bounds upper bounds for the constraints.
#' @param tolerance tolerance for equality checking. Defaults to 10^-8.
#' @return A list containing a new constraints matrix "Amat", a vector "bvec" holding lower bounds
#' and a number meq telling how many are equality constraints.
mosek2qp  = function(constraint_matrix, lower_bounds, upper_bounds, tolerance = 10^-8) {

  # First we isolate the equality constraints.
  meq_indices = sapply(1:length(lower_bounds), function(i) isTRUE(all.equal(
    lower_bounds[i], upper_bounds[i],tolerance = tolerance)))
  bvec = lower_bounds[meq_indices]
  meq  = sum(meq_indices)
  Amat = constraint_matrix[meq_indices,]

  # Now we must isolate the lower and upper bounds.
  nmeq_indices = !meq_indices
  nmeq_mat = constraint_matrix[!meq_indices,]

  # From this we obtain the lower bounds:
  lower_indices = which(lower_bounds[nmeq_indices] != -Inf)
  bvec = c(bvec, (lower_bounds[nmeq_indices])[lower_indices])
  Amat = rbind(Amat, nmeq_mat[lower_indices,])

  # And the upper bounds:
  upper_indices = which(upper_bounds[nmeq_indices] != Inf)
  bvec = c(bvec, -upper_bounds[nmeq_indices][upper_indices])
  Amat = rbind(Amat, -nmeq_mat[upper_indices,])
  list(Amat = Amat, bvec = bvec, meq = meq)
}

