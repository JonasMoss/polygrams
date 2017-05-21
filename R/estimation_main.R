#' Estimates a polygram with supplied \code{data}, \code{m}, \code{p} and split points. Several
#' hard constraints are supported.
#' @export
#' @param formula The data. Can be supplied as vector of values or as a \code{formula} object,
#' see \code{details} for how to use the formula.
#' @param s The vector of split points. If a natural number \code{K}, \code{s} is chosen as
#' \code{(1:K)/(K+1)}. Defaults to \code{length(data)^(1/3)}. If \code{0} or \code{NULL}, a Bernstein
#' density is computed.
#' @param m The degree of the Bernstein polynomial. Can be either a number or a vector. If a number, it
#' is repeated |s|+1 times; if a vector, it must be of length |s|+1, where the kth element
#' specifies the degree of the Bernstein polynomial at the kth bin. Defaults to 3.
#' @param p Connectedness order.  Can be either a number or a vector. If a number, it
#' is repeated |s| times; if a vector, it must be of length |s|, where the kth element
#' specifies the degree of connectedness at the kth split. Defaults to m-1.
#' @param support A vector \code{c(l, u)} specifying the support of the density.
#' Defaults to \code{c(0, 1)} provided it contains the support of the data. Use
#' this option when you know something about the support of your data.
#' @param method Optimisation method. "quadprog" (default) and "mosek" supported. Note that mosek requires
#' a license.
#' @return A \code{polygram} object.
#' @details The \code{formula} argument takes formulas such as \code{x ~ convex + decreasing},
#' with the semantic that the data is on the left hand side and constraints on the right hand side.
#' Supported constraints are c("increasing", "decreasing", "symmetric", "concave", "convex"), but
#' this list will be larger soon.
#' @examples
#'
#' # A convex and decreasing density.
#' x = datasets::sunspot.month
#' plot(polygram(x ~ convex + decreasing, support = c(0, 270), s = 10))
#'
#' # A crazy density.
#' set.seed(1984)
#' x = runif(11)
#' m = c(3,1,3,3,2,2,3,3,1,3)
#' p = c(1,1,2,1,0,1,2,1,1)
#' obj = polygram(x ~ symmetric, m = m, p = p)
#' plot(obj, bins = FALSE, main = "Symmetric polygram")

polygram = function(formula, s = NULL, m = NULL, p = NULL, support = NULL,
                    method = "quadprog") {

  ## ==========================================================================
  ## This first section handles checks and fills in default values.
  ## ==========================================================================

  if("formula" %in% class(formula)) {
    # Stashes all the data into a separate class.
    response = head(all.vars(formula), n = 1)
    data = eval(parse(text = response),
         envir = attr(formula, ".Environment"))

  } else {
    tryCatch({
    data = as.numeric(formula)
    formula = formula ~ NULL
    }, error = function(e) {
      print("Supply either a valid formula or some valid (vector) data.")
    })
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

  # Here we manipulate the support and the data. The data is molded into the
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

  ## ==========================================================================
  ## The second objective is to fill in the constraints. This is done by
  ## finding the appropriate terms in the formula object.
  ## ==========================================================================

  constraint_list = formulas_to_constraint_list(formula, ms, s, ps)
  constraint_matrix = constraint_list$constraint
  lower_bounds      = constraint_list$lower
  upper_bounds      = constraint_list$upper

  qobj = polygram_objective_matrix(ms, s)        # The objective matrix.
  dvec = -polygram_objective_vector(data, ms, s) # The objective vector.

  if (method == "mosek" | method == "rmosek" | method == "mosek") {

    if (!requireNamespace("Rmosek", quietly = TRUE)) {
      stop("Rmosek needed for this function to work. Please install it.",
           call. = FALSE)
    }

    qobj = as.mosek_mat(qobj)

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
    Dmat = qobj
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

