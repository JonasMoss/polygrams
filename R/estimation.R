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
#' @param data an optional vector containing the data used in the formula.
#' @param nfold the number of folds in cross-validation. Defaults to 5.
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
                    data = NULL,
                    method = c("quadprog", "Rmosek", "constrOptim")) {

  ## ==========================================================================
  ## This first section handles checks and fills in default values.
  ## ==========================================================================

  method = match.arg(method)
  if("formula" %in% class(formula)) {
    # Stashes all the data into a separate class.
    response = head(all.vars(formula), n = 1)
    data_ = eval(parse(text = response),
         envir = attr(formula, ".Environment"))
    flag_formula = TRUE

  } else {
    tryCatch({
    data_ = as.numeric(formula)
    formula = formula ~ NULL
    flag_formula = FALSE
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
    s = floor(length(data_)^(1/3))
  }

  #if (p >= m) stop("p must be smaller than m.")

  # Here we manipulate the support and the data. The data is molded into the
  # [0,1]-paradigm by dividing removing the lower support and dividing by the
  # length of the support. The same is done with the vector of splits.

  if(!is.null(data)) {
    data_ = data
  }

  if(!is.null(data_) & (NA %in% data_)) {
      data_ = data_[!is.na(data_)]
  }


  if(is.null(support)) {
    if(is.null(data_)) {
      support = c(0,1)
    } else {
      support = c(min(c(0, data_)), max(1, data_))
    }
  }


  if(support[1] > support[2]) {
    warning(paste0("The support vector appears to be flipped: c(",support[1],support[2],").
                   We will flip it back for you."))
    support = c(support[2], support[1])

  }

  if(!is.null(data_)) {

    if (NA %in% data_) {
      data_ = data_[!is.na(data_)]
    }

    if(max(data_) > support[2]) {
      warning("Support does not contain the max data_. Support adjust accordingly.")
      support[2] = max(data_)
    } else if(min(data_) < support[1]) {
      warning("Support does not contain the min data_. Support adjust accordingly.")
      support[1] = min(data_)
    }
    data_  = (data_ - support[1])/(support[2] - support[1])

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
  dvec = -polygram_objective_vector(data_, ms, s, c(0, 1)) # The objective vector.

  if (method == "Rmosek") {

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

    r = Rmosek::mosek(qo1, opts = list(soldetail = 1))
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

  } else if (method == "quadprog") {

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
  } else if (method == "constrOptim") {

    # The linear constrained are obtained from the quadprog procedures.
    qp_list = mosek2qp(constraint_matrix, lower_bounds, upper_bounds)
    Dmat = qobj
    Amat = qp_list$Amat
    bvec = qp_list$bvec
    meq  = qp_list$meq

    # The starting value is obtained by running a quadprog iteration

    # In addition, quadprog needs positivity constraints:
    bvec = c(bvec, positive = rep(0,sum(ms+1)))
    Amat = rbind(Amat, positive = diag(sum(ms+1)))

    r = quadprog::solve.QP(Dmat = Dmat, dvec = -dvec, Amat = t(Amat), bvec = bvec, meq = meq)
    v = r$sol

    v = split(v, rep(1:(length(ms)), ms + 1))
    names(v) = NULL
    attr(v, "loss") = 2*r$value
  }

  class(v) = c("polygram")
  attr(v, "m") = ms
  attr(v, "K") = K
  attr(v, "s") = if (is.null(s_org)) rep(1,0) else s_org
  attr(v, "p") = ps
  attr(v, "support") = support
  attr(v, "method") = method
  attr(v, "call") = sys.call()
  v

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

#' Calculates the penalty matrix for s, m and p.
#'
#' @param s A vector of splits. Sould not include end points.
#' @param m The orders of the Bernstein polynomial.
#' @param p Derivative order.
#' @return A matrix

polygram_penalty_matrix = function(s, m, p) {

  K = length(s) + 1

  # The penalty matrix is a direct sum of K square matrices
  # of size m[k] x m[k]. Each of these matrices multplied with
  # a scalar after production.

  if(length(m) == 1) {
    m = rep(m, K)
  }

  s_aug = c(0, s, 1)

  Sigma = Vectorize(function(nu, eta, mk, p) {
    lower_bound_nu  = max(0, nu + p - mk)
    lower_bound_eta = max(0, eta + p - mk)
    upper_bound_nu  = min(nu, p)
    upper_bound_eta = min(eta, p)

    if ((lower_bound_nu  > upper_bound_nu) |
        (lower_bound_eta > upper_bound_eta)) {
      return(0)
    }

    is = lower_bound_nu:upper_bound_nu
    js = lower_bound_eta:upper_bound_eta

    summands = function(i, j) {
      first  = (-1)^(j+i)*choose(mk - p, nu - i)*choose(mk - p, eta - j)
      second = choose(p, i)*choose(p, j)
      third  = choose(2*(mk-p), nu - i + eta -j)
      first*second/third
    }

    sum(outer(is, js, FUN = summands))
  })

  matrix_list = lapply(1:K, function(k) {
    mk = m[k]
    wk = s_aug[k+1] - s_aug[k]
    mk_weight_a = (choose(mk, p)*factorial(p))^2*(mk - p + 1)^2
    mk_weight_b = (2*(mk - p) + 1)*wk^(2*p+1)
    mk_matrix = outer(0:mk, 0:mk, Sigma, mk = mk, p = p)
    mk_matrix*mk_weight_a/mk_weight_b
  })


  A = do.call(direct_sum, matrix_list)

  return(A)
}

#' Transforms a symmetric matrix into the form demanded by Rmosek.
#'
#' @param mat A matrix.
#' @return A matrix in triangular sparse triplet form.

as.mosek_mat = function(mat) {

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix needed for this function to work. Please install it.",
         call. = FALSE)
  }

  new = mat*0
  new[upper.tri(mat, diag = TRUE)] = mat[upper.tri(mat, diag = TRUE)]
  new = as(new, "dtTMatrix")
  mat_ = list()
  mat_$i = attr(new, "j") + 1
  mat_$j = attr(new, "i") + 1
  mat_$v = attr(new, "x")
  mat_
}


#' Calculates the objective matrix for use in a polygram.
#'
#' @param ms Vector of Bernstein orders.
#' @param s The vector of splits.

polygram_objective_matrix = function(ms, s) {
  # The objective matrix is a weighted direct sum of |s| + 1
  # matrices on the same form: A_

  A_ = function(m) {

    FUN = Vectorize(function(eta,nu) {
      (m+1)/(choose(2*m+1,m+1))*choose(eta+nu,eta)*choose(2*m-eta-nu,m-eta)
    })

    outer(0:m, 0:m, FUN)
  }

  # The next step is to do the direct sum of the weighted matrices.
  s_aug = c(0,s,1)

  weights = 1/sapply(1:(length(s_aug) - 1), function(i) {
    s_aug[i+1] - s_aug[i]
  })

  K = length(s) + 1

  A = do.call(direct_sum, lapply(1:K, function(k) weights[k]*A_(ms[k])))

  return(A)
}

