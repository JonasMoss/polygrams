#' A helper function for \code{polygram}. Creates a list of constraints out of
#' a \code{formula}, \code{ms}, \code{ps} and \code{s}.
#' @param form a formula
#' @param ms segrees of Bernstein densities
#' @param s split points
#' @param ps connectivity constraints
#' @return A list containing a 1.) a constraint matrix, 2.) lower bounds and 3.)
#' upper bounds.
#' @details
#' Currently supported options are:
#'\enumerate{
#'  \item "convex", "concave"; both together produces straight line.
#'  \item "decreasing", "increasing"; both together produces uniform.
#'  \item "symmetric"
#'  \item mean: Takes two "arguments", one in \code{c("=", ">=", "=<")} and a numeric
#'            in support.
#'  \item moment: Takes three "arguments", one in \code{c("=", ">=", "=<")}, a numeric
#'            in support, and a natural specifying the moment.
#'  \item quantile: Takes three "arguments", one in \code{c("=", ">=", "=<")}, a numeric
#'  \item quantile: Takes three "arguments", one in \code{c("=", ">=", "=<")}, a numeric
#'            in support, and a float in (0,1) specifying the quantile.
#'  }

formulas_to_constraint_list = function(form, ms, s, ps) {
  ## First the simplex and connectivity constraints.
  simplex_constraint_list = list(constraint = rep(1, sum(ms + 1)),lower = 1,upper = 1)
  connectivity_constraint_list = get_connectivity_constraints(ms = ms, s = s, ps = ps,directions = NULL)

  constraint_list = list(simplex_constraint_list      = simplex_constraint_list,
                         connectivity_constraint_list = connectivity_constraint_list)

  terms = attr(terms(form), "term.labels")
  if(length(terms) != 0) {
    geq_indices = sapply(terms, function(term) grepl(">=", term))
    leq_indices = sapply(terms, function(term) grepl("<=", term))
    eq_indices  = sapply(terms, function(term) grepl("=", term)) & !geq_indices & !leq_indices
    no_indices  = !geq_indices & !leq_indices & !eq_indices

    ## Now we will handle the terms with no (in)equalities.
    simple_terms = terms[no_indices]
    simple_terms_constraint_list = simple_terms_to_constraint_list(simple_terms, ms, s, ps)

    ## Now the (in)equality constraints.

    ## This must be combined:
    constraint_list = c(constraint_list,
                        list(simple_terms_constraint_list = simple_terms_constraint_list))


  }

    constraint_matrix = do.call(rbind, lapply(constraint_list, function(object) object$constraint))
    lower_bounds      = do.call(c,     lapply(constraint_list, function(object) object$lower))
    upper_bounds      = do.call(c,     lapply(constraint_list, function(object) object$upper))

    return(list(constraint_matrix = constraint_matrix,
                lower_bounds      = lower_bounds,
                upper_bounds      = upper_bounds))
}

#' A helper function for \code{formulas_to_constraint_list }. Handles the simple terms.
#' @param options a list of options
#' @param ms segrees of Bernstein densities
#' @param s split points
#' @param ps connectivity constraints
#' @return A list containing a 1.) a constraint matrix, 2.) lower bounds and 3.)
#' upper bounds.

simple_terms_to_constraint_list = function(simple_terms, ms, s, ps) {
  alternatives = c("convex", "concave", "symmetric", "decreasing", "increasing")
  all_allowed = all(simple_terms %in% alternatives)

  if(!all_allowed) {
    error_string = paste(alternatives, collapse = ", ")
    error_string = paste0("A single constraint term must be one of: ", error_string)
    stop(error_string)
  }

  constraint_list = list()

  if("symmetric" %in% simple_terms) {
    constraint_list$symmetric = get_symmetry_constraints(ms = ms, s = s, symmetric = TRUE)
  }

  if("convex" %in% simple_terms) {
    constraint_list$convex = get_shape_constraints(ms = ms, s = s, shape = "convex")
  }

  if("concave" %in% simple_terms) {
    constraint_list$concave = get_shape_constraints(ms = ms, s = s, shape = "concave")
  }

  if("increasing" %in% simple_terms) {
    constraint_list$increasing = get_monotonicity_constraints(ms = ms, s = s, monotone = "increasing")
  }

  if("decreasing" %in% simple_terms) {
    constraint_list$decreasing = get_monotonicity_constraints(ms = ms, s = s, monotone = "decreasing")
  }


  constraint_matrix = do.call(rbind, lapply(constraint_list, function(object) object$constraint))
  lower_bounds      = do.call(c,     lapply(constraint_list, function(object) object$lower))
  upper_bounds      = do.call(c,     lapply(constraint_list, function(object) object$upper))

  return(list(constraint_matrix = constraint_matrix,
              lower_bounds      = lower_bounds,
              upper_bounds      = upper_bounds))
}
