#' Performs k-fold cross-validation on the input formula.
#'
#' @export
#' @param formula The formula to work on, same semantics as for \code{polygram}.
#' Raw data can also be supplied.
#' @param nfolds Support of the data. Defaults to [0,1].
#' @param candidates a list of candidate polygrams. Each element of the list
#' must contain a vector (or integer) of splits, a vector m of Bernstein
#' degrees, and a vector p of derivative orders. If an integer N is supplied,
#' it takes the candidate list to be s = 0:N and m = 3, p = 2.
#' @param data an optional vector of data to be used in the formula.
#' @param m optional list of ms, only relevant if candidates includes only
#' splits.
#' @param p optional list of ps, only relevant if candidates includes only
#' splits.
#' @param ... additional arguments passed to \code{polygram}.
#' @return A \code{cv.polygram} object.

cv.polygrams = function(formula, nfolds = 10, candidates = NULL,
                        support = NULL, data = NULL, m = NULL, p = NULL, ...) {

  if(is.numeric(candidates) | is.integer(candidates)){

    if(length(candidates) == 1) {
      candidates = 0:ceiling(candidates)
    } else if(length(candidates) == 0) {
      stop("candidates cannot be empty.")
    }

    if(is.null(m)) {
      m = as.list(rep(3, length(candidates)))
    } else if (length(m) == 1) {
      m = as.list(rep(m, length(candidates)))
    }

    if(is.null(p)) {
      p = lapply(m, function(m) m - 1)
    } else if (length(p) == 1) {
      p = as.list(rep(p, length(candidates)))
    }

    candidates = as.list(candidates)
    for(i in 1:length(candidates)) {
      candidates[[i]] = list(s = candidates[[i]], m = m[[i]], p = p[[i]])
    }

  }

  if("formula" %in% class(formula)) {
    # Stashes all the data into a separate class.
    response = head(all.vars(formula), n = 1)
    data_ = eval(parse(text = response),
                envir = attr(formula, ".Environment"))

  } else {
    tryCatch({
      data_ = as.numeric(formula)
      formula = formula ~ NULL
    }, error = function(e) {
      print("Supply either a valid formula or some valid (vector) data.")
    })
  }

  if(!is.null(data)) {
    data_ = data
  }
  if(!is.null(data_) & (NA %in% data_)) {
    data_ = data_[!is.na(data_)]
  }

  n_data_ = length(data_)
  foldid = sample(rep(seq(nfolds), length = n_data_))
  cross_validates = sapply(candidates, function(param_list) {
    s = param_list$s
    m = param_list$m
    p = param_list$p

    results = sapply(seq(nfolds), function(i) {
      tryCatch({
      current_indices = (foldid == i)
      obj = polygram(formula, data = data_[!current_indices], s = s,
                     m = m, p = p, support = support)
      pen = penalty(obj) - 2*mean(dpolygram(data_[current_indices], obj))
      pen
      }, error = function(e) NA)
    })

    results
  })

  results = colMeans(cross_validates, na.rm = TRUE)
  #attr(results, "vars") = apply(cross_validates, 2, var)
  attr(results, "candidates") = candidates
  attr(results, "nfolds") = nfolds
  attr(results, "support") = support
  attr(results, "ellipses") = list(...)
  attr(results, "class") = c("cv.polygram","numeric")

  results
}

print.cv.polygram = function(x, ...) {
  cat(x)
}

cv.polygrams_select = function(formula, nfolds = 10, candidates = NULL,
                               support = NULL, data = NULL, m = NULL, p = NULL, ...) {
  models = cv.polygrams(formula, nfolds = nfolds, candidates = candidates,
            support = support, data = data, m = m, p = p, ...)
  selected = attr(models,"candidates")[[which.min(models)]]

  polygram_object = polygram(formula, s = selected$s, m = selected$m,
                             p = selected$p, support = support,
                             data = data, ...)
  return(polygram_object)
}
