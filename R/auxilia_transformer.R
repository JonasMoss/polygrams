#' Changes the split vector \code{s} of a polygram with the supplied
#' \code{s_new}.
#'
#' @export
#' @param polygram_object The supplied polygram_object
#' @param s_new The new vector of split points. Must satisfy \code{s} $\subset$
#' \code{s_new}.
#' @return A \code{polygram} object.
#' @details Returns an equivalent polygram with the new split points. The new
#' polygram is equivalent to the old one in the sense that the distribution
#' functions are equal. It does not attempt to find inverse solutions, so
#' \code{s} $\subset$ \code{s_new} is required. Furthermore, this inclusion
#' must be exact. (No different floating points arithmetics allowed.)
#'
#' @examples

polygram_transform = function(polygram_object, s_new) {
  # First we test if s %in% s_new
}
