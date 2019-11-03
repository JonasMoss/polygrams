#' @export
plot.polygram = function(polygram_object, derivative = 0, rug = TRUE, bins = TRUE, xlab = NULL,
                         ylab = NULL, ylim = NULL, main = NULL, type = NULL, bty = NULL,...) {

  support = attr(polygram_object, "support")
  s       = c(attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  s_aug_  = c(0,(s - support[1])/(support[2] - support[1]),1)
  xx = seq(0, 1, by = 0.001)
  xx = xx*(support[2]-support[1]) + support[1]

  K = attr(polygram_object, "K")
  m = attr(polygram_object, "m")
  p = attr(polygram_object, "p")

  if(!is.integer(derivative) & !is.numeric(derivative)) {
    stop("Option 'derivative' must be an integer greater than -2.")
  } else if (derivative < -1) {
    stop("Option 'derivative' must be an integer greater than -2.")
  }

  yy_hist = dxpolygram(xx, polygram_object, p = derivative)

  if(is.null(main)) {
    if(all(m == m[1])) {
      main = paste0("Bernstein polygram (m = ", m[1], ", K = ", K,")")
    } else {
      main = paste0("Bernstein polygram (max(m) = ", max(m), ", K = ", K,")")
    }
  }
  if(is.null(ylim)) ylim = c(min(c(yy_hist,0)), max(yy_hist))
  if(is.null(xlab)) xlab = "x"
  if(is.null(ylab)) ylab = "Density"
  if(is.null(type)) type = "l"
  if(is.null(bty)) bty = "l"

  if(any(p == -1) | any(m == 0)) {

    plot(xx, yy_hist, bty = bty, ylim = ylim, xlab = xlab, ylab = ylab,
         main = main, type = "n")


    diffs = sapply(1:(length(s_aug_)-1), function(i) s_aug_[i+1] - s_aug_[i])
    eps = min(diffs)/2

    for (i in 1:(length(s)+1)) {
      xx_ = seq(s_aug_[i], s_aug_[i+1] - eps, by = eps/2)
      xx_ = xx*(support[2]-support[1]) + support[1]
      lines(xx_, dxpolygram(xx_, polygram_object, p = derivative), ...)
    }

  } else {
    plot(xx, yy_hist, bty = "l", ylim = ylim, xlab = xlab, ylab = ylab,
         main = main, type = "l", ...)

  }

  if(bins) {
    abline(v = s_aug, col = "grey", lty = 3)
  }
}

#' @export
lines.polygram = function(polygram_object, derivative = 0, eps = NULL, rug = FALSE, bins = FALSE, ...) {

  support = attr(polygram_object, "support")
  s       = c(attr(polygram_object, "s"))
  s_aug   = c(support[1], attr(polygram_object, "s"), support[2])
  s_aug_  = c(0,(s - support[1])/(support[2] - support[1]),1)
  xx = seq(0, 1, by = 0.001)
  xx = xx*(support[2]-support[1]) + support[1]

  K = attr(polygram_object, "K")
  m = attr(polygram_object, "m")
  p = attr(polygram_object, "p")

  if(!is.integer(derivative) & !is.numeric(derivative)) {
    stop("Option 'derivative' must be an integer greater than -2.")
  } else if (derivative < -1) {
    stop("Option 'derivative' must be an integer greater than -2.")
  }

  yy_hist = dxpolygram(xx, polygram_object, p = derivative)

  if(any(p == -1) | any(m == 0)) {

    eps = 0.001

    for (i in 1:(length(s)+1)) {
      xx_ = seq(s_aug_[i] + eps, s_aug_[i+1], by = 0.0001)
      xx_ = xx*(support[2]-support[1]) + support[1]
      lines(xx_, dxpolygram(xx_, polygram_object, p = derivative), ...)
    }

  } else {
    lines(xx, yy_hist, ...)

  }

  if(bins) {
    abline(v = s_aug, col = "grey", lty = 3)
  }

}

#' @export
points.polygram = function(polygram_object, derivative = 0, eps = NULL, rug = FALSE, bins = FALSE, ...) {
  lines.polygram(polygram_object, derivative = derivative, eps = eps,
                 rug = rug, bins = bins, type = "p", ...)
}
