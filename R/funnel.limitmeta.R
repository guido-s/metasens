#' Funnel plot for limit meta-analysis
#' 
#' Draws a funnel plot in the active graphics window.
#' 
#' A funnel plot is drawn in the active graphics window. In addition
#' this function adds the adjusted effect estimate as well as a
#' nonlinear regression line (also called adjusted regression line) if
#' argument \code{line} is \code{TRUE}. The adjusted regression line
#' is representing the dependence of the treatment effect estimate on
#' the standard error across studies. The adjusted regression line is
#' only plotted in addition to the adjusted treatment effect if
#' argument \code{method.adjust="beta0"} (default) has been used in
#' the \code{\link{limitmeta}} function.
#' 
#' If argument \code{shrunken} is \code{TRUE} the shrunken effect
#' estimates are also plotted. Lines are connecting original and
#' shrunken effect estimates.
#' 
#' Internally, R function \code{\link{funnel.meta}} is called to
#' create a funnel plot. For more information see help page of the
#' \code{\link[meta]{funnel}} function.
#' 
#' @param x An object of class \code{limitmeta}.
#' @param pch The plotting symbol used for individual studies.
#' @param cex The magnification to be used for plotting symbol.
#' @param col A vector with colour of plotting symbols.
#' @param bg A vector with background colour of plotting symbols (only
#'   used if \code{pch} in \code{21:25}).
#' @param lwd The line width for confidence intervals (see
#'   \code{\link[meta]{funnel}}).
#' @param show.ci.adjust A logical indicating whether to show the
#'   confidence interval of the adjusted estimate.
#' @param pch.adjust The plotting symbol used for the adjusted effect
#'   estimate.
#' @param cex.adjust The magnification to be used for the plotting
#'   symbol of the adjusted effect estimate.
#' @param col.adjust Colour of plotting symbol for adjusted effect
#'   estimate.
#' @param bg.adjust Background colour of plotting symbol for adjusted
#'   effect estimate.
#' @param line A logical indicating whether adjusted regression line
#'   should be plotted.
#' @param xmin.line Minimal value for the adjusted regression line (on
#'   x-axis).
#' @param xmax.line Maximum value for the adjusted regression line (on
#'   x-axis).
#' @param lty.line Line type of the adjusted regression line.
#' @param col.line Color of the adjusted regression line.
#' @param lwd.line The line width of the adjusted regression line.
#' @param shrunken A logical indicating whether shrunken treatment
#'   estimates should be plotted.
#' @param pch.shrunken The plotting symbol used for shrunken effect
#'   estimates.
#' @param cex.shrunken The magnification to be used for the plotting
#'   symbol of the shrunken effect estimates.
#' @param col.shrunken Colour of plotting symbol for shrunken effect
#'   estimates.
#' @param bg.shrunken Background colour of plotting symbol for
#'   shrunken effect estimates.
#' @param lty.connect Line type for line connecting original and
#'   shrunken treatment estimates.
#' @param lwd.connect The line width of the connecting lines.
#' @param col.connect Color of the connecting lines.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param \dots Additional arguments for \code{\link[meta]{funnel}}
#'   function.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}, Gerta
#'   Rücker \email{ruecker@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{limitmeta}}, \code{\link[meta]{funnel}}
#'
#' @keywords hplot
#' 
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'               data = Moore1998, sm = "OR", method = "Inverse")
#' 
#' print(summary(limitmeta(m1)), digits = 2)
#' funnel(limitmeta(m1))
#' 
#' # Print results on log scale
#' #
#' print(summary(limitmeta(m1)), digits = 2, backtransf = FALSE)
#' funnel(limitmeta(m1), backtransf = FALSE)
#'
#' @method funnel limitmeta
#' @export
#' @export funnel.limitmeta
#'
#' @importFrom meta funnel
#' @importFrom graphics curve points polygon segments


funnel.limitmeta <- function(x,
                             ##
                             pch = 21,
                             cex = 1,
                             col = "black",
                             bg = "darkgray",
                             ##
                             lwd = 1,
                             ##
                             show.ci.adjust = FALSE,
                             pch.adjust = 18,
                             cex.adjust = 1.5,
                             col.adjust = "gray",
                             bg.adjust = "gray",
                             ##
                             line = TRUE,
                             xmin.line,
                             xmax.line,
                             lty.line = 1,
                             lwd.line = lwd,
                             col.line = "gray",
                             ##
                             shrunken = FALSE,
                             pch.shrunken = 22,
                             cex.shrunken = 1,
                             col.shrunken = "black",
                             bg.shrunken = "white",
                             ##
                             lty.connect = 1,
                             lwd.connect = 0.8,
                             col.connect = "black",
                             ##
                             backtransf = x$backtransf,
                             ...) {
  
  
  meta:::chkclass(x, "limitmeta")
  ##
  meta:::chklogical(show.ci.adjust)
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  ##
  minTE <- min(TE, na.rm = TRUE)
  maxTE <- max(TE, na.rm = TRUE)
  x.incr <- (maxTE - minTE) / 1000
  ##
  TE.adjust <- x$TE.adjust
  lower.adjust <- x$lower.adjust
  upper.adjust <- x$upper.adjust
  ##
  tau <- x$tau
  alpha.r <- x$alpha.r
  beta.r <- x$beta.r
  ##
  sm <- x$sm
  
  
  if (alpha.r < 0) {
    if (missing(xmin.line))
      xmin.line <- minTE
    if (missing(xmax.line))
      xmax.line <- TE.adjust - x.incr
  }
  if (alpha.r > 0) {
    if (missing(xmin.line))
      xmin.line <- TE.adjust + x.incr
    if (missing(xmax.line))
    xmax.line <- maxTE
  }
  
  
  if (backtransf & meta:::is.relative.effect(sm)) {
    TE <- exp(TE)
    TE.limit <- exp(TE.limit)
    TE.adjust <- exp(TE.adjust)
    lower.adjust <- exp(lower.adjust)
    upper.adjust <- exp(upper.adjust)
  }
  
  
  ##
  ## Generate funnel plot
  ##
  f1 <- funnel(x$x, pch = pch, cex = cex, col = col, bg = bg, lwd = lwd,
               backtransf = backtransf, ...)
  
  
  ##
  ## Add line for adjustment method beta0
  ##
  if (line) {
    if (x$method.adjust == "beta0") {
      if (backtransf & meta:::is.relative.effect(sm)) {
        curve(sqrt((log(x) - beta.r)^2 / alpha.r^2 - tau^2),
              from = exp(xmin.line), to = exp(xmax.line),
              lty = lty.line, col = col.line, lwd = lwd.line, add = TRUE)
      }
      else {
        curve(sqrt((x - beta.r)^2 / alpha.r^2 - tau^2),
              from = xmin.line, to = xmax.line,
              lty = lty.line, col = col.line, lwd = lwd.line, add = TRUE)
      }
    }
  }
  
  
  ##
  ## Add adjusted treatment effect
  ##
  if (!show.ci.adjust)
    points(TE.adjust, 0,
           pch = pch.adjust, cex = cex.adjust,
           col = col.adjust, bg = bg.adjust)
  else {
    if (!is.null(f1$ylim))
      y.incr <- (f1$ylim[2] - f1$ylim[1]) / 50
    else
      y.incr <- 0.01
    ##
    polygon(c(lower.adjust, TE.adjust, upper.adjust, TE.adjust),
            0 + c(0, -y.incr, 0, y.incr),
            col = col.adjust, border = bg.adjust)
    }
  
  
  ##
  ## Add lines
  ##
  if (shrunken)
    segments(TE, seTE, TE.limit, seTE.limit,
             lty = lty.connect, lwd = lwd.connect, col = col.connect)
  
  
  ##
  ## Plot studies again
  ##
  points(TE, seTE, pch = pch, cex = cex, col = col, bg = bg)
  
  
  ##
  ## Add shrunken estimates
  ##
  if (shrunken)
    points(TE.limit, seTE.limit,
           pch = pch.shrunken, cex = cex.shrunken, col = col.shrunken, bg = bg.shrunken)
  
  
  invisible(NULL)
}
