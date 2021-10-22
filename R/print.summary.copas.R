#' Print detailed results of Copas selection model
#' 
#' @description
#' Print method for objects of class \code{summary.copas}.
#' 
#' 
#' This function prints the following information:
#' 
#' Range of gamma0 values used (see \code{help(copas)});
#' 
#' Range of gamma1 values used (see \code{help(copas)});
#' 
#' Largest SE of all studies in meta-analysis;
#' 
#' Range of probability publishing trial with largest SE;
#' 
#' The next table gives details relating to the summary of the contour plot.
#' Specifically, it gives details from fitting a straight line to each
#' treatment-contour in the contour plot. Column 1 (headed level) shows the
#' treatment-contours; column 2 (nobs) shows the number of observations used by
#' the contour plot command within the \code{copas} function to plot this
#' contour line; column 3 (adj.r.square) shows the adjusted r-square from
#' fitting a straight line to this contour; columns 4 & 5 show the slope and
#' its standard error from fitting a straight line to this contour.
#'
#' Next, the printout of \code{summary.copas} is shown.
#' 
#' @aliases print.summary.copas
#' 
#' @param x An object of class \code{summary.copas}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE} (default), results are printed as odds ratios rather than
#'   log odds ratio, for example.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   deviations and standard errors, see \code{print.default}.
#' @param ... Additional arguments (ignored).
#'
#' @author James Carpenter \email{James.Carpenter@@lshtm.ac.uk}, Guido
#' Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{copas}}, \code{\link{plot.copas}},
#'   \code{\link{summary.copas}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Fleiss1993bin, package = "meta")
#' 
#' # Perform meta analysis, effect measure is odds ratio (OR)
#' #
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data=Fleiss1993bin, sm="OR")
#' 
#' # Print summary of Copas analysis
#' #
#' summary(copas(m1), level=0.95)
#'
#' @method print summary.copas
#' @export
#' @export print.summary.copas
#'
#' @importFrom meta gs


print.summary.copas <- function(x, backtransf = x$backtransf,
                                legend = TRUE,
                                digits = gs("digits"),
                                digits.se = gs("digits.se"),
                                ...) {
  
  chkclass(x, "summary.copas")
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summarycopas"
  ##
  warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  else
    chklogical(backtransf)
  ##
  chklogical(legend)
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.se, min = 0, length = 1)
  
  
  crtitle(x)
  ##  
  print.copas(x$x,
              header = FALSE, legend = FALSE,
              digits = digits, backtransf = backtransf,
              ...)
  ##
  res <- cbind(c("range of gamma0: ", "range of gamma1: "),
               format(c(round(x$x$gamma0.range[1], digits),
                        round(x$x$gamma1.range[1], digits))),
               ##
               format(c(round(x$x$gamma0.range[2], digits),
                        round(x$x$gamma1.range[2], digits))))
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("", "min", "max"))
  cat("\n")
  prmatrix(res, quote = FALSE, right = TRUE)
  ##  
  cat("\nLargest standard error (SE):",
      max(round(x$x$seTE, digits.se)), "\n\n")
  ##
  cat("Range of probability publishing trial with largest SE:\n")
  ##
  res <-
    matrix(format(round(range(pnorm(x$x$gamma0 + x$x$gamma1 / max(x$x$seTE))),
                        digits)), nrow = 1)
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("min", "max"))
  ##
  prmatrix(res, quote = FALSE, right = TRUE)
  ##
  cat("\nCalculation of orthogonal line:\n\n")
  ##
  res <-
    as.matrix(data.frame(x$x$regr)[ , c("levels", "nobs",
                                        "adj.r.squareds",
                                        "slopes", "se.slopes")])
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("level", "nobs",
                          "adj.r.square",
                          "slope", "se.slope"))
  prmatrix(res, quote = FALSE, right = TRUE)
  ##
  if (legend) {
    cat("\n Legend:\n")
    cat(" p.publ - Probability of publishing study with largest SE\n")
    cat(" p.trt  - P-value for test of overall treatment effect\n")
    cat(" p.rsb  - P-value for test of residual selection bias\n")
    cat(" N      - Estimated number of unpublished studies\n")
  }
  
  invisible(NULL)
}
