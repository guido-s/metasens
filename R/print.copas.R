#' Print method for Copas selection model
#' 
#' Print method for objects of class \code{copas}.
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
#' @aliases print.copas
#' 
#' @param x An object of class \code{copas}.
#' @param sign.rsb The significance level for the test of residual selection bias.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.se Minimal number of significant digits for standard
#'   deviations and standard errors, see \code{print.default}.
#' @param ... other arguments to the function will be ignored (this
#'   option included only to conform with R standards)
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
#' data(Fleiss93)
#' 
#' # Perform meta analysis, effect measure is odds ratio (OR)
#' #
#' m1 <- metabin(event.e, n.e, event.c, n.c, data = Fleiss93, sm = "OR")
#' 
#' # Perform Copas analysis
#' #
#' cop1 <- copas(m1)
#' cop1
#' @export print.copas
#' @export
#'
#' @importFrom meta gs
#' @importFrom stats pnorm


print.copas <- function(x, sign.rsb = x$sign.rsb,
                        backtransf = x$backtransf,
                        digits = gs("digits"),
                        digits.se = gs("digits.se"),
                        ...) {
  
  meta:::chkclass(x, "copas")
  ##
  chklevel <- meta:::chklevel
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.copas"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  if (is.null(sign.rsb))
    sign.rsb <- 0.1
  else
    chklevel(sign.rsb)
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  else
    chklogical(backtransf)
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.se, min = 0, single = TRUE)
  
  
  meta:::crtitle(x)
  
  cat("Copas selection model analysis\n\n")
  
  
  res <- cbind(c("range of gamma0: ", "range of gamma1: "),
               format(c(round(x$gamma0.range[1], digits),
                        round(x$gamma1.range[1], digits))),
               ##
               format(c(round(x$gamma0.range[2], digits),
                        round(x$gamma1.range[2], digits))))
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("", "min", "max"))
  ##
  prmatrix(res, quote = FALSE, right = TRUE)
  
  
  cat("\nLargest standard error (SE):", max(round(x$seTE, digits.se)), "\n\n")
  ##
  cat("Range of probability publishing trial with largest SE:\n")
  ##
  res <- matrix(format(round(range(pnorm(x$gamma0 + x$gamma1 / max(x$seTE))),
                             digits)), nrow = 1)
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("min", "max"))
  ##
  prmatrix(res, quote = FALSE, right = TRUE)

  cat("\n\nCalculation of orthogonal line:\n\n")
  ##
  res <- as.matrix(data.frame(x$regr)[ ,c("levels", "nobs",
                                          "adj.r.squareds",
                                          "slopes", "se.slopes")])
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("level", "nobs",
                          "adj.r.square",
                          "slope", "se.slope"))
  prmatrix(res, quote = FALSE, right = TRUE)
  
  cat("\n\n")
  print(summary(x, sign.rsb = sign.rsb),
        header = FALSE,
        digits = digits, backtransf = backtransf,
        ...)
  
  invisible(NULL)
}
