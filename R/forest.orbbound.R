#' Forest plot for \code{orbbound} object (bound for outcome reporting
#' bias)
#' 
#' Draws a forest plot in the active graphics window (using grid
#' graphics system).
#' 
#' A forest plot, also called confidence interval plot, is drawn in
#' the active graphics window.
#' 
#' For relative effect measures, e.g., 'RR', 'OR', and 'HR', the
#' column labeled "Maximum bias" contains the relative bias, e.g. a
#' value of 1.10 means a maximum overestimation by 10 percent. If
#' \code{backtransf=FALSE} for these summary measures, maximum bias is
#' instead printed as absolute bias.
#'
#' Internally, R function \code{\link{forest.meta}} is called to
#' create a forest plot. For more information see help page of the
#' \code{\link{forest.meta}} function.
#'
#' @param x An object of class \code{orbbound}.
#' @param comb.fixed A logical indicating whether sensitivity analysis
#'   for fixed effect model should be plotted.
#' @param comb.random A logical indicating whether sensitivity
#'   analysis for random effects model should be plotted.
#' @param text.fixed A character string used in the plot to label
#'   subgroup with results for fixed effect model.
#' @param text.random A character string used in the plot to label
#'   subgroup with results for random effects model.
#' @param smlab A label printed at top of figure. If only results for
#'   either fixed effect or random effects model is plotted, text
#'   indicates which model was used.
#' @param leftcols A character vector specifying (additional) columns
#'   to be plotted on the left side of the forest plot or a logical
#'   value (see \code{\link{forest.meta}} help page for details).
#' @param leftlabs A character vector specifying labels for
#'   (additional) columns on left side of the forest plot (see
#'   \code{\link{forest.meta}} help page for details).
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param \dots Additional arguments for \code{\link{forest.meta}}
#'   function.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{orbbound}}, \code{\link{print.orbbound}}
#'
#' @keywords hplot
#'
#' @examples
#' data(Fleiss93, package = "meta")
#' 
#' m1 <- metabin(event.e, n.e, event.c, n.c, data = Fleiss93, sm = "OR")
#' 
#' orb1 <- orbbound(m1, k.suspect = 1:5)
#' print(orb1, digits = 2)
#' forest(orb1, xlim = c(0.7, 1.5))
#' \dontrun{forest(orb1, backtransf = FALSE)}
#'
#' @method forest orbbound
#' @export
#' @export forest.orbbound
#'
#' @importFrom meta forest metagen


forest.orbbound <- function(x,
                            comb.fixed = x$x$comb.fixed,
                            comb.random = x$x$comb.random,
                            text.fixed = "FE model",
                            text.random = "RE model",
                            smlab = NULL,
                            leftcols = c("studlab", "maxbias"),
                            leftlabs = c("Missing\nstudies", "Maximum\nbias"),
                            backtransf = x$backtransf,
                            digits = max(3, .Options$digits - 3),
                            ...) {
  
  meta:::chkclass(x, "orbbound")
  
  k <- x$x$k
  sm <- x$x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.limitmeta"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  
  
  if (length(comb.fixed) == 0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random) == 0)
    comb.random <- TRUE
  
  
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT", "PLN"))
      sm.lab <- "proportion"
  }
  else 
    if (meta:::is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  
  
  ci.lab <- paste(round(100 * x$x$level.comb, 1), "%-CI", sep = "")
  
  
  TE.fixed   <- x$fixed$TE
  seTE.fixed <- rep(x$fixed$seTE, length(TE.fixed))
  ##
  TE.random   <- x$random$TE
  seTE.random <- rep(x$random$seTE, length(TE.random))
  
  
  if (comb.fixed & comb.random) {
    TE <- c(TE.fixed, TE.random)
    seTE <- c(seTE.fixed, seTE.random)
    FEvsRE <- c(rep(text.fixed, length(TE.fixed)),
                rep(text.random, length(TE.random)))
  }
  if (comb.fixed & !comb.random) {
    TE <- TE.fixed
    seTE <- seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed effect model"
  }
  if (!comb.fixed & comb.random) {
    TE <- TE.random
    seTE <- seTE.random
    if (is.null(smlab))
      smlab <- "Random effects model"
  }
  if (!comb.fixed & !comb.random) {
    warning("No forest plot generated as both arguments 'comb.fixed' and 'comb.random' are FALSE")
    return(invisible(NULL))
  }
  
  
  if (comb.fixed & comb.random)
    m1 <- metagen(TE, seTE, sm = sm.lab,
                  byvar = FEvsRE, print.byvar = FALSE,
                  warn = FALSE)
  else
    m1 <- metagen(TE, seTE, sm = sm.lab, warn = FALSE)
  ##
  if (comb.fixed & comb.random) {
    m1$studlab <- c(x$k.suspect, x$k.suspect)
    m1$maxbias <- c(x$maxbias, x$maxbias)
    m1$npft.ma <- c(1 / mean(1 / x$x$n), 1 / mean(1 / x$x$n))
  }
  else {
    m1$studlab <- x$k.suspect
    m1$maxbias <- x$maxbias
    m1$npft.ma <- 1 / mean(1 / x$x$n)
  }
  
  
  if (backtransf)
    m1$maxbias <- meta:::backtransf(m1$maxbias, sm, "mean",
                                    m1$npft.ma, warn = FALSE)
  ##
  m1$maxbias <- format(round(m1$maxbias, digits))
  
  forest(m1,
         comb.fixed = FALSE, comb.random = FALSE,
         hetstat = FALSE,
         leftcols = leftcols,
         leftlabs = leftlabs,
         smlab = smlab,
         just.studlab = "center",
         weight.study = "same",
         ...)
  
  
  invisible(NULL)
}
