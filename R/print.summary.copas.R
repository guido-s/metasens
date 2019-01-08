#' Print method for summary of Copas selection model
#' 
#' Print method for objects of class \code{summary.copas}.
#' 
#' This function prints a summary of a Copas analysis, performed using
#' the function \code{copas}. It complements the graphical summary of
#' the results, generated using \code{plot.copas}.
#' 
#' Specifically it prints a table where the:
#' 
#' first column corresponds to the x-axis in plots 3 & 4 from
#' \code{plot.copas};
#' 
#' second column corresponds to the treatment effect displayed in plot 3 from
#' \code{plot.copas};
#' 
#' third and fourth columns give the confidence intervals for this treatment
#' effect,
#' 
#' fifth colum gives the p-value for an overall treatment effect,
#' 
#' sixth column gives the p-value for residual publication bias (the y-axis of
#' plot 4 from \code{plot.copas} (see help(plot.copas) under plot 4 for a
#' further explanation of this p-value))
#' 
#' seventh column gives an approximate estimate of the number of studies the
#' model suggests remain unpublished if the probability of publishing the study
#' with the largest SE is as in column 1.
#' 
#' Below this is displayed the results of the Copas analysis for the smallest
#' degree of selection for which the p-value for evidence of residual selection
#' bias exceeds \code{sign.rsb} (default: 0.1). This is simply extracted from
#' the corresponding row in the table above.
#' 
#' Lastly, the usual random effects estimate (based on the DerSimonian-Laird
#' method) and 95\% confidence interval is printed.
#' 
#' @aliases print.summary.copas
#' 
#' @param x An object of class \code{summary.copas}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
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
#' m1 <- metabin(event.e, n.e, event.c, n.c, data=Fleiss93, sm="OR")
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


print.summary.copas <- function(x,
                                backtransf = x$backtransf,
                                digits = gs("digits"),
                                digits.pval = max(gs("digits.pval"), 2),
                                digits.prop=gs("digits.prop"),
                                scientific.pval=gs("scientific.pval"),
                                big.mark=gs("big.mark"),
                                header = TRUE,
                                ...) {
  
  
  meta:::chkclass(x, "summary.copas")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.copas"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
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
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.prop, min = 0, single = TRUE)
  chklogical(scientific.pval)
  
  
  sm <- x$sm
  relative <- meta:::is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  else
    sm.lab <- sm
  
  
  if (backtransf & relative) {
    x$slope$TE    <- exp(x$slope$TE)
    x$slope$lower <- exp(x$slope$lower)
    x$slope$upper <- exp(x$slope$upper)
    ##
    x$adjust$TE    <- exp(x$adjust$TE)
    x$adjust$lower <- exp(x$adjust$lower)
    x$adjust$upper <- exp(x$adjust$upper)
    ##
    x$random$TE    <- exp(x$random$TE)
    x$random$lower <- exp(x$random$lower)
    x$random$upper <- exp(x$random$upper)
  }  
  
  
  TE.adj    <- round(x$adjust$TE, digits)
  lowTE.adj <- round(x$adjust$lower, digits)
  uppTE.adj <- round(x$adjust$upper, digits)
  ##
  TE.random    <- round(x$random$TE, digits)
  lowTE.random <- round(x$random$lower, digits)
  uppTE.random <- round(x$random$upper, digits)
  ##  
  TE.slope    <- round(x$slope$TE, digits)
  lowTE.slope <- round(x$slope$lower, digits)
  uppTE.slope <- round(x$slope$upper, digits)
  ##
  publprob <- round(x$publprob, digits.prop)
  pval.TE  <- round(x$slope$p, digits.pval)
  pval.rsb <- round(x$pval.rsb, digits.pval)
  
  
  res <- cbind(c(formatN(publprob, digits.prop, ""),
                 "", "Copas model (adj)","Random effects model"),
               formatN(c(TE.slope, NA, TE.adj, TE.random),
                       digits, "", big.mark = big.mark),
               formatCI(formatN(c(lowTE.slope, NA, lowTE.adj, lowTE.random),
                                digits, "NA", big.mark = big.mark),
                        formatN(c(uppTE.slope, NA, uppTE.adj, uppTE.random),
                                digits, "NA", big.mark = big.mark)),
               formatPT(c(pval.TE, NA, x$adjust$p, x$random$p),
                        digits = digits.pval,
                        scientific = scientific.pval),
               formatPT(c(pval.rsb, NA, x$pval.rsb.adj, NA),
                        digits = digits.pval,
                        scientific = scientific.pval),
               formatN(c(round(x$N.unpubl), NA,
                         round(x$N.unpubl.adj), NA),
                       0, "", big.mark = big.mark)
               )
  ##
  res[meta:::rmSpace(res) == "--"] <- ""
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("publprob", sm.lab, x$ci.lab,
                          "pval.treat", "pval.rsb", "N.unpubl"))
  
  if (header)
    meta:::crtitle(x)
  
  cat("Summary of Copas selection model analysis:\n\n")
  ##
  prmatrix(res, quote = FALSE, right = TRUE, ...)
  ##
  cat("\n",
      "Significance level for test of residual selection bias:",
      ifelse(is.null(x$sign.rsb), 0.1, x$sign.rsb), "\n")
  ##
  cat("\n Legend:\n")
  cat(" publprob   - Probability of publishing study with largest standard error\n")
  cat(" pval.treat - P-value for hypothesis of overall treatment effect\n")
  cat(" pval.rsb   - P-value for hypothesis that no selection remains unexplained\n")
  cat(" N.unpubl   - Approximate number of unpublished studies suggested by model\n")
  
  invisible(NULL)
}
