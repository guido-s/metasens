#' Print results of Copas selection model
#' 
#' @description
#' Print method for objects of class \code{copas}.
#' 
#' This function prints the main results of a Copas analysis,
#' performed using the function \code{copas}. It complements the
#' graphical summary of the results, generated using
#' \code{\link{plot.copas}}.
#' 
#' Specifically it prints a table where the:
#' 
#' first column corresponds to the x-axis in plots 3 & 4 from
#' \code{plot.copas};
#' 
#' second column corresponds to the treatment effect displayed in plot
#' 3 from \code{plot.copas};
#' 
#' third and fourth columns give the confidence intervals for this
#' treatment effect,
#' 
#' fifth colum gives the p-value for an overall treatment effect,
#' 
#' sixth column gives the p-value for residual publication bias (the
#' y-axis of plot 4 from \code{plot.copas} (see
#' \code{\link{plot.copas}} under plot 4 for a further explanation of
#' this p-value))
#' 
#' seventh column gives an approximate estimate of the number of
#' studies the model suggests remain unpublished if the probability of
#' publishing the study with the largest SE is as in column 1.
#' 
#' Below this is displayed the results of the Copas analysis (Adjusted
#' estimate) for the smallest degree of selection for which the
#' p-value for evidence of residual selection bias exceeds
#' \code{sign.rsb} (default: 0.1). This is simply extracted from the
#' corresponding row in the table above.
#' 
#' Lastly, the unadjusted random effects estimate and 95\% confidence
#' interval is printed.
#' 
#' @aliases print.copas
#' 
#' @param x An object of class \code{copas}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If \code{backtransf =
#'   TRUE} (default), results are printed as odds ratios rather than
#'   log odds ratio, for example.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.prop Minimal number of significant digits for
#'   proportions, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance \eqn{\tau^2}, see \code{print.default}.
#' @param digits.tau Minimal number of significant digits for
#'   \eqn{\tau}, the square root of the between-study variance
#'   \eqn{\tau^2}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param legend A logical indicating whether a legend should be
#'   printed.
#' @param text.tau2 Text printed to identify between-study variance
#'   \eqn{\tau^2}.
#' @param text.tau Text printed to identify \eqn{\tau}, the square
#'   root of the between-study variance \eqn{\tau^2}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param ... Additional arguments (ignored).
#'
#' @author James Carpenter \email{James.Carpenter@@lshtm.ac.uk}, Guido
#' Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
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
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin, sm = "OR")
#' 
#' # Perform Copas analysis
#' #
#' cop1 <- copas(m1)
#' cop1
#'
#' @method print copas
#' @export
#' @export print.copas


print.copas <- function(x,
                        backtransf = x$backtransf,
                        digits = gs("digits"),
                        digits.pval = max(gs("digits.pval"), 2),
                        digits.prop=gs("digits.prop"),
                        digits.tau2 = gs("digits.tau2"),
                        digits.tau = gs("digits.tau"),
                        scientific.pval=gs("scientific.pval"),
                        big.mark=gs("big.mark"),
                        header = TRUE, legend = TRUE,
                        text.tau2 = gs("text.tau2"),
                        text.tau = gs("text.tau"),
                        ...) {
  
  
  chkclass(x, "copas")
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.copas"
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
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.prop, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.tau, min = 0, length = 1)
  chklogical(scientific.pval)
  ##
  chklogical(header)
  chklogical(legend)
  ##
  chkchar(text.tau2, length = 1)
  chkchar(text.tau, length = 1)
  
  
  sm <- x$sm
  relative <- is.relative.effect(sm)
  ##
  if (!backtransf & relative)
    sm.lab <- paste("log", sm, sep = "")
  else
    sm.lab <- sm
  
  
  if (backtransf & relative) {
    x$TE.slope    <- exp(x$TE.slope)
    x$lower.slope <- exp(x$lower.slope)
    x$upper.slope <- exp(x$upper.slope)
    ##
    x$TE.adjust    <- exp(x$TE.adjust)
    x$lower.adjust <- exp(x$lower.adjust)
    x$upper.adjust <- exp(x$upper.adjust)
    ##
    x$TE.random    <- exp(x$TE.random)
    x$lower.random <- exp(x$lower.random)
    x$upper.random <- exp(x$upper.random)
  }  
  
  
  TE.adj    <- round(x$TE.adjust, digits)
  lowTE.adj <- round(x$lower.adjust, digits)
  uppTE.adj <- round(x$upper.adjust, digits)
  tau2.adj  <- round(x$tau.adjust^2, digits.tau2)
  tau.adj   <- round(x$tau.adjust, digits.tau)
  ##
  TE.random    <- round(x$TE.random, digits)
  lowTE.random <- round(x$lower.random, digits)
  uppTE.random <- round(x$upper.random, digits)
  tau2         <- round(replaceNULL(x$tau, NA)^2, digits.tau2)
  tau          <- round(replaceNULL(x$tau), digits.tau)
  ##  
  TE.slope    <- round(x$TE.slope, digits)
  lowTE.slope <- round(x$lower.slope, digits)
  uppTE.slope <- round(x$upper.slope, digits)
  tau2.slope  <- round(x$tau.slope^2, digits.tau2)
  tau.slope   <- round(x$tau.slope, digits.tau)
  ##
  publprob <- round(x$publprob, digits.prop)
  ## 
  ci.lab <- paste(round(100 * x$level.ma, 1), "%-CI", sep = "")
  
  
  res <- cbind(c(formatN(publprob, digits.prop, ""),
                 "", "Adjusted estimate", "Unadjusted estimate"),
               formatN(c(TE.slope, NA, TE.adj, TE.random),
                       digits, "", big.mark = big.mark),
               formatCI(formatN(c(lowTE.slope, NA, lowTE.adj, lowTE.random),
                                digits, "NA", big.mark = big.mark),
                        formatN(c(uppTE.slope, NA, uppTE.adj, uppTE.random),
                                digits, "NA", big.mark = big.mark)),
               formatPT(c(tau2.slope, NA, tau2.adj, tau2),
                        digits = digits.tau2, lab.NA = "",
                        big.mark = big.mark),
               formatPT(c(tau.slope, NA, tau.adj, tau),
                        digits = digits.tau, lab.NA = "",
                        big.mark = big.mark),
               formatPT(c(x$pval.slope, NA, x$pval.adjust, x$pval.random),
                        digits = digits.pval,
                        scientific = scientific.pval),
               formatPT(c(x$pval.rsb, NA, x$pval.rsb.adj, NA),
                        digits = digits.pval,
                        scientific = scientific.pval),
               formatN(c(round(x$N.unpubl), NA,
                         round(x$N.unpubl.adj), NA),
                       0, "", big.mark = big.mark)
               )
  ##
  res[rmSpace(res) == "--"] <- ""
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("p.publ", sm.lab, ci.lab,
                          text.tau2, text.tau,
                          "p.trt", "p.rsb", "N"))
  
  if (header)
    crtitle(x)
  ##
  cat("Copas selection model analysis\n\n")
  ##
  prmatrix(res, quote = FALSE, right = TRUE, ...)
  ##
  cat("\nSignificance level for test of residual selection bias:",
      ifelse(is.null(x$sign.rsb), 0.1, x$sign.rsb), "\n")
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
