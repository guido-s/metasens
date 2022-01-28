#' Print results for limit meta-analysis
#' 
#' Print method for objects of class \code{limitmeta}.
#' 
#' This function prints the main results of a limit meta-analysis
#' (Rücker et al., 2011).
#' 
#' @aliases print.limitmeta
#' 
#' @param x An object of class \code{limitmeta}.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param header A logical indicating whether information on title of
#'   meta-analysis, comparison and outcome should be printed at the
#'   beginning of the printout.
#' @param pscale A numeric giving scaling factor for printing of
#'   single event probabilities, i.e. if argument \code{sm} is equal
#'   to \code{"PLOGIT"}, \code{"PLN"}, \code{"PRAW"}, \code{"PAS"}, or
#'   \code{"PFT"}.
#' @param irscale A numeric defining a scaling factor for printing of
#'   rates, i.e. if argument \code{sm} is equal to \code{"IR"},
#'   \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}.
#' @param irunit A character specifying the time unit used to
#'   calculate rates, e.g. person-years.
#' @param digits.stat Minimal number of significant digits for z- or
#'   t-value, see \code{print.default}.
#' @param digits.pval Minimal number of significant digits for p-value
#'   of overall treatment effect, see \code{print.default}.
#' @param digits.Q Minimal number of significant digits for
#'   heterogeneity statistic Q, see \code{print.default}.
#' @param digits.tau2 Minimal number of significant digits for
#'   between-study variance, see \code{print.default}.
#' @param digits.I2 Minimal number of significant digits for I-squared
#'   and Rb statistic, see \code{print.default}.
#' @param scientific.pval A logical specifying whether p-values should
#'   be printed in scientific notation, e.g., 1.2345e-01 instead of
#'   0.12345.
#' @param big.mark A character used as thousands separator.
#' @param print.Rb A logical specifying whether heterogeneity
#'   statistic Rb should be printed.
#' @param warn.backtransf A logical indicating whether a warning
#'   should be printed if backtransformed proportions and rates are
#'   below 0 and backtransformed proportions are above 1.
#' @param \dots Additional arguments (ignored).
#'
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#'
#' @seealso \code{\link{limitmeta}}, \code{\link{summary.limitmeta}},
#'   \code{\link{print.summary.limitmeta}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'   data = Moore1998, sm = "OR", method = "Inverse")
#' 
#' print(limitmeta(m1), digits = 2)
#'
#' @method print limitmeta
#' @export


print.limitmeta <- function(x,
                            backtransf = x$backtransf,
                            digits = gs("digits"),
                            header = TRUE,
                            pscale = x$x$pscale,
                            irscale = x$x$irscale,
                            irunit = x$x$irunit,
                            digits.stat = gs("digits.stat"),
                            digits.pval = gs("digits.pval"),
                            digits.Q = gs("digits.Q"),
                            digits.tau2 = gs("digits.tau2"),
                            digits.I2 = gs("digits.I2"),
                            scientific.pval = gs("scientific.pval"),
                            big.mark = gs("big.mark"),
                            print.Rb = gs("print.Rb"),
                            warn.backtransf = FALSE,
                            ...) {
  
  
  ##
  ##
  ## (1) Check for summary.limitmeta object
  ##
  ##
  chkclass(x, "limitmeta")
  
  
  ##
  ##
  ## (2) Check and set other arguments
  ##
  ##
  chknumeric(digits, min = 0, length = 1)
  chknumeric(digits.tau2, min = 0, length = 1)
  chknumeric(digits.stat, min = 0, length = 1)
  chknumeric(digits.pval, min = 1, length = 1)
  chknumeric(digits.Q, min = 0, length = 1)
  chknumeric(digits.I2, min = 0, length = 1)
  chklogical(backtransf)
  chklogical(scientific.pval)
  chklogical(print.Rb)
  chklogical(warn.backtransf)
  ##
  is.prop <- is.prop(x$sm)
  is.rate <- is.rate(x$sm)
  ##
  if (!is.prop)
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, length = 1)
  else
    pscale <- 1
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate)
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, length = 1)
  else
    irscale <- 1
  if (!backtransf & irscale != 1) {
    warning("Argument 'irscale' set to 1 as argument 'backtransf' is FALSE.")
    irscale <- 1
  }
  ##
  ## Additional arguments / checks for metacont objects
  ##
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.limitmeta"
  ##
  warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  ci.lab <- paste(round(100 * x$level.ma, 1), "%-CI", sep = "")
  ##  
  k <- x$k
  sm <- x$sm
  ##
  if (is.null(x$df.Q))
    df.Q <- k - 1
  else
    df.Q <- x$df.Q
  ##
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    else if (is.prop) {
      if (pscale == 1)
        sm.lab <- "proportion"
      else
        sm.lab <- "events"
    }
    else if (is.rate) {
      if (irscale == 1)
        sm.lab <- "rate"
      else
        sm.lab <- "events"
    }
  }
  else
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  
  
  ##
  ##
  ## (4) Set and backtransform results of meta-analysis
  ##
  ##
  TEs <- c(x$TE.adjust, x$TE.random)
  lower <- c(x$lower.adjust, x$lower.random)
  upper <- c(x$upper.adjust, x$upper.random)
  ##
  if (backtransf) {
    if (inherits(x$x, "metarate"))
      harmonic.mean <- 1 / mean(1 / x$x$time)
    else
      harmonic.mean <- 1 / mean(1 / x$x$n)
    ##
    TEs   <- backtransf(TEs, sm, "mean",
                               harmonic.mean, warn = warn.backtransf)
    lower <- backtransf(lower, sm, "lower",
                               harmonic.mean, warn = warn.backtransf)
    upper <- backtransf(upper, sm, "upper",
                               harmonic.mean, warn = warn.backtransf)
  }
  ##
  ## Apply argument 'pscale' to proportions and 'irscale' to rates
  ##
  if (is.prop | is.rate) {
    if (is.prop)
      scale <- pscale
    else if (is.rate)
      scale <- irscale
    ##
    TEs    <- scale * TEs
    lower <- scale * lower
    upper <- scale * upper
  }
  ##
  ## Round and round ...
  ##
  TEs   <- round(TEs, digits)
  lower <- round(lower, digits)
  upper <- round(upper, digits)
  ##
  pvals <- c(x$pval.adjust, x$pval.random)
  statistics <- round(c(x$statistic.adjust, x$statistic.random), digits.stat)
  ##
  I2 <- round(100 * x$x$I2, digits.I2)
  lowI2 <- round(100 * x$x$lower.I2, digits.I2)
  uppI2 <- round(100 * x$x$upper.I2, digits.I2)
  ##
  Rb <- round(100 * x$x$Rb, digits.I2)
  lowRb <- round(100 * x$x$lower.Rb, digits.I2)
  uppRb <- round(100 * x$x$upper.Rb, digits.I2)
  ##
  G2 <- round(100 * x$G.squared, digits.I2)
  
  
  ##
  ##
  ## (5) Print result for meta-analysis
  ##
  ##
  if (header)
    crtitle(x)
  ##
  res <- cbind(c("Adjusted estimate",
                 "Unadjusted estimate"),
               formatN(TEs, digits, "NA", big.mark = big.mark),
               formatCI(formatN(lower, digits, "NA", big.mark = big.mark),
                        formatN(upper, digits, "NA", big.mark = big.mark)),
               formatN(statistics, digits.stat, big.mark = big.mark),
               formatPT(pvals, digits = digits.pval,
                        scientific = scientific.pval)
               )
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("Random effects model",
                          sm.lab, ci.lab, "z", "pval"))
  ##
  cat("Result of limit meta-analysis:\n\n")
  ##
  prmatrix(res, quote = FALSE, right = TRUE)
  ##  
  if (!is.na(x$tau)) {
    ##
    cat(paste("\nQuantifying heterogeneity:\n",
              ##
              formatPT(x$tau^2,
                       lab = TRUE, labval = gs("text.tau2"),
                       digits = digits.tau2,
                       big.mark = big.mark,
                       lab.NA = "NA"),
              ##
              paste("; ",
                    "I^2 = ",
                    if (is.nan(I2)) "NA" else paste(formatN(I2, digits.I2),
                                                    "%", sep = ""),
                    ifelse(k > 2 & !(is.na(lowI2) | is.na(uppI2)),
                           paste(" ",
                                 formatCI(paste(formatN(lowI2, digits.I2),
                                                "%", sep = ""),
                                          paste(formatN(uppI2, digits.I2),
                                                "%", sep = "")),
                                 sep = ""),
                           ""),
                    "; ",
                    "G^2 = ",
                    if (is.nan(G2)) "NA" else paste(formatN(G2, digits.I2),
                                                    "%", sep = ""),
                    sep = ""),
              if (print.Rb)
                paste(";\n ",
                      "Rb = ",
                      if (is.nan(Rb)) "NA" else paste(formatN(Rb, digits.I2),
                                                      "%", sep = ""),
                      ifelse(k > 2 & !(is.na(lowRb) | is.na(uppRb)),
                             paste(" ",
                                   formatCI(paste(formatN(lowRb, digits.I2),
                                                  "%", sep = ""),
                                            paste(formatN(uppRb, digits.I2),
                                                  "%", sep = "")),
                                   sep = ""),
                             ""),
                      sep = ""),
              "\n", sep = ""))
  }
  
  Qd <- function(Q, df)
    data.frame(Q = round(Q, 2), df = df,
               p = 1 - pchisq(Q, df = df))
  ##
  Qds <- rbind(Qd(x$Q, k - 1),
               Qd(x$Q.small, 1),
               Qd(x$Q.resid, k - 2))
  Qds$Q  <- format(Qds$Q)
  Qds$df <- format(Qds$df)
  Qds$p  <- formatPT(Qds$p, digits = digits.pval, scientific = scientific.pval)
  ##
  Qd1 <- Qds[1, ]
  Qd2 <- Qds[2, ]
  Qd3 <- Qds[3, ]
  ##
  dimnames(Qd1) <- list("", c("Q", "d.f.", "p-value"))
  cat("\nTest of heterogeneity:\n")
  prmatrix(Qd1, quote = FALSE, right = TRUE, ...)
  ##
  dimnames(Qd2) <- list("", c("Q-Q'", "d.f.", "p-value"))
  cat("\nTest of small-study effects:\n")
  prmatrix(Qd2, quote = FALSE, right = TRUE, ...)
  ##
  dimnames(Qd3) <- list("", c("Q'", "d.f.", "p-value"))
  cat("\nTest of residual heterogeneity beyond small-study effects:\n")
  prmatrix(Qd3, quote = FALSE, right = TRUE, ...)
  
  
  imeth <- charmatch(tolower(x$method.adjust),
                     c("beta0", "betalim", "mulim"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("Argument 'method.adjust' should be \"beta0\", \"betalim\", or \"mulim\"")
  ##
  method.adjust.detail <- c("- expectation (beta0)",
                            "- including bias parameter (beta-lim)",
                            "- excluding bias parameter (mu-lim)")[imeth]
  ##
  cat("\nDetails on adjustment method:\n",
      method.adjust.detail,
      "\n", sep = "")
  
  
  invisible(NULL)
}
