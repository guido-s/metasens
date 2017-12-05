print.summary.limitmeta <- function(x,
                                    backtransf = x$backtransf,
                                    digits = gs("digits"),
                                    header = TRUE,
                                    pscale = x$x$pscale,
                                    irscale = x$x$irscale,
                                    irunit = x$x$irunit,
                                    digits.zval = gs("digits.zval"),
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
  meta:::chkclass(x, "summary.limitmeta")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  formatCI <- meta:::formatCI
  
  
  ##
  ##
  ## (2) Check and set other arguments
  ##
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.Q, min = 0, single = TRUE)
  chknumeric(digits.I2, min = 0, single = TRUE)
  chklogical(backtransf)
  chklogical(scientific.pval)
  chklogical(print.Rb)
  chklogical(warn.backtransf)
  ##
  is.prop <- meta:::is.prop(x$sm)
  is.rate <- meta:::is.rate(x$sm)
  ##
  if (!is.prop)
    pscale <- 1
  if (!is.null(pscale))
    chknumeric(pscale, single = TRUE)
  else
    pscale <- 1
  if (!backtransf & pscale != 1) {
    warning("Argument 'pscale' set to 1 as argument 'backtransf' is FALSE.")
    pscale <- 1
  }
  if (!is.rate)
    irscale <- 1
  if (!is.null(irscale))
    chknumeric(irscale, single = TRUE)
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
  fun <- "print.summary.limitmeta"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  
  
  ##
  ##
  ## (3) Some additional settings
  ##
  ##
  ci.lab <- paste(round(100 * x$level.comb, 1), "%-CI", sep = "")
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
    if (meta:::is.relative.effect(sm))
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
    TEs   <- meta:::backtransf(TEs, sm, "mean",
                               harmonic.mean, warn = warn.backtransf)
    lower <- meta:::backtransf(lower, sm, "lower",
                               harmonic.mean, warn = warn.backtransf)
    upper <- meta:::backtransf(upper, sm, "upper",
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
  zvals <- round(c(x$zval.adjust, x$zval.random), digits.zval)
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
    meta:::crtitle(x)
  ##
  res <- cbind(c("Adjusted estimate",
                 "Unadjusted estimate"),
               formatN(TEs, digits, "NA", big.mark = big.mark),
               formatCI(formatN(lower, digits, "NA", big.mark = big.mark),
                        formatN(upper, digits, "NA", big.mark = big.mark)),
               formatN(zvals, digits.zval, big.mark = big.mark),
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
