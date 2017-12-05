print.orbbound <- function(x,
                           comb.fixed = x$x$comb.fixed,
                           comb.random = x$x$comb.random,
                           header = TRUE, backtransf = x$backtransf,
                           digits = gs("digits"),
                           digits.zval = gs("digits.zval"),
                           digits.pval = max(gs("digits.pval"), 2),
                           digits.tau2 = gs("digits.tau2"),
                           scientific.pval = gs("scientific.pval"),
                           big.mark = gs("big.mark"),
                           ...) {
  
  
  meta:::chkclass(x, "orbbound")
  ##
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatCI <- meta:::formatCI
  formatN <- meta:::formatN
  formatPT <- meta:::formatPT
  rmSpace <- meta:::rmSpace
  
  
  k <- x$x$k
  k.suspect <- x$k.suspect
  sm <- x$x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.orbbound"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg = "backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  ##
  chknumeric(digits, min = 0, single = TRUE)
  chknumeric(digits.zval, min = 0, single = TRUE)
  chknumeric(digits.pval, min = 1, single = TRUE)
  chknumeric(digits.tau2, min = 0, single = TRUE)
  ##
  chklogical(header)
  chklogical(backtransf)
  chklogical(scientific.pval)
  
  
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
  
  
  if (length(comb.fixed) == 0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random) == 0)
    comb.random <- TRUE
  
  
  ci.lab <- paste(round(100 * x$x$level.comb, 1), "%-CI", sep = "")
  
  
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  maxbias <- x$maxbias
  
  
  if (backtransf) {
    ##
    npft.ma <- 1 / mean(1 / x$x$n)
    ##
    TE.fixed    <- meta:::backtransf(TE.fixed, sm, "mean",
                                     npft.ma, warn = comb.fixed)
    lowTE.fixed <- meta:::backtransf(lowTE.fixed, sm, "lower",
                                     npft.ma, warn = comb.fixed)
    uppTE.fixed <- meta:::backtransf(uppTE.fixed, sm, "upper",
                                     npft.ma, warn = comb.fixed)
    ##
    TE.random <- meta:::backtransf(TE.random, sm, "mean",
                                   npft.ma, warn = comb.random)
    lowTE.random <- meta:::backtransf(lowTE.random, sm, "lower",
                                      npft.ma, warn = comb.random)
    uppTE.random <- meta:::backtransf(uppTE.random, sm, "upper",
                                      npft.ma, warn = comb.random)
    ##
    maxbias <- meta:::backtransf(maxbias, sm, "mean", npft.ma, warn = FALSE)
  }
  ##
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits)
  ##
  maxbias <- round(maxbias, digits)
  
  
  if (header)
    meta:::crtitle(x$x)
  
  
  if (comb.fixed | comb.random) {
    tau2 <- formatPT(x$tau^2,
                     lab = TRUE, labval = "tau^2",
                     digits = digits.tau2,
                     lab.NA = "NA",
                     big.mark = big.mark)
    
    cat("\n        Sensitivity Analysis for Outcome Reporting Bias (ORB)\n\n")
    
    cat(paste("Number of studies combined: k=", x$x$k, "\n", sep = ""))
    cat(paste("Between-study variance: ", tau2, "\n\n", sep = ""))
  }
  
  if (comb.fixed) {
    if (comb.random & x$tau == 0)
      cat("Fixed effect model / Random effects model\n\n")
    else
      cat("Fixed effect model\n\n")
    
    res <- cbind(k.suspect,
                 formatN(maxbias, digits, big.mark = big.mark),
                 formatN(TE.fixed, digits, "NA", big.mark = big.mark),
                 formatCI(formatN(lowTE.fixed, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.fixed, digits, "NA",
                                  big.mark = big.mark)),
                 formatN(round(zTE.fixed, digits = digits.zval),
                         digits.zval, big.mark = big.mark),
                 formatPT(pTE.fixed, digits = digits.pval,
                          scientific = scientific.pval))
    
    zlab <- "z"
    
    names(res) <- c("k.suspect", "maxbias",
                    sm.lab, ci.lab, zlab, "p-value")
    
    dimnames(res) <- list(rep("", length(k.suspect)),
                          c("k.suspect", "maxbias",
                            sm.lab, ci.lab, zlab, "p-value"))
    
    prmatrix(res, quote = FALSE, right = TRUE)
  }
  
  
  if (comb.random & (x$tau != 0 | !comb.fixed)) {
    cat("\nRandom effects model\n\n")
    
    res <- cbind(k.suspect,
                 formatN(maxbias, digits, big.mark = big.mark),
                 formatN(TE.random, digits, "NA", big.mark = big.mark),
                 formatCI(formatN(lowTE.random, digits, "NA",
                                  big.mark = big.mark),
                          formatN(uppTE.random, digits, "NA",
                                  big.mark = big.mark)),
                 formatN(round(zTE.random, digits = digits.zval),
                         digits.zval, big.mark = big.mark),
                 formatPT(pTE.random, digits = digits.pval,
                          scientific = scientific.pval))
    
    zlab <- "z"
    
    names(res) <- c("k.suspect", "maxbias",
                    sm.lab, ci.lab, zlab, "p-value")
    
    dimnames(res) <- list(rep("", length(k.suspect)),
                          c("k.suspect", "maxbias",
                            sm.lab, ci.lab, zlab, "p-value"))
    
    prmatrix(res, quote = FALSE, right = TRUE)
  }
  
  
  if (comb.fixed | comb.random) {
    ## Print information on summary method:
    meta:::catmeth(method = x$x$method,
                   method.tau = if (comb.random) x$x$method.tau else "",
                   sm = sm,
                   k.all = 666)
  }
  
  
  invisible(NULL)
}
