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
