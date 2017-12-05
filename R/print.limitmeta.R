print.limitmeta <- function(x,
                            sortvar,
                            backtransf = x$backtransf,
                            digits = gs("digits"),
                            big.mark = gs("big.mark"),
                            ...) {
  
  
  meta:::chkclass(x, "limitmeta")
  ##
  bt <- meta:::backtransf
  chklogical <- meta:::chklogical
  chknumeric <- meta:::chknumeric
  formatN <- meta:::formatN
  formatCI <- meta:::formatCI
  ##
  sm <- x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.limitmeta"
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
  
  
  ci.lab <- paste(round(100 * x$level, 1), "%-CI", sep = "")
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  
  
  k.all <- length(TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("Arguments 'x' and 'sortvar' have different length")
  
  
  ci.TE <- ci(TE, seTE, level = x$level)
  lowTE <- ci.TE$lower
  uppTE <- ci.TE$upper
  ##
  ci.limit <- ci(TE.limit, seTE.limit, level = x$level)
  lowTE.limit <- ci.limit$lower
  uppTE.limit <- ci.limit$upper
  
  
  if (backtransf) {
    ##
    npft.ma <- 1 / mean(1 / x$x$n)
    ##
    TE    <- bt(TE, sm, "mean", npft.ma, warn = TRUE)
    lowTE <- bt(lowTE, sm, "lower", npft.ma, warn = TRUE)
    uppTE <- bt(uppTE, sm, "upper", npft.ma, warn = TRUE)
    ##
    TE.limit <- bt(TE.limit, sm, "mean", npft.ma, warn = TRUE)
    lowTE.limit <- bt(lowTE.limit, sm, "lower", npft.ma, warn = TRUE)
    uppTE.limit <- bt(uppTE.limit, sm, "upper", npft.ma, warn = TRUE)
  }
  ##
  TE    <- round(TE, digits)
  lowTE <- round(lowTE, digits)
  uppTE <- round(uppTE, digits)
  ##
  TE.limit    <- round(TE.limit, digits)
  lowTE.limit <- round(lowTE.limit, digits)
  uppTE.limit <- round(uppTE.limit, digits)
  
  
  res <- cbind(" ",
               formatN(TE, digits, "NA", big.mark = big.mark),
               formatCI(formatN(lowTE, digits, "NA",
                                big.mark = big.mark),
                        formatN(uppTE, digits, "NA",
                                big.mark = big.mark)),
               "  ",
               formatN(TE.limit, digits, "NA", big.mark = big.mark),
               formatCI(formatN(lowTE.limit, digits, "NA",
                                big.mark = big.mark),
                        formatN(uppTE.limit, digits, "NA",
                                big.mark = big.mark)))
  ##
  dimnames(res) <-
    list(x$studlab, c("", sm.lab, ci.lab, "", sm.lab, ci.lab))
  ##
  cat("Results for individual studies (left: original data; right: shrunken estimates)\n\n")
  prmatrix(res[order(sortvar), ], quote = FALSE, right = TRUE)
  
  
  cat("\n")
  
  
  print(summary(x),
        header = FALSE,
        digits = digits, backtransf = backtransf, big.mark = big.mark,
        ...)
  
  
  invisible(NULL)
}
