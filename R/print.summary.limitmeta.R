#' Print detailed results for limit meta-analysis
#' 
#' @description
#' Print method for objects of class \code{summary.limitmeta}.
#' 
#' This function prints the main results of a limit meta-analysis
#' (Rücker et al., 2011) as well as the following study information:
#' 
#' \itemize{
#' \item Effect estimate with confidence interval
#' \item Shrunken effect estimates with confidence interval
#' }
#' 
#' @aliases print.summary.limitmeta
#' 
#' @param x An object of class \code{summary.limitmeta}
#' @param sortvar An optional vector used to sort the individual
#'   studies (must be of same length as \code{x$TE}).
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param big.mark A character used as thousands separator.
#' @param truncate An optional vector used to truncate the printout of
#'   results for individual studies (must be a logical vector of same
#'   length as \code{x$TE} or contain numerical values).
#' @param text.truncate A character string printed if study results
#'   were truncated from the printout.
#' @param \dots Additional arguments which are
#'   passed on to \code{print.limitmeta} called internally.
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#'
#' @seealso \code{\link{limitmeta}}, \code{\link{summary.limitmeta}}
#' 
#' @keywords print
#' 
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'   data = Moore1998, sm = "OR", method = "Inverse")
#' 
#' print(summary(limitmeta(m1)), digits = 2)
#'
#' @method print summary.limitmeta
#' @export


print.summary.limitmeta <- function(x,
                                    sortvar,
                                    backtransf = x$backtransf,
                                    digits = gs("digits"),
                                    big.mark = gs("big.mark"),
                                    truncate,
                                    text.truncate = "*** Output truncated ***",
                                    ...) {
  
  
  chkclass(x, "summary.limitmeta")
  ##
  sm <- x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.limitmeta"
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
  
  
  sm.lab <- sm
  ##
  if (backtransf) {
    if (sm == "ZCOR")
      sm.lab <- "COR"
    if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT", "PLN"))
      sm.lab <- "proportion"
  }
  else 
    if (is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep = "")
  
  
  ci.lab <- paste(round(100 * x$level, 1), "%-CI", sep = "")
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  ##
  k.all <- length(TE)
  
  
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  error <-
    try(sortvar <- catch("sortvar", mc, x, sfsp),
        silent = TRUE)
  if (inherits(error, "try-error")) {
    sortvar <- catch("sortvar", mc, x$data, NULL)
    if (isCol(x$data, ".subset"))
      sortvar <- sortvar[x$data$.subset]
  }
  sort <- !is.null(sortvar)
  if (sort && (length(sortvar) != k.all))
    stop("Number of studies in object 'x' and ",
         "argument 'sortvar' have different length.")
  if (!sort)
    sortvar <- 1:k.all
    ##
  ## Catch 'truncate' from meta-analysis object:
  ##
  missing.truncate <- missing(truncate)
  if (!missing.truncate) {
    truncate <- catch("truncate", mc, x, sfsp)
    ##
    if (is.null(truncate))
      truncate <- catch("truncate", mc, x$data, sfsp)
    ##
    if (length(truncate) > k.all)
      stop("Length of argument 'truncate' is too long.",
           call. = FALSE)
    else if (length(truncate) < k.all) {
      if (is.numeric(truncate)) {
        if (any(is.na(truncate)) | max(truncate) > k.all | min(truncate) < 0)
          stop("Numeric values in argument 'truncate' must be between 1 and ",
               k.all, ".",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[truncate] <- TRUE
        truncate <- truncate2
      }
      else if (is.character(truncate)) {
        if (any(!(truncate %in% x$studlab)))
          stop("At least one value of argument 'truncate' does not ",
               "match a study label.",
               call. = FALSE)
        truncate2 <- rep(FALSE, k.all)
        truncate2[x$studlab %in% truncate] <- TRUE
        truncate <- truncate2
      }
      else
        stop("Argument 'truncate' must contain integers or study labels if ",
             "length differs from number of treatment effects.",
             call. = FALSE)
    }
  }
  
  
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
    TE    <- backtransf(TE, sm, "mean", npft.ma, warn = TRUE)
    lowTE <- backtransf(lowTE, sm, "lower", npft.ma, warn = TRUE)
    uppTE <- backtransf(uppTE, sm, "upper", npft.ma, warn = TRUE)
    ##
    TE.limit <- backtransf(TE.limit, sm, "mean", npft.ma, warn = TRUE)
    lowTE.limit <- backtransf(lowTE.limit, sm, "lower", npft.ma, warn = TRUE)
    uppTE.limit <- backtransf(uppTE.limit, sm, "upper", npft.ma, warn = TRUE)
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
  cat("Results for individual studies\n(left: original data; right: shrunken estimates)\n\n")
  ##
  if (!missing.truncate) {
    sortvar <- sortvar[truncate]
    res <- res[truncate, , drop = FALSE]
  }
  ##
  prmatrix(res[order(sortvar), , drop = FALSE],
           quote = FALSE, right = TRUE)
  if (!missing.truncate)
    cat(text.truncate, "\n")
  cat("\n")
  
  x.limit <- x
  class(x.limit) <- class(x)[-1]
  print.limitmeta(x,
                  header = FALSE,
                  digits = digits, backtransf = backtransf,
                  big.mark = big.mark,
                  ...)
  
  
  invisible(NULL)
}
