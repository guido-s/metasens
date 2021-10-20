#' Imputation methods for missing binary data
#' 
#' @description
#' Imputation methods for the meta-analysis of binary outcomes with
#' missing data.
#' 
#' @param x An object of class \code{metabin}.
#' @param miss.e Number of missing observations in experimental group.
#' @param miss.c Number of missing observations in control group.
#' @param IMOR.e IMOR in experimental group (see Details).
#' @param IMOR.c IMOR in control group (see Details).
#' @param method.miss A character string indicating which method is
#'   used to impute missing values. Either \code{"GH"}, \code{"IMOR"},
#'   \code{"0"}, \code{"1"}, \code{"pc"}, \code{"pe"}, \code{"p"},
#'   \code{"b"}, or \code{"w"}, can be abbreviated (see Details).
#' @param small.values A character string specifying whether small
#'   treatment effects indicate a beneficial (\code{"good"}) or
#'   harmful (\code{"bad"}) effect, can be abbreviated (see Details).
#' @param fixed A logical indicating whether a fixed effect
#'   meta-analysis should be conducted.
#' @param random A logical indicating whether a random effects
#'   meta-analysis should be conducted.
#' @param prediction A logical indicating whether a prediction
#'   interval should be printed.
#'
#' @details
#' This function provides several imputation methods to deal with
#' missing data in the meta-analysis of binary outcomes (Gamble &
#' Hollis, 2005; Higgins et al., 2008). In order to utilise these
#' methods, the number of observations with missing outcomes must be
#' provided for the experimental and control group (arguments
#' \code{miss.e} and \code{miss.c}).
#'
#' The following imputation methods for missing binary data are available.
#' \tabular{ll}{
#' \bold{Argument}\tab \bold{Method} \cr 
#' \code{method.miss = "GH"}\tab Method by Gamble & Hollis (2005) \cr
#' \code{method.miss = "IMOR"}\tab Based on group-specific IMORs \cr
#' \code{method.miss = "0"}\tab Imputed as no events, (i.e., 0) \cr
#' \code{method.miss = "1"}\tab Imputed as events (i.e., 1) \cr
#' \code{method.miss = "pc"}\tab Based on observed risk in control group \cr
#' \code{method.miss = "pe"}\tab Based on observed risk in
#'   experimental group \cr
#' \code{method.miss = "p"}\tab Based on group-specific risks \cr
#' \code{method.miss = "b"}\tab Best case scenario for experimental group \cr
#' \code{method.miss = "w"}\tab Worst case scenario for experimental group
#' }
#'
#' The method by Gamble & Hollis (2005) is based on uncertainty
#' intervals for individual studies resulting from best and worst case
#' scenarios taking the missing data into account. The uncertainty
#' intervals are used to calculate (inflated) standard errors which
#' are considered in a generic inverse variance meta-analysis instead
#' of the standard errors from the complete case meta-analysis.
#'
#' All other methods are based on the Informative Missingness Odds
#' Ratio (IMOR) which is defined as the odds of an event in the
#' missing group over the odds of an event in the observed group
#' (Higgins et al., 2008). For example, an IMOR of 2 means that the
#' odds for an event is assumed to be twice as likely for missing
#' observations. For \code{method.miss = "IMOR"}, the IMORs in the
#' experimental (argument \code{IMOR.e}) and control group (argument
#' \code{IMOR.c}) must be specified by the user. For all other
#' methods, the input for arguments \code{IMOR.e} and \code{IMOR.c} is
#' ignored as these values are determined by the respective imputation
#' method (see Table 2 in Higgins et al., 2008).
#'
#' For the best and worst case scenarios (i.e., argument
#' \code{method.miss} equal to \code{"b"} or \code{"w"}), the user has
#' to specify whether the outcome is beneficial (argument
#' \code{small.values = "good"}) or harmful (\code{small.values =
#' "bad"}).
#' 
#' @return
#' An object of class \code{c("metamiss", "metagen", "meta")} with
#' corresponding \code{print}, \code{summary}, and \code{forest}
#' functions. See \code{\link[meta]{metagen}} for more information.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link[meta]{metabin}}, \code{\link[meta]{metagen}}
#' 
#' @references
#' Gamble C, Hollis S (2005):
#' Uncertainty method improved on best–worst case analysis in a binary
#' meta-analysis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{58}, 579--88
#' 
#' Higgins JPT, White IR, Wood AM (2008):
#' Imputation methods for missing outcome data in meta-analysis of
#' clinical trials.
#' \emph{Clinical Trials},
#' \bold{5}, 225--39
#' 
#' @examples
#' d1 <- data.frame(author = c("Beasley", "Selman"),
#'                  resp.h = c(29, 17), fail.h = c(18, 1), drop.h = c(22, 11),
#'                  resp.p = c(20, 7), fail.p = c(14, 4), drop.p = c(34, 18))
#' m1 <- metabin(resp.h, resp.h + fail.h, resp.p, resp.p + fail.p,
#'               data = d1, studlab = author, sm = "RR", method = "I")
#' m1
#'
#' # Treat missings as no events
#' metamiss(m1, drop.h, drop.p)
#' 
#' # Assume IMORs of 2 for both experimental and control group
#' metamiss(m1, drop.h, drop.p, IMOR.e = 2)
#'
#' # Gamble & Hollis (2005)
#' d2 <- data.frame(author = c("Lefevre", "van Vugt", "van Vugt"),
#'                  year = c(2001, 2000, 1998),
#'                  para.al = c(7, 4, 49), n.al = c(155, 134, 273),
#'                  miss.al = c(9, 16, 36),
#'                  para.ma = c(0, 0, 7), n.ma = c(53, 47, 264),
#'                  miss.ma = c(2, 3, 44))
#' 
#' m2 <- metabin(para.al, n.al, para.ma, n.ma,
#'               data = d2, studlab = paste0(author, " (", year, ")"),
#'               method = "Inverse", method.tau = "DL",
#'               sm = "OR")
#' 
#' metamiss(m2, miss.al, miss.ma, method = "GH")
#' @export metamiss
#'
#' @importFrom meta metabin


metamiss <- function(x,
                     miss.e, miss.c,
                     IMOR.e, IMOR.c = IMOR.e,
                     method.miss = if (missing(IMOR.e)) "0" else "IMOR",
                     small.values = "good",
                     fixed = x$fixed,
                     random = x$random,
                     prediction = x$prediction) {
  
  
  ##
  ##
  ## (1) Check for meta object and upgrade older meta objects
  ##
  ##
  chkclass(x, "metabin")
  x <- updateversion(x)
  ##
  ##if (!is.null(x$subgroup)) {
  ##  warning("Function metamiss() does not work with subgroup analyses.",
  ##          call. = FALSE)
  ##  return(NULL)
  ##}
  
  
  ##
  ## Catch 'miss.e' and 'miss.c' from data:
  ##
  mf <- match.call()
  miss.e <- eval(mf[[match("miss.e", names(mf))]],
                 x$data, enclos = sys.frame(sys.parent()))
  ##
  miss.c <- eval(mf[[match("miss.c", names(mf))]],
                 x$data, enclos = sys.frame(sys.parent()))
  ##
  if (is.null(miss.e))
    stop("Argument 'miss.e' missing.", call. = FALSE)
  ##
  if (is.null(miss.c))
    stop("Argument 'miss.c' missing.", call. = FALSE)
  ##
  if (isCol(x$data, ".subset")) {
    miss.e <- miss.e[x$data$.subset]
    miss.c <- miss.c[x$data$.subset]
  }
  ##
  event.e <- x$event.e
  event.c <- x$event.c
  n.e <- x$n.e
  n.c <- x$n.c
  ##
  chklength(miss.e, length(event.e),
            "metamiss",
            text = paste("Length of argument 'miss.e' must be equal to",
                         "number of studies in meta-analysis."))
  chklength(miss.c, length(event.e),
            "metamiss",
            text = paste("Length of argument 'miss.c' must be equal to",
                         "number of studies in meta-analysis."))
  ##
  incr <- 0.5 * (event.e == 0 | event.c == 0 | n.e == event.e | n.c == event.c)
  
  
  mm <- c("gh", "imor", "0", "1", "pc", "pe", "p", "b", "w")
  ##
  method.miss <- setchar(as.character(method.miss), mm)
  if (method.miss == "imor")
    method.miss <- "IMOR"
  if (method.miss == "gh")
    method.miss <- "GH"
  ##
  small.values <- setchar(small.values, c("good", "bad"))
  
  
  if (method.miss == "GH") {
    lower <- metabin(event.e, n.e + miss.e, event.c + miss.c, n.c + miss.c,
                     sm = x$sm, method.tau.ci = "")$lower
    upper <- metabin(event.e + miss.e, n.e + miss.e, event.c, n.c + miss.c,
                     sm = x$sm, method.tau.ci = "")$upper
    ##
    seTE <- TE.seTE.ci(lower, upper)$seTE
    ##
    res <- metagen(x$TE, seTE,
                   ##
                   studlab = x$studlab,
                   exclude = x$exclude,
                   ##
                   data = x$data,
                   ##
                   sm = x$sm,
                   ##
                   level = x$level, level.ma = x$level.ma,
                   fixed = fixed, random = random,
                   ##
                   hakn = x$hakn, method.tau = x$method.tau,
                   method.tau.ci = x$method.tau.ci,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   tau.common = x$tau.common,
                   ##
                   prediction = prediction, level.predict = x$level.predict,
                   ##
                   backtransf = x$backtransf,
                   title = x$title, complab = x$complab, outclab = x$outclab,
                   ##
                   label.e = x$label.e, label.c = x$label.c,
                   label.right = x$label.right, label.left = x$label.left,
                   ##
                   subgroup = x$subgroup, subgroup.name = x$subgroup.name,
                   print.subgroup.name = x$print.subgroup.name,
                   sep.subgroup = x$sep.subgroup,
                   ##
                   control = x$control)
    ##
    res$lower <- lower
    res$upper <- upper
    ##
    res$method.miss <- method.miss
    res$small.values <- small.values
    ##
    res$event.e <- event.e
    res$miss.e <- miss.e
    res$n.e <- n.e + miss.e
    ##
    res$event.c <- event.c
    res$miss.c <- miss.c
    res$n.c <- n.c + miss.c
  }
  else {
    ##
    p.e <- (event.e + incr) / (n.e + 1 * incr)
    p.c <- (event.c + incr) / (n.c + 1 * incr)
    k.All <- length(p.e)
    ##
    if (method.miss == "IMOR") {
      chknumeric(IMOR.e, min = 0)
      chknumeric(IMOR.c, min = 0)
      ##
      if (length(IMOR.e) == 1)
        IMOR.e <- rep(IMOR.e, k.All)
      if (length(IMOR.c) == 1)
        IMOR.c <- rep(IMOR.c, k.All)
      ##
      txt1 <- "Argument 'IMOR."
      txt2 <- paste("' must be of same length as number of",
                    "studies in meta-analysis or a single number.")
      txt.e <- paste0(txt1, "e", txt2)
      txt.c <- paste0(txt1, "c", txt2)
      chklength(IMOR.e, k.All, text = txt.e)
      chklength(IMOR.c, k.All, text = txt.c)
    }
    else if (method.miss == "0") {
      IMOR.e <- 0
      IMOR.c <- 0
    }
    ##
    if (method.miss == "1") {
      IMOR.e <- 1e8
      IMOR.c <- 1e8
    }
    ##
    if (method.miss == "pc") {
      IMOR.e <- p.c * (1 - p.e) / ((1 - p.c) * p.e)
      IMOR.c <- 1
    }
    ##
    if (method.miss == "pe") {
      IMOR.e <- 1
      IMOR.c <- p.e * (1 - p.c) / ((1 - p.e) * p.c)
    }
    ##
    if (method.miss == "p") {
      IMOR.e <- 1
      IMOR.c <- 1
    }
    ##
    if (method.miss == "b") {
      if (small.values == "good") {
        IMOR.e <- 0
        IMOR.c <- 1e8
      }
      else {
        IMOR.e <- 1e8
        IMOR.c <- 0
      }
    }
    ##
    if (method.miss == "w") {
      if (small.values == "good") {
        IMOR.e <- 1e8
        IMOR.c <- 0
      }
      else {
        IMOR.e <- 0
        IMOR.c <- 1e8
      }
    }
    ##
    pmiss.e <- miss.e / (n.e + miss.e)
    pmiss.c <- miss.c / (n.c + miss.c)
    ##
    p.star.e <-
      p.e * (1 - pmiss.e) +
      p.e * IMOR.e * pmiss.e / (p.e * IMOR.e + 1 - p.e)
    ##
    p.star.c <-
      p.c * (1 - pmiss.c) +
      p.c * IMOR.c * pmiss.c / (p.c * IMOR.c + 1 - p.c)
    
    
    var.p.star.e <-
      p.e * (1 - p.e) / n.e *
      (1 - pmiss.e + pmiss.e * IMOR.e / (p.e * IMOR.e + 1 - p.e)^2)^2 +
      pmiss.e * (1 - pmiss.e) / (n.e + miss.e) *
      (p.e * (1 - p.e) * (IMOR.e - 1) / (p.e * IMOR.e + 1 - p.e))^2
    ##
    var.p.star.c <-
      p.c * (1 - p.c) / n.c *
      (1 - pmiss.c + pmiss.c * IMOR.c / (p.c * IMOR.c + 1 - p.c)^2)^2 +
      pmiss.c * (1 - pmiss.c) / (n.c + miss.c) *
      (p.c * (1 - p.c) * (IMOR.c - 1) / (p.c * IMOR.c + 1 - p.c))^2
    
    
    if (x$sm == "RD") {
      TE.miss <- p.star.e - p.star.c
      varTE.miss <- var.p.star.e + var.p.star.c
    }
    else if (x$sm == "RR") {
      TE.miss <- log(p.star.e / p.star.c)
      varTE.miss <-
        var.p.star.e / p.star.e^2 +
        var.p.star.c / p.star.c^2
    }
    else if (x$sm == "OR") {
      TE.miss <- log((p.star.e / (1 - p.star.e)) / (p.star.c / (1 - p.star.c)))
      varTE.miss <-
        var.p.star.e / (p.star.e * (1 - p.star.e))^2 +
        var.p.star.c / (p.star.c * (1 - p.star.c))^2
    }
    
    
    res <- metagen(TE.miss, sqrt(varTE.miss),
                   ##
                   studlab = x$studlab,
                   exclude = x$exclude,
                   ##
                   data = x$data,
                   ##
                   sm = x$sm,
                   ##
                   level = x$level, level.ma = x$level.ma,
                   fixed = fixed, random = random,
                   ##
                   hakn = x$hakn, method.tau = x$method.tau,
                   method.tau.ci = x$method.tau.ci,
                   tau.preset = x$tau.preset, TE.tau = x$TE.tau,
                   tau.common = x$tau.common,
                   ##
                   prediction = prediction, level.predict = x$level.predict,
                   ##
                   backtransf = x$backtransf,
                   title = x$title, complab = x$complab, outclab = x$outclab,
                   ##
                   label.e = x$label.e, label.c = x$label.c,
                   label.right = x$label.right, label.left = x$label.left,
                   ##
                   subgroup = x$subgroup, subgroup.name = x$subgroup.name,
                   print.subgroup.name = x$print.subgroup.name,
                   sep.subgroup = x$sep.subgroup,
                   ##
                   control = x$control)
    ##
    res$event.e <- event.e
    res$miss.e <- miss.e
    res$n.e <- n.e + miss.e
    ##
    res$event.c <- event.c
    res$miss.c <- miss.c
    res$n.c <- n.c + miss.c
    ##
    res$IMOR.e <- IMOR.e
    res$IMOR.c <- IMOR.c
    ##
    res$method.miss <- method.miss
    res$small.values <- small.values
    ##
    res$incr <- incr
    res$p.e <- p.e
    res$p.c <- p.c
    ##
    res$pmiss.e <- pmiss.e
    res$pmiss.c <- pmiss.c
    ##
    res$p.star.e <- p.star.e
    res$p.star.c <- p.star.c
    ##
    res$var.p.star.e <- var.p.star.e
    res$var.p.star.c <- var.p.star.c
  }

  
  class(res) <- c("metamiss", class(res))
  ##
  res
}
