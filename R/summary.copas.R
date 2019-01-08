#' Summary method for Copas selection model
#' 
#' Summary method for objects of class \code{copas}.
#' 
#' This function complements the graphical summary of the results of a
#' Copas selection model, generated using \code{plot.copas}.
#'
#' @aliases summary.copas
#'
#' @param object An object of class \code{copas}.
#' @param level The level used to calculate confidence intervals
#'   (between 0 and 1).
#' @param sign.rsb The significance level for the test of residual
#'   selection bias (between 0 and 1).
#' @param ... other arguments to the function will be ignored (this
#'   option included only to conform with R standards)
#'
#' @return An object of class "summary.copas" with corresponding print
#'   function. The object is a list containing the following
#'   components:
#'
#' \item{slope}{Results for points on orthogonal line (a list with
#'   elements TE, seTE, lower, upper, z, p, level).}
#' \item{publprob}{Vector of probabilities of publishing the smallest
#'   study.}
#' \item{pval.rsb}{P-values for tests on presence of residual
#'   selection bias}
#' \item{N.unpubl}{Approximate number of studies the model suggests
#'   remain unpublished}
#' \item{adjust}{Result of Copas selection model adjusted for
#'   selection bias (a list with elements TE, seTE, lower, upper, z,
#'   p, level).}
#' \item{sign.rsb}{The significance level for the test of residual
#'   selection bias.}
#' \item{pval.rsb.adj}{P-value for test on presence of residual
#'   selection bias for adjusted effect given in \code{adjust}.}
#' \item{N.unpubl.adj}{Approximate number of studies the model
#'   suggests remain unpublished for adjusted effect given in
#'   \code{adjust}}
#' \item{random}{Results for usual random effects model (a list with
#'   elements TE, seTE, lower, upper, z, p, level).}
#' \item{sm}{A character string indicating underlying summary
#'   measure.}
#' \item{ci.lab}{Label for confidence interval.}
#' \item{title}{Title of meta-analysis / systematic review.}
#' \item{complab}{Comparison label.} \item{outclab}{Outcome label.}
#' \item{version}{Version of R package metasens used to create
#'   object.}
#'
#' @author James Carpenter \email{James.Carpenter@@lshtm.ac.uk}, Guido
#' Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{copas}}, \code{\link{plot.copas}},
#' \code{\link[meta]{metabias}}, \code{\link[meta]{metagen}}
#' 
#' @examples
#' data(Fleiss93)
#' 
#' # Perform meta analysis, effect measure is odds ratio (OR)
#' #
#' m1 <- metabin(event.e, n.e, event.c, n.c, data = Fleiss93, sm = "OR")
#' 
#' # Print summary of Copas analysis
#' #
#' summary(copas(m1), level = 0.95)
#'
#' @method summary copas
#' @export
#' @export summary.copas
#'
#' @importFrom meta ci


summary.copas <- function(object, level = 0.95,
                          sign.rsb = object$sign.rsb, ...) {
  
  meta:::chkclass(object, "copas")
  
  seTE <- object$seTE
  TE.random <- object$TE.random
  seTE.random <- object$seTE.random
  gamma0.slope <- object$gamma0.slope
  gamma1.slope <- object$gamma1.slope
  TE.slope <- object$TE.slope
  seTE.slope <- object$seTE.slope
  publprob <- object$publprob
  pval.rsb <- object$pval.rsb
  N.unpubl <- object$N.unpubl
  ##
  ci.random <- ci(TE.random, seTE.random, level)
  ##
  if (is.null(sign.rsb))
    sign.rsb <- 0.1
  else
    meta:::chklevel(sign.rsb)
  
  
  ci.lab <- paste(round(100 * level, 1), "%-CI", sep = "")
  
  
  ord <- rev(order(publprob)) 
  pom <- publprob[ord]
  ##
  TE.slope <- TE.slope[ord]
  seTE.slope <- seTE.slope[ord]
  pval.rsb <- pval.rsb[ord]
  N.unpubl <- N.unpubl[ord]
  ##
  ci.slope <- ci(TE.slope, seTE.slope, level)
  
  
  ##
  ## Copas estimate adjusted for selection bias (added by sc, 24.09.2007):
  ##
  tres <- data.frame(seq = seq(along = pval.rsb),
                     cumsum = cumsum(pval.rsb <= sign.rsb),
                     diff = seq(along = pval.rsb) - cumsum(pval.rsb <= sign.rsb))
  pval.rsb.sign.all <- all(tres$diff == 0)
  pval.rsb.sign <- ifelse(sum(tres$diff == 0) > 0, TRUE, FALSE)
  ##
  if (pval.rsb.sign.all) {
    TE.adj <- NA
    seTE.adj <- NA
    pval.rsb.adj <- NA
    N.unpubl.adj <- NA
  }
  else {
    if(pval.rsb.sign) {
      sel.adj <- tres$seq[tres$diff > 0][1]
      TE.adj <- TE.slope[sel.adj]
      seTE.adj <- seTE.slope[sel.adj]
      pval.rsb.adj <- pval.rsb[sel.adj]
      N.unpubl.adj <- N.unpubl[sel.adj]
    }
    else {
      TE.adj <- TE.slope[1]
      seTE.adj <- seTE.slope[1]
      pval.rsb.adj <- pval.rsb[1]
      N.unpubl.adj <- N.unpubl[1]
    }
  }
  ##
  adjust <- ci(TE.adj, seTE.adj, level)
  
  res <- list(slope = ci.slope,
              publprob = pom,
              pval.rsb = pval.rsb,
              N.unpubl = N.unpubl,
              adjust = adjust,
              sign.rsb = sign.rsb,
              pval.rsb.adj = pval.rsb.adj,
              N.unpubl.adj = N.unpubl.adj,
              random = ci.random,
              sm = object$sm,
              ci.lab = ci.lab
              )
  
  class(res) <- c("summary.copas")
  
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  
  res$backtransf <- object$backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  res
}
