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
#' @param ... other arguments to the function will be ignored (this
#'   option included only to conform with R standards)
#'
#' @return An object of class "summary.copas" with corresponding print
#'   function. The object is a list containing the following
#'   components:
#'
#' \item{slope}{Results for points on orthogonal line (a list with
#'   elements TE, seTE, lower, upper, statistic, p, level).}
#' \item{publprob}{Vector of probabilities of publishing the smallest
#'   study.}
#' \item{pval.rsb}{P-values for tests on presence of residual
#'   selection bias}
#' \item{N.unpubl}{Approximate number of studies the model suggests
#'   remain unpublished}
#' \item{adjust}{Result of Copas selection model adjusted for
#'   selection bias (a list with elements TE, seTE, lower, upper,
#'   statistic, p, level).}
#' \item{sign.rsb}{The significance level for the test of residual
#'   selection bias.}
#' \item{pval.rsb.adj}{P-value for test on presence of residual
#'   selection bias for adjusted effect given in \code{adjust}.}
#' \item{N.unpubl.adj}{Approximate number of studies the model
#'   suggests remain unpublished for adjusted effect given in
#'   \code{adjust}}
#' \item{random}{Results for usual random effects model (a list with
#'   elements TE, seTE, lower, upper, statistic, p, level).}
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
#' data(Fleiss1993bin, package = "meta")
#' 
#' # Perform meta analysis, effect measure is odds ratio (OR)
#' #
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin, sm = "OR")
#' 
#' # Print summary of Copas analysis
#' #
#' summary(copas(m1, level.ma = 0.95))
#'
#' @method summary copas
#' @export
#' @export summary.copas


summary.copas <- function(object, ...) {
  
  
  chkclass(object, "copas")
  
  
  ci.random <- list(TE = object$TE.random,
                    seTE = object$seTE.random,
                    lower = object$lower.random,
                    upper = object$upper.random,
                    statistic = object$statistic.random,
                    p = object$pval.random,
                    level = object$level.ma)
  ##
  ci.slope <- list(TE = object$TE.slope,
                   seTE = object$seTE.slope,
                   lower = object$lower.slope,
                   upper = object$upper.slope,
                   statistic = object$statistic.slope,
                   p = object$pval.slope,
                   level = object$level.ma)
  ##
  ci.adjust <- list(TE = object$TE.adjust,
                    seTE = object$seTE.adjust,
                    lower = object$lower.adjust,
                    upper = object$upper.adjust,
                    statistic = object$statistic.adjust,
                    p = object$pval.adjust,
                    level = object$level.ma)
  ## 
  ci.lab <- paste(round(100 * object$level.ma, 1), "%-CI", sep = "")
  
  
  res <- list(slope = ci.slope,
              publprob = object$publprob,
              pval.rsb = object$pval.rsb,
              N.unpubl = object$N.unpubl,
              adjust = ci.adjust,
              sign.rsb = object$sign.rsb,
              pval.rsb.adj = object$pval.rsb.adj,
              N.unpubl.adj = object$N.unpubl.adj,
              random = ci.random,
              sm = object$sm,
              ci.lab = ci.lab,
              x = object
              )
  
  class(res) <- c("summary.copas")
  
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  
  res$backtransf <- object$backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  res
}
