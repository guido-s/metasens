#' Limit meta-analysis
#' 
#' Implementation of the limit meta-analysis method by Rücker et al.
#' (2011) to adjust for bias in meta-analysis.
#' 
#' This function provides the method by Rücker et al. (2011) to
#' estimate an effect estimate adjusted for bias in meta-analysis. The
#' underlying model is an extended random effects model that takes
#' account of possible small study effects by allowing the treatment
#' effect to depend on the standard error:
#' 
#' theta(i) = beta + sqrt(SE(i)^2 + tau^2)(epsilon(i) + alpha),
#' 
#' where epsilon(i) follows a standard normal distribution. Here
#' theta(i) is the observed effect in study i, beta the global mean,
#' SE(i) the within-study standard error, and tau^2 the between-study
#' variance. The parameter alpha represents the bias introduced by
#' small-study effects. On the one hand, alpha can be interpreted as
#' the expected shift in the standardized treatment effect if
#' precision is very small. On the other hand, theta(adj) = beta +
#' tau*alpha is interpreted as the limit treatment effect for a study
#' with infinite precision (corresponding to SE(i) = 0).
#' 
#' Note that as alpha is included in the model equation, beta has a
#' different interpretation as in the usual random effects model. The
#' two models agree only if alpha=0. If there are genuine small-study
#' effects, the model includes a component making the treatment effect
#' depend on the standard error. The expected treatment effect of a
#' study of infinite precision, beta + tau*alpha, is used as an
#' adjusted treatment effect estimate.
#' 
#' The maximum likelihood estimates for alpha and beta can be
#' interpreted as intercept and slope in linear regression on a
#' so-called generalised radial plot, where the x-axis represents the
#' inverse of sqrt(SE(i)^2 + tau^2) and the y-axis represents the
#' treatment effect estimates, divided by sqrt(SE(i)^2 + tau^2).
#' 
#' Two further adjustments are available that use a shrinkage
#' procedure. Based on the extended random effects model, a limit
#' meta-analysis is defined by inflating the precision of each study
#' with a common factor. The limit meta-analysis yields shrunken
#' estimates of the study-specific effects, comparable to empirical
#' Bayes estimates.  Based on the extended random effects model, we
#' obtain three different treatment effect estimates that are adjusted
#' for small-study effects: \itemize{ \item an estimate based on the
#' expectation of the extended random effects model, beta0 = beta +
#' tau*alpha (\code{method.adjust="beta0"}) \item the extended random
#' effects model estimate of the limit meta-analysis, including bias
#' parameter (\code{method.adjust="betalim"}) \item the usual random
#' effects model estimate of the limit meta-analysis, excluding bias
#' parameter (\code{method.adjust="mulim"}) }
#' 
#' See Rücker, Schwarzer et al. (2011), Section 7, for the definition
#' of G^2 and the three heterogeneity statisticics \code{Q},
#' \code{Q.small}, and \code{Q.resid}.
#' 
#' For comparison, the original random effects meta-analysis is always
#' printed in the sensitivity analysis.
#' 
#' @param x An object of class \code{meta}.
#' @param method.adjust A character string indicating which adjustment
#'   method is to be used. One of \code{"beta0"}, \code{"betalim"}, or
#'   \code{"mulim"}, can be abbreviated.
#' @param level The level used to calculate confidence intervals for
#'   individual studies.
#' @param level.ma The level used to calculate confidence intervals
#'   for pooled estimates.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=FALSE}, results for the odds ratio are printed
#'   as log odds ratios rather than odds ratio, for example.
#' @param title Title of meta-analysis / systematic review.
#' @param complab Comparison label.
#' @param outclab Outcome label.
#' 
#' @return An object of class \code{"limitmeta"} with corresponding
#'   \code{print}, \code{summary} and \code{funnel} function. The
#'   object is a list containing the following components:
#' 
#' \item{x, level, level.ma,method.adjust,title, complab,
#'   outclab}{As defined above.}
#' \item{TE, seTE}{Estimated treatment effect and standard error of
#'   individual studies.}
#' \item{TE.limit, seTE.limit}{Shrunken estimates and standard error
#'   of individual studies.}
#' \item{studlab}{Study labels.}
#' \item{TE.random, seTE.random}{Unadjusted overall treatment effect
#'   and standard error (random effects model).}
#' \item{lower.random, upper.random}{Lower and upper confidence
#'   interval limits (random effects model).}
#' \item{statistic.random, pval.random}{Statistic and corresponding
#'   p-value for test of overall treatment effect (random effects
#'   model).}
#' \item{w.random}{Weight of individual studies (in random effects
#'   model).}
#' \item{tau}{Square-root of between-study variance.}
#' \item{TE.adjust, seTE.adjust}{Adjusted overall effect and standard
#'   error (random effects model).}
#' \item{lower.adjust, upper.adjust}{Lower and upper confidence
#'   interval limits for adjusted effect estimate (random effects
#'   model).}
#' \item{statistic.adjust, pval.adjust}{Statistic and corresponding
#'   p-value for test of overall treatment effect for adjusted
#'   estimate (random effects model).}
#' \item{alpha.r}{Intercept of the linear regression line on the
#'   generalised radial plot, here interpreted as bias parameter in an
#'   extended random effects model. Represents the expected shift in
#'   the standardized treatment effect if precision is very small.}
#' \item{beta.r}{Slope of the linear regression line on the
#'   generalised radial plot.} \item{Q}{Heterogeneity statistic.}
#' \item{Q.small}{Heterogeneity statistic for small study effects.}
#' \item{Q.resid}{Heterogeneity statistic for residual heterogeneity
#'   beyond small study effects.}
#' \item{G.squared}{Heterogeneity statistic G^2 (ranges from 0 to
#'   100\%).}
#' \item{k}{Number of studies combined in meta-analysis.}
#' \item{call}{Function call.} \item{version}{Version of R package
#'   metasens used to create object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#' Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{funnel.limitmeta}}, \code{\link{print.limitmeta}}
#' 
#' @references
#'
#' Rücker G, Carpenter JR, Schwarzer G (2011):
#' Detecting and adjusting for small-study effects in meta-analysis.
#' \emph{Biometrical Journal},
#' \bold{53}, 351--68
#' 
#' Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011):
#' Treatment-effect estimates adjusted for small-study effects via a limit
#' meta-analysis.
#' \emph{Biostatistics},
#' \bold{12}, 122--42
#' 
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'   data = Moore1998, sm = "OR", method = "Inverse")
#' 
#' print(limitmeta(m1), digits = 2)
#' @export limitmeta


limitmeta <- function(x,
                      method.adjust = "beta0",
                      level = x$level, level.ma = x$level.ma,
                      backtransf = x$backtransf,
                      title = x$title, complab = x$complab,
                      outclab = x$outclab) {
  
  chkclass(x, "meta")
  ##
  method.adjust <-
    setchar(method.adjust, c("beta0", "betalim", "mulim"))
  
  
  TE <- x$TE
  seTE <- x$seTE
  tau <- x$tau
  w.random <- x$w.random
  k <- x$k
  Q <- x$Q
  sm <- x$sm
  ##
  seTE.tau  <- sqrt(1 / w.random)
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  lower.random <- x$lower.random
  upper.random <- x$upper.random
  statistic.random <- x$statistic.random
  pval.random <- x$pval.random
  
  
  ##
  ## Radial plot, slope best fit (beta-F)
  ##
  reg.f <- radialregression(TE, seTE, k)
  
  
  ##
  ## Generalized radial plot, slope best fit (beta-R)
  ##
  reg.r <- radialregression(TE, seTE.tau, k)
  ##
  alpha.r <- reg.r$intercept
  beta.r  <- reg.r$slope
  
  
  ##
  ## Conduct limit meta-analysis
  ##
  TE.limit <- beta.r + sqrt(tau^2 / seTE.tau^2) * (TE - beta.r)
  seTE.limit <- seTE / 1 # 1 == "Infinity"
  ##
  m.lim <- metagen(TE.limit, seTE.limit, sm = sm, method.tau.ci = "")
  ##
  reg.l <- radialregression(m.lim$TE, m.lim$seTE, k)
  
  
  ##
  ##
  ## Conduct adjustment methods
  ##
  ##
  if (method.adjust == "beta0") {
    ##
    ## Expectation (beta-0)
    ##
    TE.adjust <- as.vector(beta.r + tau * alpha.r)
    ##
    var.beta <-
      as.vector(1 / var(sqrt(w.random)) / (k - 1))
    var.alpha <-
      as.vector(mean(w.random) / var(sqrt(w.random)) / (k - 1))
    cov.alpha.beta <-
      as.vector(- mean(sqrt(w.random)) / var(sqrt(w.random)) / (k - 1))
    ##
    seTE.adjust <- sqrt(var.beta + tau^2 * var.alpha + 2 * tau * cov.alpha.beta)
  }
  else {
    if (method.adjust == "mulim") {
      ##
      ## Limit radial plot, slope through origin (mu-lim)
      ##
      TE.adjust   <- m.lim$TE.common
      seTE.adjust <- m.lim$seTE.common
    }
    else if (method.adjust == "betalim") {
      ##
      ## Limit radial plot, slope best fit (beta-lim)
      ##
      TE.adjust   <- as.vector(reg.l$slope)
      seTE.adjust <- as.vector(reg.l$se.slope)
    }
  }
  ##
  ci.adjust <- ci(TE.adjust, seTE.adjust, level = level.ma)
  ##
  lower.adjust <- ci.adjust$lower
  upper.adjust <- ci.adjust$upper
  statistic.adjust <- ci.adjust$statistic
  pval.adjust <- ci.adjust$p
  ##
  if (inherits(x, c("metaprop"))) {
    statistic.adjust <- NA
    pval.adjust <- NA
  }
  
  
  ##
  ## Only recalculate RE confidence interval if argument 'level.ma'
  ## is not missing
  ##
  if (!missing(level.ma)) {
    ci.r <- ci(TE.random, seTE.random, level = level.ma)
    ##
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  
  
  ## Ruecker et al. (2011), Biostatistics, pages 133 - 4
  ##
  Q.resid <- reg.f$sigma^2 * (k - 2)
  Q.small <- Q - Q.resid
  ##
  G.squared = 1 - reg.l$r.squared
  
  
  res <- list(TE = TE,
              seTE = seTE,
              ##
              TE.limit = TE.limit,
              seTE.limit = seTE.limit,
              ##
              studlab = x$studlab,
              ##
              TE.random = TE.random,
              seTE.random = seTE.random,
              lower.random = lower.random,
              upper.random = upper.random,
              statistic.random = statistic.random,
              pval.random = pval.random,
              w.random = w.random,
              ##
              TE.adjust = TE.adjust,
              seTE.adjust = seTE.adjust,
              lower.adjust = lower.adjust,
              upper.adjust = upper.adjust,
              statistic.adjust = statistic.adjust,
              pval.adjust = pval.adjust,
              ##
              alpha.r = alpha.r,
              beta.r = beta.r,
              ##
              Q = Q,
              df.Q = x$df.Q,
              Q.small = Q.small,
              Q.resid = Q.resid,
              G.squared = G.squared,
              tau = tau,
              ##
              level = level,
              level.ma = level.ma,
              ##
              k = k,
              sm = sm,
              method.adjust = method.adjust,
              ##
              title = title,
              complab = complab,
              outclab = outclab,
              ##
              call = match.call(),
              x = x)
  
  
  res$backtransf <- backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  class(res) <- "limitmeta"
  
  res
}
