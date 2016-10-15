\name{limitmeta}
\alias{limitmeta}

\title{Limit meta-analysis}

\description{
  Implementation of the limit meta-analysis method by Rücker et
  al. (2011) to adjust for bias in meta-analysis.
}

\usage{limitmeta(x, method.adjust="beta0",
          level=x$level, level.comb=x$level.comb,
          backtransf=x$backtransf,
          title=x$title, complab=x$complab, outclab=x$outclab)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{method.adjust}{A character string indicating which adjustment
    method is to be used. One of \code{"beta0"}, \code{"betalim"}, or
    \code{"mulim"}, can be abbreviated.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If
    \code{backtransf=FALSE}, results for the odds ratio are printed as
    log odds ratios rather than odds ratio, for example.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
}

\details{
  This function provides the method by Rücker et al. (2011) to
  estimate an effect estimate adjusted for bias in meta-analysis. The
  underlying model is an extended random effects model that takes
  account of possible small study effects by allowing the treatment
  effect to depend on the standard error:

  theta(i) = beta + sqrt(SE(i)^2 + tau^2)(epsilon(i) + alpha),

  where epsilon(i) follows a standard normal distribution. Here
  theta(i) is the observed effect in study i, beta the global mean,
  SE(i) the within-study standard error, and tau^2 the between-study
  variance. The parameter alpha represents the bias introduced by
  small-study effects. On the one hand, alpha can be interpreted as
  the expected shift in the standardized treatment effect if precision
  is very small. On the other hand, theta(adj) = beta + tau*alpha is
  interpreted as the limit treatment effect for a study with infinite
  precision (corresponding to SE(i) = 0).

  Note that as alpha is included in the model equation, beta has a
  different interpretation as in the usual random effects model. The
  two models agree only if alpha=0. If there are genuine small-study
  effects, the model includes a component making the treatment effect
  depend on the standard error. The expected treatment effect of a
  study of infinite precision, beta + tau*alpha, is used as an
  adjusted treatment effect estimate.

  The maximum likelihood estimates for alpha and beta can be
  interpreted as intercept and slope in linear regression on a
  so-called generalised radial plot, where the x-axis represents the
  inverse of sqrt(SE(i)^2 + tau^2) and the y-axis represents the
  treatment effect estimates, divided by sqrt(SE(i)^2 + tau^2).

  Two further adjustments are available that use a shrinkage
  procedure. Based on the extended random effects model, a limit
  meta-analysis is defined by inflating the precision of each study
  with a common factor. The limit meta-analysis yields shrunken
  estimates of the study-specific effects, comparable to empirical
  Bayes estimates.  Based on the extended random effects model, we
  obtain three different treatment effect estimates that are adjusted
  for small-study effects:
  \itemize{
   \item an estimate based on the expectation of the extended random effects model, beta0 = beta + tau*alpha  (\code{method.adjust="beta0"})
   \item the extended random effects model estimate of the limit meta-analysis, including bias parameter (\code{method.adjust="betalim"})
   \item the usual random effects model estimate of the limit meta-analysis, excluding bias parameter (\code{method.adjust="mulim"})
  }
  
  See Rücker, Schwarzer et al. (2011), Section 7, for the definition
  of G^2 and the three heterogeneity statisticics \code{Q},
  \code{Q.small}, and \code{Q.resid}.
  
  For comparison, the original random effects meta-analysis is always
  printed in the sensitivity analysis.
}

\value{
  An object of class \code{"limitmeta"} with corresponding
  \code{print}, \code{summary} and \code{funnel} function. The object
  is a list containing the following components:

  \item{x, level, level.com}{}
  \item{method.adjust}{}
  \item{title, complab, outclab}{As defined above.}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{TE.limit, seTE.limit}{Shrunken estimates and standard error of individual studies.}
  \item{studlab}{Study labels.}
  \item{TE.random, seTE.random}{Unadjusted overall treatment effect and
    standard error (random effects model).}
  \item{lower.random, upper.random}{Lower and upper confidence interval limits
    (random effects model).}
  \item{zval.random, pval.random}{z-value and corresponding p-value
    for test of overall treatment effect (random effects model).}
  \item{w.random}{Weight of individual studies (in random effects model).}
  \item{tau}{Square-root of between-study variance.}
  \item{TE.adjust, seTE.adjust}{Adjusted overall effect and
    standard error (random effects model).}
  \item{lower.adjust, upper.adjust}{Lower and upper confidence
    interval limits for adjusted effect estimate (random effects
    model).}
  \item{zval.adjust, pval.adjust}{z-value and corresponding p-value
    for test of overall treatment effect for adjusted estimate (random
    effects model).}
  \item{alpha.r}{Intercept of the linear regression line on the
    generalised radial plot, here interpreted as bias parameter in an
    extended random effects model. Represents the expected shift in the
    standardized treatment effect if precision is very small.}
  \item{beta.r}{Slope of the linear regression line on the generalised
    radial plot.}
  \item{Q}{Heterogeneity statistic.}
  \item{Q.small}{Heterogeneity statistic for small study effects.}
  \item{Q.resid}{Heterogeneity statistic for residual heterogeneity
    beyond small study effects.}
  \item{G.squared}{Heterogeneity statistic G^2 (ranges from 0 to 100\%).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{call}{Function call.}
  \item{version}{Version of R package metasens used to create object.}
}

\references{
  Rücker G, Carpenter JR, Schwarzer G (2011),
  Detecting and adjusting for small-study effects in meta-analysis.
  \emph{Biometrical Journal}, \bold{53}, 351--68. 

  Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011),
  Treatment-effect estimates adjusted for small-study effects via a limit meta-analysis.
  \emph{Biostatistics}, \bold{12}, 122--42.
}

\author{Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{funnel.limitmeta}}, \code{\link{print.limitmeta}}}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(summary(limitmeta(m1)), digits=2)
}
