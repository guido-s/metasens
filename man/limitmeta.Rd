\name{limitmeta}
\alias{limitmeta}

\title{Sensitivity Analysis for Bias in Meta-Analysis}

\description{
  Implementation of the method by Rücker et al. (2011) to adjust for
  bias in meta-analysis.
}

\usage{limitmeta(x, method.adjust="beta0", sm=x$sm,
          level=x$level, level.comb=x$level.comb,
          title=x$title, complab=x$complab, outclab=x$outclab)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{method.adjust}{A character string indicating which adjustment
    method is to be used. One of \code{"beta0"}, \code{"betalim"}, or
    \code{"mulim"}, can be abbreviated.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"MD"},
    \code{"SMD"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
}

\details{
  This function provides the method by Rücker et al. (2011) to
  estimate an effect estimate adjusted for bias in meta-analysis.

  TODO - MORE DETAILS ON METHOD - GENERALISED RADIAL PLOT ...

  Three different adjustment methods are available (Rücker et al., 2011):
  \itemize{
   \item Expectation (\code{method.adjust="beta0"})
   \item Including bias parameter (\code{method.adjust="betalim"})
   \item Excluding bias parameter (\code{method.adjust="mulim"})
  }

  TODO - WHY DEFAULT beta0 ???

  For comparison, the original random effects meta-analysis is always
  printed in the sensitivity analysis.
}

\value{
  An object of class \code{c("limitmeta")} with corresponding
  \code{print}, \code{summary} and \code{funnel} function. The object
  is a list containing the following components:

  \item{x, level, level.com}{}
  \item{sm, method.adjust}{}
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
  \item{alpha.r}{TODO - DESCRIBE alpha.r}
  \item{beta.r}{TODO - DESCRIBE beta.r}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{call}{Function call.}
  \item{version}{Version of R package metasens used to create object.}
}

\references{
  Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011),
  Treatment-effect estimates adjusted for small-study effects via a limit meta-analysis.
  \emph{Biostatistics}, \bold{12}, 122--142.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker
\email{ruecker@imbi.uni-freiburg.de}}

\seealso{\code{\link{funnel.limitmeta}}, \code{\link{print.limitmeta}}}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(summary(limitmeta(m1)), digits=2)
}
