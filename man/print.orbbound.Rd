\name{print.orbbound}

\alias{print.orbbound}

\title{Print method for objects of class orbbound}

\description{
  Print method for objects of class \code{orbbound}.
}

\usage{
\method{print}{orbbound}(x,
      comb.fixed=x$x$comb.fixed, comb.random=x$x$comb.random,
      header=TRUE, backtransf=x$backtransf,
      digits=max(3, .Options$digits - 3),
      ...)
}

\arguments{
  \item{x}{An object of class \code{orbbound}.}
  \item{comb.fixed}{A logical indicating whether sensitivity analysis
    for fixed effect model should be printed.}
  \item{comb.random}{A logical indicating whether sensitivity analysis
    for random effects model should be printed.}
  \item{header}{A logical indicating whether information on
    meta-analysis should be printed at top of printout.}
  \item{backtransf}{A logical indicating whether printed results
    should be back transformed. If \code{backtransf=TRUE}, results for
    \code{sm="OR"} are printed as odds ratios rather than log odds
    ratios and results for \code{sm="ZCOR"} are printed as
    correlations rather than Fisher's z transformed correlations, for
    example.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{\dots}{Additional arguments}
}

\details{
  For summary measures 'RR', 'OR', and 'HR' column labeled maxbias
  contains the relative bias, e.g. a value of 1.10 means a maximum
  overestimation by 10 percent. If \code{logscale=TRUE} for these
  summary measures, maximum bias is instead printed as absolute bias.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{orbbound}}, \code{\link{forest.orbbound}}}

\examples{
data(Fleiss93, package="meta")

meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")

orb1 <- orbbound(meta1, k.suspect=1:5)

print(orb1, digits=2)

#
# Print log odds ratios instead of odds ratios
#
print(orb1, digits=2, backtransf=FALSE)

#
# Assuming that studies are missing on the left side
#
orb1.missleft <- orbbound(meta1, k.suspect=1:5, left=TRUE)

orb1.missleft

meta2 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR", method="Inverse")

orb2 <- orbbound(meta2, k.suspect=1:5)

print(orb2, digits=2)
}

\keyword{print}
