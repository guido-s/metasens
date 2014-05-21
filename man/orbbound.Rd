\name{orbbound}
\alias{orbbound}

\title{Sensitivity Analysis for Outcome Reporting Bias (ORB)}

\description{
  Implementation of the method by Copas & Jackson (2004) to evaluate
  outcome reporting bias in meta-analysis. An upper bound for outcome
  reporting bias is estimated for a given number of studies suspected
  with outcome reporting bias.
}

\usage{orbbound(x, k.suspect=1, tau=x$tau)
}

\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{k.suspect}{Number of studies with suspected outcome reporting bias.}
  \item{tau}{Square-root of between-study variance
    tau-squared.}
}

\details{
  This function provides the method by Copas and Jackson (2004) to
  estimate an upper bound for bias for a given number of studies with
  suspected outcome reporting bias.

  Based on the upper bound of outcome reporting bias, treatment
  estimates and confidence limits adjusted for bias are calculated.

  For comparison, the original meta-analysis is always considered in
  the sensitivity analysis (i.e. value 0 is always added to
  \code{k.suspect}).
}

\value{
  An object of class \code{c("orbbound")} with corresponding
  \code{print} and \code{forest} function. The object is a list
  containing the following components:
  \item{k.suspect, tau}{As defined above.}
  \item{maxbias}{Maximum bias for given values of \code{k.suspect}.}
  \item{fixed}{Adjusted treatment estimates and corresponding
    quantities for fixed effect model (a list with elements TE, seTE,
    lower, upper, z, p, level, df).}
  \item{random}{Adjusted treatment estimates and corresponding
    quantities for random effects model (a list with elements TE,
    seTE, lower, upper, z, p, level, df).}
  \item{meta}{Meta-analysis object (i.e. argument \code{x} from function call).}
  \item{call}{Function call.}
  \item{version}{Version of R package metasens used to create object.}
}

\references{
  Copas J, Jackson D (2004),
  A bound for publication bias based on the fraction of unpublished studies.
  \emph{Biometrics}, \bold{60}, 146--153.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{forest.orbbound}}, \code{\link{print.orbbound}}}

\examples{
data(Fleiss93, package="meta")

meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")

orb1 <- orbbound(meta1, k.suspect=1:5)

print(orb1, digits=2)

forest(orb1, xlim=c(0.7, 1.5))
}
