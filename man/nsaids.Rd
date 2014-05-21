\name{nsaids}
\docType{data}
\alias{nsaids}
\title{NSAIDS in acute pain}

\description{Meta-analysis on the effectiveness of topical
  non-steroidal anti-inflammatory drugs (NSAIDS) in acute pain}

\usage{data(nsaids)}

\format{A data frame with the following columns:
  \describe{
  \item{study}{study number}
  \item{Ee}{number of treatment successes (reduction in pain of at
      least 50\%) in NSAIDS group}
  \item{Ne}{number of patients in NSAIDS group}
  \item{Ec}{number of treatment successes in control group}
  \item{Nc}{number of patients in control group} }
}

\source{
 Moore RA, Tramer MR, Carroll D, Wiffen PJ, McQuay HJ (1998),
 Quantitive systematic review of topically applied non-steroidal
 anti-inflammatory drugs.
 \emph{British Medical Journal}, \bold{316}, 333--338.
}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(limitmeta(m1), digits=2)
}

\keyword{datasets}
