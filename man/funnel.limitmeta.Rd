\name{funnel.limitmeta}

\alias{funnel.limitmeta}

\title{Funnel plot}

\description{
  Draws a funnel plot in the active graphics window.
}

\usage{
\method{funnel}{limitmeta}(x,
       pch=21, cex=1, col="black", lwd=1,
       pch.adjust=18, cex.adjust=1.5, col.adjust="gray",
       xmin.line, xmax.line, 
       lty.line=1, col.line="gray", lwd.line=2,
       ...)
}

\arguments{
  \item{x}{An object of class \code{limitmeta}.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{lwd}{The line width for confidence intervals.}
  \item{pch.adjust}{The plotting symbol used for the adjusted effect estimate.}
  \item{cex.adjust}{The magnification to be used for the plotting
    symbol of the adjusted effect estimate.}
  \item{col.adjust}{Colour of plotting symbol for adjusted effect estimate.}
  \item{xmin.line}{Minimal value for line TODO - HOW DO WE CALL THE LINE ???  (on x-axis).}
  \item{xmax.line}{Maximum value for line TODO - HOW DO WE CALL THE LINE ???  (on x-axis).}
  \item{lty.line}{Line type (TODO - HOW DO WE CALL THE LINE ???).}
  \item{col.line}{Color of line TODO - HOW DO WE CALL THE LINE ???.}
  \item{lwd.line}{The line width for TODO - HOW DO WE CALL THE LINE ???.}
  \item{\dots}{Additional arguments for \code{\link{funnel.meta}}
    function.}
}

\details{
  A funnel plot is drawn in the active graphics window. Internally, R
  function \code{\link{funnel.meta}} is called to create a funnel
  plot. In addition this function adds the adjusted effect estimate as
  well as a line representing TODO - HOW DO WE DESCRIBE THE CURVE ???.

  For more information see help page of \code{\link{funnel.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{limitmeta}}, \code{\link[meta]{funnel.meta}}}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(summary(limitmeta(m1)), digits=2)

funnel(limitmeta(m1))
}

\keyword{hplot}
