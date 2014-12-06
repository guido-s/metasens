\name{funnel.limitmeta}

\alias{funnel.limitmeta}

\title{Funnel plot for limit meta-analysis}

\description{
  Draws a funnel plot in the active graphics window.
}

\usage{
\method{funnel}{limitmeta}(x,
       pch=21, cex=1, col="black", bg="darkgray",
       lwd=1,
       pch.adjust=18, cex.adjust=1.5, col.adjust="gray", bg.adjust="gray",
       line=TRUE, xmin.line, xmax.line, 
       lty.line=1, lwd.line=lwd, col.line="gray",
       shrunken=FALSE, pch.shrunken=22, cex.shrunken=1,
       col.shrunken="black", bg.shrunken="white",
       lty.connect=1, lwd.connect=0.8, col.connect="black",
       backtransf=x$backtransf, ...)
}

\arguments{
  \item{x}{An object of class \code{limitmeta}.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{bg}{A vector with background colour of plotting symbols (only
    used if \code{pch} in \code{21:25}).}
  \item{lwd}{The line width for confidence intervals (see
    \code{\link{funnel.meta}}).}
  \item{pch.adjust}{The plotting symbol used for the adjusted effect estimate.}
  \item{cex.adjust}{The magnification to be used for the plotting
    symbol of the adjusted effect estimate.}
  \item{col.adjust}{Colour of plotting symbol for adjusted effect estimate.}
  \item{bg.adjust}{Background colour of plotting symbol for adjusted effect estimate.}
  \item{line}{A logical indicating whether adjusted regression line should be plotted.}
  \item{xmin.line}{Minimal value for the adjusted regression line (on x-axis).}
  \item{xmax.line}{Maximum value for the adjusted regression line (on x-axis).}
  \item{lty.line}{Line type of the adjusted regression line.}
  \item{col.line}{Color of the adjusted regression line.}
  \item{lwd.line}{The line width of the adjusted regression line.}
  \item{shrunken}{A logical indicating whether shrunken treatment
    estimates should be plotted.}
  \item{pch.shrunken}{The plotting symbol used for shrunken effect estimates.}
  \item{cex.shrunken}{The magnification to be used for the plotting
    symbol of the shrunken effect estimates.}
  \item{col.shrunken}{Colour of plotting symbol for shrunken effect estimates.}
  \item{bg.shrunken}{Background colour of plotting symbol for shrunken effect estimates.}
  \item{lty.connect}{Line type for line connecting original and
    shrunken treatment estimates.}
  \item{lwd.connect}{The line width of the connecting lines.}
  \item{col.connect}{Color of the connecting lines.}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If \code{backtransf=TRUE}
    (default), results for \code{sm="OR"} are printed as odds ratios
    rather than log odds ratio, for example.}
  \item{\dots}{Additional arguments for \code{\link{funnel.meta}}
    function.}
}

\details{
  A funnel plot is drawn in the active graphics window. In addition
  this function adds the adjusted effect estimate as well as a
  nonlinear regression line (also called adjusted regression line) if
  argument \code{line} is \code{TRUE}. The adjusted regression line is
  representing the dependence of the treatment effect estimate on the
  standard error across studies. The adjusted regression line is only
  plotted in addition to the adjusted treatment effect if argument
  \code{method.adjust="beta0"} (default) has been used in the
  \code{\link{limitmeta}} function.

  If argument \code{shrunken} is \code{TRUE} the shrunken effect
  estimates are also plotted. Lines are connecting original and
  shrunken effect estimates.

  Internally, R function \code{\link{funnel.meta}} is called to create
  a funnel plot. For more information see help page of the
  \code{\link{funnel.meta}} function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}}

\seealso{\code{\link{limitmeta}}, \code{\link[meta]{funnel.meta}}}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(summary(limitmeta(m1)), digits=2)

funnel(limitmeta(m1))

#
# Print results on log scale
#
print(summary(limitmeta(m1)), digits=2, backtransf=FALSE)
funnel(limitmeta(m1), backtransf=FALSE)
}

\keyword{hplot}
