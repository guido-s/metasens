\name{print.limitmeta}
\alias{print.limitmeta}
\alias{summary.limitmeta}
\alias{print.summary.limitmeta}
\title{Print and summary method for limit meta-analysis}
\description{
   Print and summary method for objects of class \code{limitmeta}.
}
\usage{
\method{print}{limitmeta}(x, sortvar, logscale=FALSE,
      digits=max(3, .Options$digits - 3), header=TRUE, ...)

\method{summary}{limitmeta}(object, ...)

\method{print}{summary.limitmeta}(x, logscale=FALSE,
      digits = max(3, .Options$digits - 3),
      header=TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{limitmeta} or \code{summary.limitmeta}.}
  \item{object}{An object of class \code{limitmeta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{logscale}{A logical indicating whether results for summary
    measures 'RR', 'OR', 'HR', and 'IRR' will be printed on
    logarithmic scale.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{...}{other arguments to the function will be ignored (this
    option included only to conform with R standards)}
}
\details{

  The \code{summary.limitmeta} function prints a summary of a limit
  meta-analysis (Rücker et al., 2011); unadjusted as well as adjusted effect estimates in a random effects model are printed.

  The function \code{print.limitmeta} prints the summary information above
  together with the following study information:

  \itemize{
   \item Effect estimate with confidence interval
   \item Shrunken effect estimates with confidence interval
  }
}
\value{
  The function \code{summary.limitmeta} returns the same list as the
  function \code{limitmeta}, however class "summary.limitmeta" is
  added to the object in order to print a short summary of the limit
  meta-analysis object.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{limitmeta}}, \code{\link{funnel.limitmeta}}, \code{\link[meta]{metabias}}, \code{\link[meta]{metagen}}}

\examples{
data(nsaids)
m1 <- metabin(Ee, Ne, Ec, Nc,
              data=nsaids, sm="OR", method="Inverse")

print(limitmeta(m1), digits=2)
}
\keyword{print}
