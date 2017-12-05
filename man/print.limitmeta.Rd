\name{print.limitmeta}
\alias{print.limitmeta}
\alias{summary.limitmeta}
\alias{print.summary.limitmeta}
\title{Print and summary method for limit meta-analysis}
\description{
   Print and summary method for objects of class \code{limitmeta}.
}
\usage{
\method{print}{limitmeta}(x, sortvar, backtransf=x$backtransf,
      digits=gs("digits"), big.mark=gs("big.mark"),
      ...)

\method{summary}{limitmeta}(object, ...)

\method{print}{summary.limitmeta}(x, backtransf=x$backtransf,
      digits=gs("digits"),
      header=TRUE,
      pscale=x$x$pscale,
      irscale=x$x$irscale, irunit=x$x$irunit,
      digits.zval=gs("digits.zval"), digits.pval=gs("digits.pval"),
      digits.Q=gs("digits.Q"),
      digits.tau2=gs("digits.tau2"), digits.I2=gs("digits.I2"),
      scientific.pval=gs("scientific.pval"),
      big.mark=gs("big.mark"),
      print.Rb=gs("print.Rb"),
      warn.backtransf=FALSE, ...)
}

\arguments{
  \item{x}{An object of class \code{limitmeta} or \code{summary.limitmeta}.}
  \item{object}{An object of class \code{limitmeta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{backtransf}{A logical indicating whether results should be
    back transformed in printouts and plots. If \code{backtransf=TRUE}
    (default), results for \code{sm="OR"} are printed as odds ratios
    rather than log odds ratio, for example.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{digits.zval}{Minimal number of significant digits for z- or
    t-value, see \code{print.default}.}
  \item{digits.pval}{Minimal number of significant digits for p-value
    of overall treatment effect, see \code{print.default}.}
  \item{digits.Q}{Minimal number of significant digits for
    heterogeneity statistic Q, see \code{print.default}.}
  \item{digits.tau2}{Minimal number of significant digits for
    between-study variance, see \code{print.default}.}
  \item{digits.I2}{Minimal number of significant digits for I-squared
    and Rb statistic, see \code{print.default}.}
  \item{print.Rb}{A logical specifying whether heterogeneity statistic
    Rb should be printed.}
  \item{scientific.pval}{A logical specifying whether p-values should
    be printed in scientific notation, e.g., 1.2345e-01 instead of
    0.12345.}
  \item{big.mark}{A character used as thousands separator.} 
  \item{pscale}{A numeric giving scaling factor for printing of single
    event probabilities, i.e. if argument \code{sm} is equal to
    \code{"PLOGIT"}, \code{"PLN"}, \code{"PRAW"}, \code{"PAS"}, or
    \code{"PFT"}.}
  \item{irscale}{A numeric defining a scaling factor for printing of
    rates, i.e. if argument \code{sm} is equal to \code{"IR"},
    \code{"IRLN"}, \code{"IRS"}, or \code{"IRFT"}.}
  \item{irunit}{A character specifying the time unit used to calculate
    rates, e.g. person-years.}
  \item{warn.backtransf}{A logical indicating whether a warning should
    be printed if backtransformed proportions and rates are below 0
    and backtransformed proportions are above 1.}
  \item{\dots}{In \code{print.limitmeta}, additional arguments are
    passed on to \code{print.summary.limitmeta} called internally;
    otherwise, this argument is ignored.}
}

\details{
  The \code{summary.limitmeta} function prints a summary of a limit
  meta-analysis (Rücker et al., 2011); unadjusted as well as adjusted
  effect estimates in a random effects model are printed.

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
