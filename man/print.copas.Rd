\name{print.copas}
\alias{print.copas}
\alias{summary.copas}
\alias{print.summary.copas}
\title{Print and summary method for Copas selection model}
\description{
   Print and summary method for objects of class \code{copas}.
}
\usage{
\method{print}{copas}(x, sign.rsb=0.1, digits=max(3, .Options$digits - 3), ...)

\method{summary}{copas}(object, level=0.95, sign.rsb=0.1, ...)

\method{print}{summary.copas}(x,
                              digits = max(3, .Options$digits - 3),
                              header=TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{copas} or \code{summary.copas}.}
  \item{object}{An object of class \code{copas}.}
  \item{level}{The level used to calculate confidence intervals (between
    0 and 1).}
  \item{sign.rsb}{The significance level for the test of residual
    selection bias (between 0 and 1).}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{header}{A logical indicating whether information on title of
    meta-analysis, comparison and outcome should be printed at the
    beginning of the printout.}
  \item{...}{other arguments to the function will be ignored (this
    option included only to conform with R standards)}
}
\details{
  The \code{summary.copas} function prints a summary of a Copas
  analysis, performed using the function \code{copas}. It complements
  the graphical summary of the results, generated using
  \code{plot.copas}.

  Specifically it prints a table where the:

  first column corresponds to
  the x-axis in plots 3 & 4 from \code{plot.copas};

  second column corresponds to
  the treatment effect displayed in plot 3 from \code{plot.copas};

  third and fourth columns give
  the confidence intervals for this treatment effect,

  fifth colum gives
  the p-value for an overall treatment effect,

  sixth column gives
  the p-value for residual publication bias (the y-axis of plot 4 from
  \code{plot.copas} (see help(plot.copas) under plot 4 for a further
  explanation of this p-value))

  seventh column gives
  an approximate estimate of the number of studies the model suggests
  remain unpublished if the probability of publishing the study with the
  largest SE is as in column 1.

  Below this is displayed the results of the Copas analysis for the
  smallest degree of selection for which the p-value for evidence of
  residual selection bias exceeds \code{sign.rsb} (default: 0.1). This
  is simply extracted from the corresponding row in the table above.

  Lastly, the usual random effects estimate (based on the
  DerSimonian-Laird method) and 95\% confidence interval
  is printed.

  The function \code{print.copas} prints the summary information above
  together with the following information:

  Range of gamma0 values used (see \code{help(copas)});

  Range of gamma1 values used (see \code{help(copas)});

  Largest SE of all studies in meta-analysis;

  Range of probability publishing trial with largest SE;

  The next table gives details relating to the summary of the contour
  plot. Specifically, it gives details from fitting a straight line
  to each treatment-contour in the contour plot. Column 1 (headed
  level) shows the treatment-contours; column 2 (nobs) shows the number
  of observations used by the contour plot command within the
  \code{copas} function to plot this contour line; column 3
  (adj.r.square) shows the adjusted r-square from fitting a straight
  line to this contour; columns 4 & 5 show the slope and its standard
  error from fitting a straight line to this contour.
}
\value{
  A list is returned by the function \code{summary.copas} with the
  following elements:
  \item{slope}{Results for points on orthogonal line (a list with
    elements TE, seTE, lower, upper, z, p, level).}
  \item{publprob}{Vector of probabilities of publishing the smallest
    study.}
  \item{pval.rsb}{P-values for tests on presence of residual selection
    bias}
  \item{N.unpubl}{Approximate number of studies the model suggests
    remain unpublished}
  \item{adjust}{Result of Copas selection model adjusted for selection
    bias (a list with elements TE, seTE, lower, upper, z, p, level).}
  \item{sign.rsb}{The significance level for the test of residual
    selection bias.}
  \item{pval.rsb.adj}{P-value for test on presence of residual selection
    bias for adjusted effect given in \code{adjust}.}
  \item{N.unpubl.adj}{Approximate number of studies the model suggests
    remain unpublished for adjusted effect given in \code{adjust}}
  \item{random}{Results for usual random effects model (a list with
    elements TE, seTE, lower, upper, z, p, level).}
  \item{sm}{A character string indicating underlying summary measure.}
  \item{ci.lab}{Label for confidence interval.}
  \item{title}{Title of meta-analysis / systematic review.}
  \item{complab}{Comparison label.}
  \item{outclab}{Outcome label.}
  \item{version}{Version of R package copas used to create object.}
}
\author{James Carpenter \email{James.Carpenter@lshtm.ac.uk}, Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{copas}}, \code{\link{plot.copas}} \code{\link[meta]{metabias}}, \code{\link[meta]{metagen}}}

\examples{
# Load data
#
data(Fleiss93)
#
# Perform meta analysis, effect measure is odds ratio (OR)
#
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")
#
# Perform Copas analysis
#
cop1 <- copas(meta1)
summary(cop1, level=0.95)
}
\keyword{print}
