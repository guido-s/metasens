\name{forest.orbbound}

\alias{forest.orbbound}

\title{Forest plot}

\description{
  Draws a forest plot in the active graphics window (using grid
  graphics system).
}

\usage{
\method{forest}{orbbound}(x,
       comb.fixed=x$meta$comb.fixed,
       comb.random=x$meta$comb.random,
       text.fixed="FE model",
       text.random="RE model",
       smlab=NULL,
       leftcols=c("studlab", "maxbias"),
       leftlabs=c("Missing\nstudies", "Maximum\nbias"),
       logscale=x$meta$logscale,
       digits=max(3, .Options$digits - 3),
       ...)
}

\arguments{
  \item{x}{An object of class \code{orbbound}.}
  \item{comb.fixed}{A logical indicating whether sensitivity analysis
    for fixed effect model should be plotted.}
  \item{comb.random}{A logical indicating whether sensitivity analysis
    for random effects model should be plotted.}
  \item{text.fixed}{A character string used in the plot to label
    subgroup with results for fixed effect model.}
  \item{text.random}{A character string used in the plot to label
    subgroup with results for random effects model.}
  \item{smlab}{A label printed at top of figure. If only results for
    either fixed effect or random effects model is plotted, text
    indicates which model was used.}
  \item{leftcols}{A character vector specifying (additional) columns
    to be plotted on the left side of the forest plot or a logical
    value (see \code{\link{forest.meta}} help page for details).}
  \item{leftlabs}{A character vector specifying labels for
    (additional) columns on left side of the forest plot (see
    \code{\link{forest.meta}} help page for details).}
  \item{logscale}{A logical indicating whether results for summary
    measures 'RR', 'OR', and 'HR' will be printed on logarithmic
    scale.}
  \item{digits}{Minimal number of significant digits, see
    \code{print.default}.}
  \item{\dots}{Additional arguments for \code{\link{forest.meta}}
    function.}
}

\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window.

  For summary measures 'RR', 'OR', and 'HR' column labeled "Maximum
  bias" contains the relative bias, e.g. a value of 1.10 means a
  maximum overestimation by 10 percent. If \code{logscale=TRUE} for
  these summary measures, maximum bias is instead printed as absolute
  bias.

  For more information see help page of \code{\link{forest.meta}}
  function.
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{orbbound}}, \code{\link{print.orbbound}}}

\examples{
data(Fleiss93, package="meta")

meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")

orb1 <- orbbound(meta1, k.suspect=1:5)

print(orb1, digits=2)

forest(orb1, xlim=c(0.7, 1.5))
}

\keyword{hplot}
