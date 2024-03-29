% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metamiss.R
\name{metamiss}
\alias{metamiss}
\title{Imputation methods for missing binary data}
\usage{
metamiss(
  x,
  miss.e,
  miss.c,
  IMOR.e,
  IMOR.c = IMOR.e,
  method.miss = if (missing(IMOR.e)) "0" else "IMOR",
  small.values = "good",
  common = x$common,
  random = x$random,
  prediction = x$prediction,
  warn.deprecated = gs("warn.deprecated"),
  fixed
)
}
\arguments{
\item{x}{An object of class \code{metabin}.}

\item{miss.e}{Number of missing observations in experimental group.}

\item{miss.c}{Number of missing observations in control group.}

\item{IMOR.e}{IMOR in experimental group (see Details).}

\item{IMOR.c}{IMOR in control group (see Details).}

\item{method.miss}{A character string indicating which method is
used to impute missing values. Either \code{"GH"}, \code{"IMOR"},
\code{"0"}, \code{"1"}, \code{"pc"}, \code{"pe"}, \code{"p"},
\code{"b"}, or \code{"w"}, can be abbreviated (see Details).}

\item{small.values}{A character string specifying whether small
treatment effects indicate a beneficial (\code{"good"}) or
harmful (\code{"bad"}) effect, can be abbreviated (see Details).}

\item{common}{A logical indicating whether a common effect
meta-analysis should be conducted.}

\item{random}{A logical indicating whether a random effects
meta-analysis should be conducted.}

\item{prediction}{A logical indicating whether a prediction
interval should be printed.}

\item{warn.deprecated}{A logical indicating whether warnings should
be printed if deprecated arguments are used.}

\item{fixed}{Deprecated argument (replaced by 'common').}
}
\value{
An object of class \code{c("metamiss", "metagen", "meta")} with
corresponding \code{print}, \code{summary}, and \code{forest}
functions. See \code{\link[meta]{metagen}} for more information.
}
\description{
Imputation methods for the meta-analysis of binary outcomes with
missing data.
}
\details{
This function provides several imputation methods to deal with
missing data in the meta-analysis of binary outcomes (Gamble &
Hollis, 2005; Higgins et al., 2008). In order to utilise these
methods, the number of observations with missing outcomes must be
provided for the experimental and control group (arguments
\code{miss.e} and \code{miss.c}).

The following imputation methods for missing binary data are available.
\tabular{ll}{
\bold{Argument}\tab \bold{Method} \cr 
\code{method.miss = "GH"}\tab Method by Gamble & Hollis (2005) \cr
\code{method.miss = "IMOR"}\tab Based on group-specific IMORs \cr
\code{method.miss = "0"}\tab Imputed as no events, (i.e., 0) \cr
\code{method.miss = "1"}\tab Imputed as events (i.e., 1) \cr
\code{method.miss = "pc"}\tab Based on observed risk in control group \cr
\code{method.miss = "pe"}\tab Based on observed risk in
  experimental group \cr
\code{method.miss = "p"}\tab Based on group-specific risks \cr
\code{method.miss = "b"}\tab Best case scenario for experimental group \cr
\code{method.miss = "w"}\tab Worst case scenario for experimental group
}

The method by Gamble & Hollis (2005) is based on uncertainty
intervals for individual studies resulting from best and worst case
scenarios taking the missing data into account. The uncertainty
intervals are used to calculate (inflated) standard errors which
are considered in a generic inverse variance meta-analysis instead
of the standard errors from the complete case meta-analysis.

All other methods are based on the Informative Missingness Odds
Ratio (IMOR) which is defined as the odds of an event in the
missing group over the odds of an event in the observed group
(Higgins et al., 2008). For example, an IMOR of 2 means that the
odds for an event is assumed to be twice as likely for missing
observations. For \code{method.miss = "IMOR"}, the IMORs in the
experimental (argument \code{IMOR.e}) and control group (argument
\code{IMOR.c}) must be specified by the user. For all other
methods, the input for arguments \code{IMOR.e} and \code{IMOR.c} is
ignored as these values are determined by the respective imputation
method (see Table 2 in Higgins et al., 2008).

For the best and worst case scenarios (i.e., argument
\code{method.miss} equal to \code{"b"} or \code{"w"}), the user has
to specify whether the outcome is beneficial (argument
\code{small.values = "good"}) or harmful (\code{small.values =
"bad"}).
}
\examples{
d1 <- data.frame(author = c("Beasley", "Selman"),
  resp.h = c(29, 17), fail.h = c(18, 1), drop.h = c(22, 11),
  resp.p = c(20, 7), fail.p = c(14, 4), drop.p = c(34, 18))
m1 <- metabin(resp.h, resp.h + fail.h, resp.p, resp.p + fail.p,
  data = d1, studlab = author, sm = "RR", method = "I")
m1

# Treat missings as no events
metamiss(m1, drop.h, drop.p)

# Assume IMORs of 2 for both experimental and control group
metamiss(m1, drop.h, drop.p, IMOR.e = 2)

# Gamble & Hollis (2005)
d2 <- data.frame(author = c("Lefevre", "van Vugt", "van Vugt"),
  year = c(2001, 2000, 1998),
  para.al = c(7, 4, 49), n.al = c(155, 134, 273),
  miss.al = c(9, 16, 36),
  para.ma = c(0, 0, 7), n.ma = c(53, 47, 264),
  miss.ma = c(2, 3, 44))

m2 <- metabin(para.al, n.al, para.ma, n.ma,
  data = d2, studlab = paste0(author, " (", year, ")"),
  method = "Inverse", method.tau = "DL",
  sm = "OR")

metamiss(m2, miss.al, miss.ma, method = "GH")

}
\references{
Gamble C, Hollis S (2005):
Uncertainty method improved on best–worst case analysis in a binary
meta-analysis.
\emph{Journal of Clinical Epidemiology},
\bold{58}, 579--88

Higgins JPT, White IR, Wood AM (2008):
Imputation methods for missing outcome data in meta-analysis of
clinical trials.
\emph{Clinical Trials},
\bold{5}, 225--39
}
\seealso{
\code{\link[meta]{metabin}}, \code{\link[meta]{metagen}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
