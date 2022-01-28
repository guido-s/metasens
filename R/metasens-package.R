#' metasens: Brief overview of methods and general hints
#' 
#' R package \bold{metasens} provides advanced statistical methods to
#' model and adjust bias in meta-analysis and supports Schwarzer et
#' al. (2015), Chapter 5 "Small-Study Effects in Meta-Analysis"
#' \url{https://link.springer.com/book/10.1007/978-3-319-21416-0}.
#'
#' @name metasens-package
#' 
#' @docType package
#'
#' @details
#'
#' R package \bold{metasens} is an add-on package for \bold{meta}
#' providing the following meta-analysis methods:
#' \itemize{
#' \item Copas selection model (function \code{\link{copas}})
#'   described in Copas & Shi (2001) and evaluated in Schwarzer et
#'   al., 2010);
#' \item limit meta-analysis (\code{\link{limitmeta}}) by Rücker et
#'   al. (2011);
#' \item upper bound for outcome reporting bias
#'   (\code{\link{orbbound}}) described in Copas & Jackson (2004);
#' \item imputation methods for missing binary data
#'   (\code{\link{metamiss}}) described in Gamble & Hollis (2005) and
#'   Higgins et al. (2008).
#' }
#' 
#' Furthermore, functions and datasets from \bold{metasens} are
#' utilised in Schwarzer et al. (2015), Chapter 5
#' "Small-Study Effects in Meta-Analysis",
#' \url{https://link.springer.com/book/10.1007/978-3-319-21416-0}.
#' 
#' Type \code{help(package = "metasens")} for a listing of R functions
#' available in \bold{metasens}.
#' 
#' Type \code{citation("metasens")} on how to cite \bold{metasens} in
#' publications.
#' 
#' To report problems and bugs
#' \itemize{
#' \item type \code{bug.report(package = "metasens")} if you do not
#'   use RStudio,
#' \item send an email to Guido Schwarzer
#'   \email{sc@imbi.uni-freiburg.de} if you use RStudio.
#' }
#' 
#' The development version of \bold{metasens} is available on GitHub
#' \url{https://github.com/guido-s/metasens}.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @references
#' Copas J, Jackson D (2004):
#' A bound for publication bias based on the fraction of unpublished
#' studies.
#' \emph{Biometrics},
#' \bold{60}, 146--53
#' 
#' Copas JB, Shi JQ (2001):
#' A sensitivity analysis for publication bias in systematic reviews.
#' \emph{Statistical Methods in Medical Research},
#' \bold{10}, 251--65
#' 
#' Gamble C, Hollis S (2005):
#' Uncertainty method improved on best–worst case analysis in a binary
#' meta-analysis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{58}, 579--88
#' 
#' Higgins JPT, White IR, Wood AM (2008):
#' Imputation methods for missing outcome data in meta-analysis of
#' clinical trials.
#' \emph{Clinical Trials},
#' \bold{5}, 225--39
#' 
#' Rücker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M (2011):
#' Treatment-effect estimates adjusted for small-study effects via a limit
#' meta-analysis.
#' \emph{Biostatistics},
#' \bold{12}, 122--42
#' 
#' Schwarzer G, Carpenter J, Rücker G (2010):
#' Empirical evaluation suggests Copas selection model preferable to
#' trim-and-fill method for selection bias in meta-analysis.
#' \emph{Journal of Clinical Epidemiology},
#' \bold{63}, 282--8
#'
#' Schwarzer G, Carpenter JR, Rücker G (2015):
#' \emph{Meta-Analysis with R (Use-R!)}.
#' Springer International Publishing, Switzerland
#'
#' @keywords package
#'
#' @importFrom meta ci forest funnel gs metabias metabin metagen
#' @importFrom stats coefficients dnorm lm loess optim pchisq pnorm
#'   qlogis qnorm qt sd update var vcov


NULL
