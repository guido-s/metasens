#' Sensitivity Analysis for Outcome Reporting Bias (ORB)
#' 
#' Implementation of the method by Copas & Jackson (2004) to evaluate
#' outcome reporting bias in meta-analysis. An upper bound for outcome
#' reporting bias is estimated for a given number of studies suspected
#' with outcome reporting bias.
#' 
#' This function provides the method by Copas and Jackson (2004) to
#' estimate an upper bound for bias for a given number of studies with
#' suspected outcome reporting bias.
#' 
#' Based on the upper bound of outcome reporting bias, treatment
#' estimates and confidence limits adjusted for bias are calculated.
#' 
#' For comparison, the original meta-analysis is always considered in
#' the sensitivity analysis (i.e. value 0 is always added to
#' \code{k.suspect}).
#' 
#' @param x An object of class \code{meta}.
#' @param k.suspect Number of studies with suspected outcome reporting
#'   bias.
#' @param tau Square-root of between-study variance tau-squared.
#' @param left A logical indicating whether the cause of any selection
#'   bias is due to missing studies on the left or right of the funnel
#'   plot: left hand side if \code{left=TRUE}, right hand side if
#'   \code{left=FALSE}. If not set, the linear regression test for
#'   funnel plot asymmetry (i.e., function \code{metabias(...,
#'   meth="linreg")}) is used to determine whether studies are missing
#'   on the left or right hand side.
#' @param backtransf A logical indicating whether results should be
#'   back transformed in printouts and plots. If
#'   \code{backtransf=TRUE} (default), results for \code{sm="OR"} are
#'   printed as odds ratios rather than log odds ratio, for example.
#' 
#' @return An object of class \code{c("orbbound")} with corresponding
#'   \code{print} and \code{forest} function. The object is a list
#'   containing the following components:
#' \item{k.suspect, tau}{As defined above.}
#' \item{maxbias}{Maximum bias for given values of \code{k.suspect}.}
#' \item{common}{Adjusted treatment estimates and corresponding
#'   quantities for common effect model (a list with elements TE, seTE,
#'   lower, upper, z, p, level, df).}
#' \item{random}{Adjusted treatment estimates and corresponding
#'   quantities for random effects model (a list with elements TE,
#'   seTE, lower, upper, z, p, level, df).}
#' \item{left}{Whether selection bias expected on left or right}
#'   \item{x}{Meta-analysis object (i.e. argument \code{x} from
#'   function call).}
#' \item{call}{Function call.}
#' \item{version}{Version of R package metasens used to create
#'   object.}
#' 
#' @author Guido Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{forest.orbbound}}, \code{\link{print.orbbound}}
#' 
#' @references
#' 
#' Copas J, Jackson D (2004):
#' A bound for publication bias based on the fraction of unpublished
#' studies.
#' \emph{Biometrics},
#' \bold{60}, 146--53
#' 
#' @examples
#' data(Fleiss1993bin, package = "meta")
#' 
#' m1 <- metabin(d.asp, n.asp, d.plac, n.plac, data = Fleiss1993bin, sm = "OR")
#' 
#' orb1 <- orbbound(m1, k.suspect = 1:5)
#' print(orb1, digits = 2)
#' forest(orb1, xlim = c(0.75, 1.5))
#' 
#' # Same result
#' #
#' orb2 <- orbbound(m1, k.suspect = 1:5, left = FALSE)
#' print(orb2, digits = 2)
#' 
#' # Assuming bias in other direction
#' #
#' orb3 <- orbbound(m1, k.suspect = 1:5, left = TRUE)
#' print(orb3, digits = 2)
#' @export orbbound


orbbound <- function(x, k.suspect = 1, tau = x$tau, left = NULL,
                     backtransf = x$backtransf) {
  
  
  ## Copas J, Jackson D. A bound for publication bias based on the
  ## fraction of unpublished studies. Biometrics 2004,
  ## Mar;60(1):146-53.
  
  
  chkclass(x, "meta")
  x <- updateversion(x)
  
  if (!(is.numeric(k.suspect)))
    stop("Argument 'k.suspect' must be a numeric vector")
  ##
  if (any(k.suspect < 0))
    stop("Negative values not allowed for argument 'k.suspect'")
  ##
  if (!is.numeric(tau) || length(tau) != 1 || tau < 0)
    stop("Argument 'tau' must be a positive numeric of length 1")
  
  
  if (is.null(left))
    left <- as.logical(sign(metabias(x, meth = "linreg", k.min = 3)$estimate[1]) == 1)
  else if (!(length(left) == 1 & is.logical(left)))
    stop("Argument 'left' must be a logical of length 1.")
  
  
  sel <- !is.na(x$seTE)

  
  if (x$hakn)
    warning("Hartung-Knapp adjustment not considered to evaluate outcome reporting bias.")
  
  
  if (min(k.suspect) != 0)
    k.suspect <- c(0, k.suspect)
  ##
  k.suspect <- sort(k.suspect)
  
  
  if (left)
    direction.bias <- -1
  else
    direction.bias <-  1
  
  
  maxbias <- (x$k + k.suspect) / x$k *
    dnorm(qnorm(x$k / (x$k + k.suspect))) *
    sum(sqrt(x$seTE[sel]^2 + tau^2)^(-1)) /
    sum(sqrt(x$seTE[sel]^2 + tau^2)^(-2))
  
  ##
  maxbias <- direction.bias * maxbias
  
  
  ci.c <- ci(x$TE.common + maxbias, x$seTE.common, level = x$level.ma)
  ci.r <- ci(x$TE.random + maxbias, x$seTE.random, level = x$level.ma)
  
  
  res <- list(maxbias = maxbias,
              k.suspect = k.suspect,
              tau = tau,
              common = ci.c,
              random = ci.r,
              left = left,
              x = x,
              call = match.call())
  
  res$backtransf <- backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  class(res) <- "orbbound"
  
  res
}
