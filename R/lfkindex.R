#' LFK Index Test for Asymmetry
#' 
#' Implementation of the LFK index test proposed by Furuya-Kanamori
#' et al. (2018) to evaluate bias in meta-analysis.
#' 
#' @param TE An object of class \code{meta} or estimated treatment
#'   effect in individual studies.
#' @param seTE Standard error of estimated treatment effect (mandatory
#'   if \code{TE} not of class \code{meta}).
#' @param data An optional data frame containing the study
#'   information.
#' @param x An object of class \code{lfkindex}.
#' @param digits Minimal number of significant digits, see
#'   \code{print.default}.
#' @param \dots Additional arguments (ignored).
#' 
#' @return
#' An object of class \code{"lfkindex"} with corresponding
#' \code{print} function. The object is a list containing the
#' following components:
#' 
#' \item{lfkindex}{LFK index.}
#' \item{interpretation}{Interpretation of value of LFK index.}
#' \item{abs.zscore}{Absolute value of z-score.}
#' \item{N, MidRank, percentile, zscore}{Quantities used to calculate
#'   LFK index.}
#' \item{TE, seTE}{Estimated treatment effect, standard error.}
#' \item{version}{Version of R package metasens used to create
#'   object.}
#' 
#' @author Gerta Rücker \email{ruecker@@imbi.uni-freiburg.de}, Guido
#'   Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{doiplot}}, \code{\link{metabias}},
#'   \code{\link{funnel.meta}}
#' 
#' @references
#' Furuya-Kanamori L, Barendregt JJ, Doi SAR (2018):
#' A new improved graphical and quantitative method for detecting bias
#' in meta-analysis.
#' \emph{International Journal of Evidence-Based Healthcare},
#' \bold{16}, 195--203
#' 
#' @examples
#' # Example from Furuya-Kanamori et al. (2018)
#' #
#' pain <- data.frame(SMD = c(-4.270, -1.710, -0.580, -0.190, 0.000),
#'                    varSMD = c(0.158,  0.076,  0.018,  0.022, 0.040))
#' 
#' lfk.pain <- lfkindex(SMD, sqrt(varSMD), data = pain)
#' lfk.pain
#'
#' doiplot(lfk.pain)
#' 
#' @export lfkindex


lfkindex <- function(TE, seTE, data = NULL) {
  
  
  ##
  ## Read data
  ##
  nulldata <- is.null(data)
  sfsp <- sys.frame(sys.parent())
  mc <- match.call()
  ##
  if (nulldata)
    data <- sfsp
  ##
  ## Catch 'TE' and 'seTE' from data:
  ##
  TE <- catch("TE", mc, data, sfsp)
  ##
  if (inherits(TE, "meta")) {
    seTE <- TE$seTE
    TE <- TE$TE
  }
  else {
    chknull(TE)
    ##
    seTE <- catch("seTE", mc, data, sfsp)
    chknull(seTE)
    ##
    chklength(seTE, length(TE), "lfkindex")
  }
  
  
  ##
  ## Calculate LFK index
  ##
  o <- order(TE)
  TE <- TE[o]
  seTE <- seTE[o]
  varTE <- seTE^2
  ##
  N <- 100 * max(varTE, na.rm = TRUE) / varTE
  MidRank <- vector(length = length(TE), mode = "numeric")
  MidRank[1] <- N[1] / 2 
  ##
  for (i in 2:length(TE)) {
    MidRank[i] <- MidRank[i - 1] + (N[i - 1] + N[i]) / 2
  }
  ##
  percentile <- (MidRank - 0.5) / sum(N)
  zscore <- qnorm(percentile)
  abs.zscore <- abs(zscore)
  ##
  TE.j <- TE[which.min(abs.zscore)]
  ##
  lfkindex <-
    5 / (2 * sum(!is.na(TE))) *
    sum(zscore +
        (max(zscore) - min(zscore)) /
        (max(TE - TE.j) - min(TE - TE.j)) *
        (TE - TE.j))
  ##
  interpretation <-
    if (abs(lfkindex) <= 1)
      "no asymmetry"
    else if (abs(lfkindex) <= 2)
      "minor asymmetry"
    else if (abs(lfkindex) > 2)
      "major asymmetry"
    else
      "unclear"
  
  
  res <- list(lfkindex = lfkindex, interpretation = interpretation,
              abs.zscore = abs.zscore,
              N = N, MidRank = MidRank,
              percentile = percentile, zscore = zscore,
              TE = TE, seTE = seTE,
              version = utils::packageDescription("metasens")$Version)
  
  
  class(res) <- "lfkindex"
  
  res
}





#' @rdname lfkindex
#' @method print lfkindex
#' @export
#' @export print.lfkindex


print.lfkindex <- function(x, digits = 2, ...) {
  
  chkclass(x, "lfkindex")
  
  cat("        LFK index test\n\n")

  cat(paste0("LFK index: ", round(x$lfkindex, digits),
             " (", x$interpretation, ")\n"))

  invisible(NULL)
}
