#' Doi plot for Asymmetry
#' 
#' Implementation of the Doi plot proposed by Furuya-Kanamori
#' et al. (2018) to evaluate bias in meta-analysis.
#' 
#' @param TE An object of class \code{lfkindex} or \code{meta} or
#'   estimated treatment effect in individual studies.
#' @param seTE Standard error of estimated treatment effect (mandatory
#'   if \code{TE} not of class \code{lfkindex} or \code{meta}).
#' @param xlim The x limits (min,max) of the plot.
#' @param ylim The y limits (min,max) of the plot.
#' @param xlab A label for the x-axis.
#' @param ylab A label for the y-axis.
#' @param lfkindex A logical indicating whether LFK index should be
#'   printed.
#' @param pos.lfkindex A character string with position of text with
#'   LFK index (see \code{\link{legend}}).
#' @param \dots Additional arguments (passed on to
#'   \code{\link{plot.default}}).
#' 
#' @author Gerta Rücker \email{gerta.ruecker@@uniklinik-freiburg.de}, Guido
#'   Schwarzer \email{guido.schwarzer@@uniklinik-freiburg.de}
#' 
#' @seealso \code{\link{lfkindex}}, \code{\link{metabias}},
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
#' @export doiplot
#'
#' @importFrom graphics legend


doiplot <- function(TE, seTE, xlim, ylim,
                    xlab = NULL, ylab = "|Z-score|",
                    lfkindex = TRUE, pos.lfkindex = "topleft",
                    ...) {
  
  
  if (!inherits(TE, "lfkindex")) {
    if (missing(seTE))
      lfk <- lfkindex(TE)
    else
      lfk <- lfkindex(TE, seTE)
  }
  else
    lfk <- TE
  ##
  if (missing(xlab)) {
    if (!is.null(lfk$x))
      xlab <- xlab(lfk$x$sm, FALSE)
    else
      xlab <- "ES"
  }
  chkchar(xlab, length = 1)
  ##
  chkchar(ylab, length = 1)
  chklogical(lfkindex)
  pos.lfkindex <-
    setchar(pos.lfkindex,
                   c("left", "center", "right",
                     "bottomleft", "bottom","bottomright",
                     "topleft", "top", "topright"))
  ##
  if (missing(xlim))
    xlim <- range(lfk$TE, na.rm = TRUE)
  ##
  if (missing(ylim))
    ylim <- c(max(lfk$abs.zscore), 0)
  
  
  plot(lfk$TE, lfk$abs.zscore, type = "b",
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       ...)
  ##
  if (lfkindex)
    legend(pos.lfkindex,
           paste("LFK index", round(lfk$lfkindex, 2)))
  
  invisible(NULL)
}
