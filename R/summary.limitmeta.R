#' Summary method for limit meta-analysis
#' 
#' Summary method for objects of class \code{limitmeta}.
#' 
#' @aliases summary.limitmeta
#' 
#' @param object An object of class \code{limitmeta}.
#' @param \dots Additional arguments (ignored).
#' 
#' @return This function returns the same list as the function
#'   \code{limitmeta}, however class "summary.limitmeta" is added to
#'   the object in order to print a detailed summary of the limit
#'   meta-analysis object.
#' 
#' @author Guido Schwarzer \email{sc@@imbi.uni-freiburg.de}
#' 
#' @seealso \code{\link{limitmeta}}, \code{\link{funnel.limitmeta}},
#'   \code{\link{print.summary.limitmeta}}
#' 
#' @examples
#' data(Moore1998)
#' m1 <- metabin(succ.e, nobs.e, succ.c, nobs.c,
#'               data = Moore1998, sm = "OR", method = "Inverse")
#' 
#' summary(limitmeta(m1))
#'
#' @method summary limitmeta
#' @export
#' @export summary.limitmeta


summary.limitmeta <- function(object, ...) {

  chkclass(object, "limitmeta")
  
  res <- object
  class(res) <- c("summary.limitmeta", class(object))
  
  res
}
