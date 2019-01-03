#' @importFrom stats dnorm pnorm


lambda <- function(x)
  dnorm(x) / pnorm(x)

