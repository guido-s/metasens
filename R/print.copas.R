print.copas <- function(x, sign.rsb=0.1,
                        digits=max(3, .Options$digits - 3),
                        ...){
  
  if (!inherits(x, "copas"))
    stop("Argument 'x' must be an object of class \"copas\"")
  
  
  meta:::crtitle(x)

  cat("Copas selection model analysis\n\n")
  
  
  res <- cbind(c("range of gamma0: ", "range of gamma1: "),
               format(c(round(x$gamma0.range[1], 2),
                        round(x$gamma1.range[1], 2))),
               ##               c(" - ", " - "),
               format(c(round(x$gamma0.range[2], 2),
                        round(x$gamma1.range[2], 2))))
  ##
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("", "min", "max"))
  ##
  prmatrix(res, quote=FALSE, right=TRUE)
  
  
  cat("\nLargest standard error (SE):", max(round(x$seTE, 4)), "\n\n")
  ##
  cat("Range of probability publishing trial with largest SE:\n")
  ##
  res <- matrix(format(round(range(pnorm(x$gamma0+x$gamma1/max(x$seTE))),
                             3)), nrow=1)
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("min", "max"))
  ##
  prmatrix(res, quote=FALSE, right=TRUE)

  cat("\n\nCalculation of orthogonal line:\n\n")
  ##
  res <- as.matrix(data.frame(x$regr)[,c("levels", "nobs",
                                         "adj.r.squareds",
                                         "slopes", "se.slopes")])
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("level", "nobs",
                          "adj.r.square",
                          "slope", "se.slope"))
  prmatrix(res, quote=FALSE, right=TRUE)
  
  cat("\n\n")
  print(summary(x, sign.rsb=0.1),
        digits=digits, header=FALSE)
  
  invisible(NULL)
}
