print.limitmeta <- function(x,
                            sortvar,
                            logscale=FALSE,
                            digits=max(3, .Options$digits - 3),
                            header=TRUE, ...){
  
  if (!inherits(x, "limitmeta"))
    stop("Argument 'x' must be an object of class \"limitmeta\"")
  
  
  format.TE <- function(TE, na=FALSE){
    TE <- meta:::rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }
  
  
  sm <- x$sm
  ##
  if (logscale & (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR"))
    sm.lab <- paste("log", sm, sep="")
  else
    sm.lab <- sm
  
  
  ci.lab <- paste(round(100*x$level, 1), "%-CI", sep="")
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  
  
  k.all <- length(TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("Arguments 'x' and 'sortvar' have different length")
  
  
  ci.TE <- meta::ci(TE, seTE, level=x$level)
  lowTE <- ci.TE$lower
  uppTE <- ci.TE$upper
  ##
  ci.limit <- meta::ci(TE.limit, seTE.limit, level=x$level)
  lower.limit <- ci.limit$lower
  upper.limit <- ci.limit$upper
  
  
  if (!logscale & (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR")){
    TE    <- exp(TE)
    lowTE <- exp(lowTE)
    uppTE <- exp(uppTE)
    ##
    TE.limit <- exp(TE.limit)
    lower.limit <- exp(lower.limit)
    upper.limit <- exp(upper.limit)
  }
  ##
  TE    <- round(TE, digits)
  lowTE <- round(lowTE, digits)
  uppTE <- round(uppTE, digits)
  ##
  TE.limit    <- round(TE.limit, digits)
  lower.limit <- round(lower.limit, digits)
  upper.limit <- round(upper.limit, digits)
  ##
  TE    <- format.TE(TE, na=TRUE)
  lowTE <- format(lowTE, trim=TRUE)
  uppTE <- format(uppTE, trim=TRUE)
  ##
  TE.limit    <- format.TE(TE.limit, na=TRUE)
  lower.limit <- format(lower.limit, trim=TRUE)
  upper.limit <- format(upper.limit, trim=TRUE)
  
  
  res <- cbind(" ",
               TE,
               meta:::p.ci(lowTE, uppTE),
               "  ",
               TE.limit,
               meta:::p.ci(lower.limit, upper.limit))
  ##
  dimnames(res) <-
    list(x$studlab, c("", sm.lab, ci.lab, "", sm.lab, ci.lab))
  ##
  cat("Results for individual studies (left: original data; right: shrunken estimates)\n\n")
  prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
  
  
  cat("\n")
  
  
  print(summary(x),
        digits=digits, header=FALSE, logscale=logscale)
  
  
  invisible(NULL)
}
