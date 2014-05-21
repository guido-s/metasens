print.summary.limitmeta <- function(x,
                                    logscale=FALSE,
                                    digits=max(3, .Options$digits - 3),
                                    header=TRUE, ...){
  
  if (!inherits(x, "summary.limitmeta"))
    stop("Argument 'x' must be an object of class \"summary.limitmeta\"")
  
  
  sm <- x$sm
  ##
  if (logscale & (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR"))
    sm.lab <- paste("log", sm, sep="")
  else
    sm.lab <- sm
  
  
  TEs <- c(x$TE.adjust, x$TE.random)
  lower <- c(x$lower.adjust, x$lower.random)
  upper <- c(x$upper.adjust, x$upper.random)


  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR" | sm == "IRR")){
    TEs   <- exp(TEs)
    lower <- exp(lower)
    upper <- exp(upper)
  }
  ##
  TEs   <- round(TEs, digits)
  lower <- round(lower, digits)
  upper <- round(upper, digits)
  ##
  TEs   <- format(TEs, trim=TRUE)
  lower <- format(lower, trim=TRUE)
  upper <- format(upper, trim=TRUE)
  ##
  pvals <- c(x$pval.adjust, x$pval.random)
  zvals <- round(c(x$zval.adjust, x$zval.random), digits)
  
  
  imeth <- charmatch(tolower(x$method.adjust),
                     c("beta0", "betalim", "mulim"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("Argument 'method.adjust' should be \"beta0\", \"betalim\", or \"mulim\"")
  ##
  method.adjust.detail <- c("- expectation (beta0)",
                            "- including bias parameter (beta-lim)",
                            "- excluding bias parameter (mu-lim)")[imeth]
  
  
  ci.lab <- paste(round(100*x$level.comb, 1),
                  "%-CI", sep="")
  
  
  res <- cbind(c("Adjusted estimate",
                 "Unadjusted estimate"),
               ifelse(TEs=="NA", "", TEs),
               meta:::p.ci(lower, upper),
               zvals,
               meta:::format.p(pvals)
               )
  ##
  res[meta:::rmSpace(res)=="--"] <- ""
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("Random effects model",
                          sm.lab, ci.lab, "z", "pval"))
  
  
  if (header)
    meta:::crtitle(x)
  
  
  cat("Result of limit meta-analysis:\n\n")
  ##
  prmatrix(res, quote=FALSE, right=TRUE)
  ##
  cat("\nDetails on adjustment method:\n",
      method.adjust.detail,
      "\n", sep="")
  
  
  invisible(NULL)
}
