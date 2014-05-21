print.orbbound <- function(x,
                           comb.fixed=x$meta$comb.fixed, comb.random=x$meta$comb.random,
                           header=TRUE, logscale=x$meta$logscale,
                           digits=max(3, .Options$digits - 3),
                           ...){
  
  if (!inherits(x, "orbbound"))
    stop("Argument 'x' must be an object of class \"orbbound\"")
  
  k <- x$meta$k
  k.suspect <- x$k.suspect
  sm <- x$meta$sm
  
  if (length(logscale)==0)
    logscale <- FALSE
  
  if (sm=="ZCOR")
    sm.lab <- "COR"
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT"))
    sm.lab <- "proportion"
  else if (!logscale & sm == "PLN")
    sm.lab <- "proportion"
  else if (logscale & (sm == "RR" | sm == "OR" | sm == "HR" | sm == "IRR"))
    sm.lab <- paste("log", sm, sep="")
  else
    sm.lab <- sm
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  
  ci.lab <- paste(round(100*x$meta$level.comb, 1), "%-CI", sep="")
  
  TE.fixed    <- x$fixed$TE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  maxbias <- x$maxbias
  
  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR" | sm == "IRR" | sm=="PLN")){
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
    ##
    maxbias <- exp(maxbias)
  }
  else if (sm=="PAS"){
    TE.fixed    <- meta:::asin2p(TE.fixed, value="mean",
                                 warn=comb.fixed)
    lowTE.fixed <- meta:::asin2p(lowTE.fixed, value="lower",
                                 warn=comb.fixed)
    uppTE.fixed <- meta:::asin2p(uppTE.fixed, value="upper",
                                 warn=comb.fixed)
    ##
    TE.random    <- meta:::asin2p(TE.random, value="mean",
                                  warn=comb.random)
    lowTE.random <- meta:::asin2p(lowTE.random, value="lower",
                                  warn=comb.random)
    uppTE.random <- meta:::asin2p(uppTE.random, value="upper",
                                  warn=comb.random)
    ##
    maxbias <- meta:::asin2p(maxbias, value="mean", warn=comb.fixed|comb.random)
  }
  else if (sm=="PLOGIT"){
    TE.fixed    <- meta:::logit2p(TE.fixed)
    lowTE.fixed <- meta:::logit2p(lowTE.fixed)
    uppTE.fixed <- meta:::logit2p(uppTE.fixed)
    ##
    TE.random <- meta:::logit2p(TE.random)
    lowTE.random <- meta:::logit2p(lowTE.random)
    uppTE.random <- meta:::logit2p(uppTE.random)
    ##
    maxbias <- meta:::logit2p(maxbias)
  }
  else if (sm=="ZCOR"){
    TE.fixed    <- meta:::z2cor(TE.fixed)
    lowTE.fixed <- meta:::z2cor(lowTE.fixed)
    uppTE.fixed <- meta:::z2cor(uppTE.fixed)
    ##
    TE.random    <- meta:::z2cor(TE.random)
    lowTE.random <- meta:::z2cor(lowTE.random)
    uppTE.random <- meta:::z2cor(uppTE.random)
    ##
    maxbias <- meta:::z2cor(maxbias)
  }
  
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed <- x$fixed$p
  zTE.fixed <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random <- x$random$p
  zTE.random <- round(x$random$z, digits)
  ##
  maxbias <- round(maxbias, digits)
  
  if (header)
    meta:::crtitle(x$meta)
  
  if (comb.fixed|comb.random){
    if (!is.na(x$tau))
      tau2 <- ifelse(x$tau^2 < 0.0001, "tau^2 < 0.0001",
                     paste("tau^2=",
                           format(round(x$tau^2, 4), 4,
                                  nsmall=4, scientific=FALSE),
                           sep=""))
    
    cat("\n        Sensitivity Analysis for Outcome Reporting Bias (ORB)\n\n")
    
    cat(paste("Number of studies combined: k=", x$meta$k, "\n", sep=""))
    cat(paste("Between-study variance: ", tau2, "\n\n", sep=""))
  }
  
  if (comb.fixed){
    if (comb.random & x$tau==0)
      cat("Fixed effect model / Random effects model\n\n")
    else
      cat("Fixed effect model\n\n")
    
    res <- cbind(k.suspect,
                 format(maxbias),
                 format(TE.fixed),
                 meta:::rmSpace(meta:::p.ci(format(lowTE.fixed),
                                            format(uppTE.fixed))),
                 format(round(zTE.fixed, 4)),
                 meta:::format.p(pTE.fixed))
    
    zlab <- "z"
    
    names(res) <- c("k.suspect", "maxbias",
                    sm.lab, ci.lab, zlab, "p.value")
    
    dimnames(res) <- list(rep("", length(k.suspect)),
                          c("k.suspect", "maxbias",
                            sm.lab, ci.lab, zlab, "p.value"))
    
    prmatrix(res, quote=FALSE, right=TRUE)
  }
  
  if (comb.random & (x$tau!=0 | !comb.fixed)){
    cat("\nRandom effects model\n\n")
    
    res <- cbind(k.suspect,
                 format(maxbias),
                 format(TE.random),
                 meta:::rmSpace(meta:::p.ci(format(lowTE.random),
                                            format(uppTE.random))),
                 format(round(zTE.random, 4)),
                 meta:::format.p(pTE.random))
    
    zlab <- "z"
    
    names(res) <- c("k.suspect", "maxbias",
                    sm.lab, ci.lab, zlab, "p.value")
    
    dimnames(res) <- list(rep("", length(k.suspect)),
                          c("k.suspect", "maxbias",
                            sm.lab, ci.lab, zlab, "p.value"))
    
    prmatrix(res, quote=FALSE, right=TRUE)
  }
  
  if (comb.fixed|comb.random){
    ## Print information on summary method:
    meta:::catmeth(method=x$meta$method,
                   method.tau=if (comb.random) x$meta$method.tau else "",
                   sm=sm,
                   k.all=666)
  }
  
  invisible(NULL)
}
