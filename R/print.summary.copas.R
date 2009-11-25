print.summary.copas <- function(x,
                                digits=max(3, .Options$digits - 3),
                                ...){
  
  sm <- x$sm
  
  if (sm == "RR" | sm == "OR" | sm == "HR"){
    x$slope$TE    <- exp(x$slope$TE)
    x$slope$lower <- exp(x$slope$lower)
    x$slope$upper <- exp(x$slope$upper)
    ##
    x$adjust$TE    <- exp(x$adjust$TE)
    x$adjust$lower <- exp(x$adjust$lower)
    x$adjust$upper <- exp(x$adjust$upper)
    ##
    x$random$TE    <- exp(x$random$TE)
    x$random$lower <- exp(x$random$lower)
    x$random$upper <- exp(x$random$upper)
  }  
  
  
  TE.adj    <- round(x$adjust$TE, digits)
  lowTE.adj <- round(x$adjust$lower, digits)
  uppTE.adj <- round(x$adjust$upper, digits)
  ##
  TE.random    <- round(x$random$TE, digits)
  lowTE.random <- round(x$random$lower, digits)
  uppTE.random <- round(x$random$upper, digits)
  
  
  TE.slope    <- round(x$slope$TE, digits)
  lowTE.slope <- round(x$slope$lower, digits)
  uppTE.slope <- round(x$slope$upper, digits)
  ##
  publprob <- x$publprob
  N.unpubl <- x$N.unpubl
  pval.TE  <- x$slope$p
  pval.rsb <- x$pval.rsb
  
  
  TE <- format(c(TE.slope, NA, TE.adj, TE.random), trim=TRUE)
  N.unpubl <- format(c(round(N.unpubl), NA,
                       round(x$N.unpubl.adj), NA), trim=TRUE)
  ##
  res <- cbind(c(format(round(publprob, 2)),
                 "", "Copas model (adj)","Random effects model"),
               ifelse(TE=="NA", "", TE),
               p.ci(format(c(lowTE.slope, NA, lowTE.adj, lowTE.random)),
                    format(c(uppTE.slope, NA, uppTE.adj, uppTE.random))),
               format.p(c(pval.TE, NA, x$adjust$p, x$random$p)),
               format.p(c(pval.rsb, NA, x$pval.rsb.adj, NA)),
               ifelse(N.unpubl=="NA", "", N.unpubl)
               )
  ##
  dimnames(res) <- list(rep("", dim(res)[[1]]),
                        c("publprob", sm, x$ci.lab,
                          "pval.treat", "pval.rsb", "N.unpubl"))
  
  
  cat("Summary of Copas selection model analysis:\n\n")
  ##
  prmatrix(res, quote=FALSE, right=TRUE, ...)
  ##
  cat("\n Legend:\n")
  cat(" publprob   - Probability of publishing the study with the largest standard error\n")
  cat(" pval.treat - P-value for hypothesis that the treatment effect is equal in both groups\n")
  cat(" pval.rsb   - P-value for hypothesis that no further selection remains unexplained\n")
  cat(" N.unpubl   - Approximate number of studies the model suggests remain unpublished\n")
  
  invisible(NULL)
}
