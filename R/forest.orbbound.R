forest.orbbound <- function(x,
                            comb.fixed=x$meta$comb.fixed,
                            comb.random=x$meta$comb.random,
                            text.fixed="FE model",
                            text.random="RE model",
                            smlab=NULL,
                            leftcols=c("studlab", "maxbias"),
                            leftlabs=c("Missing\nstudies", "Maximum\nbias"),
                            logscale=FALSE,
                            digits=max(3, .Options$digits - 3),
                            ...){

  if (!inherits(x, "orbbound"))
    stop("Argument 'x' must be an object of class \"orbbound\"")
  
  k <- x$meta$k
  sm <- x$meta$sm
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  ##
  if (length(logscale)==0)
    logscale <- FALSE
  
  if (sm=="ZCOR")
    sm.lab <- "COR"
  else if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT"))
    sm.lab <- "proportion"
  else if (!logscale & sm == "PLN")
    sm.lab <- "proportion"
  else if (logscale & (sm == "RR" | sm == "OR" | sm == "HR"))
    sm.lab <- paste("log", sm, sep="")
  else
    sm.lab <- sm
  
  ci.lab <- paste(round(100*x$meta$level.comb, 1), "%-CI", sep="")
  
  TE.fixed   <- x$fixed$TE
  seTE.fixed <- rep(x$fixed$seTE, length(TE.fixed))
  ##
  TE.random   <- x$random$TE
  seTE.random <- rep(x$random$seTE, length(TE.random))

  if (comb.fixed & comb.random){
    TE <- c(TE.fixed, TE.random)
    seTE <- c(seTE.fixed, seTE.random)
    FEvsRE <- c(rep(text.fixed, length(TE.fixed)),
                rep(text.random, length(TE.random)))
  }
  if (comb.fixed & !comb.random){
    TE <- TE.fixed
    seTE <- seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed effect model"
  }
  if (!comb.fixed & comb.random){
    TE <- TE.random
    seTE <- seTE.random
    if (is.null(smlab))
      smlab <- "Random effects model"
  }
  if (!comb.fixed & !comb.random){
    warning("No forest plot generated as both arguments 'comb.fixed' and 'comb.random' are FALSE")
    return(invisible(NULL))
  }

  if (comb.fixed & comb.random)
    m1 <- metagen(TE, seTE, sm=sm.lab,
                  byvar=FEvsRE, print.byvar=FALSE,
                  warn=FALSE)
  else
    m1 <- metagen(TE, seTE, sm=sm.lab, warn=FALSE)
  ##
  if (comb.fixed & comb.random){
    m1$studlab <- c(x$k.suspect, x$k.suspect)
    m1$maxbias <- c(x$maxbias, x$maxbias)
  }
  else{
    m1$studlab <- x$k.suspect
    m1$maxbias <- x$maxbias
  }

  if (!logscale & (sm == "RR" | sm == "OR" | sm == "HR" | sm=="PLN"))
    m1$maxbias <- format(round(exp(m1$maxbias), digits))
  else
    m1$maxbias <- format(round(m1$maxbias, digits))
  
  forest(m1,
         comb.fixed=FALSE, comb.random=FALSE,
         hetstat=FALSE,
         leftcols=leftcols,
         leftlabs=leftlabs,
         smlab=smlab,
         just.studlab="center",
         weight="same",
         ...)
  
  invisible(NULL)
}
