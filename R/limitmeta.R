limitmeta <- function(x,
                      method.adjust="beta0",
                      sm=x$sm,
                      level=x$level, level.comb=x$level.comb,
                      title=x$title, complab=x$complab, outclab=x$outclab){
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  

  TE <- x$TE
  seTE <- x$seTE
  tau <- x$tau
  w.random <- x$w.random
  k <- x$k
  ##
  seTE.tau  <- sqrt(1/w.random)
  ##
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  lower.random <- x$lower.random
  upper.random <- x$upper.random
  zval.random <- x$zval.random
  pval.random <- x$pval.random
  
  
  imeth <- charmatch(tolower(method.adjust),
                     c("beta0", "betalim", "mulim"), nomatch = NA)
  ##
  if(is.na(imeth))
    stop("Argument 'method.adjust' should be \"beta0\", \"betalim\", or \"mulim\"")
  ##
  method.adjust <- c("beta0", "betalim", "mulim")[imeth]
  
  
  ##
  ## Radial plot, slope best fit (beta-F)
  ##
  ## reg.f <- radialregression(TE, seTE, k)
  
  
  ##
  ## Generalized radial plot, slope best fit (beta-R)
  ##
  reg.r <- radialregression(TE, seTE.tau, k)
  ##
  alpha.r <- reg.r$intercept
  beta.r  <- reg.r$slope
  
  
  ##
  ## Conduct limit meta-analysis
  ##
  TE.limit <- beta.r + sqrt(tau^2/seTE.tau^2)*(TE - beta.r)
  seTE.limit <- seTE / 1 # 1 == "Infinity"
  ##
  m.lim <- metagen(TE.limit, seTE.limit, sm=sm)
  
  
  ##
  ##
  ## Conduct adjustment methods
  ##
  ##
  if (method.adjust=="beta0"){
    ##
    ## Expectation (beta-0)
    ##
    TE.adjust   <- as.vector(beta.r + tau*alpha.r)
    seTE.adjust <- as.vector(1/sd(sqrt(1/seTE^2))/sqrt(k-1))
  }
  else{
    if (method.adjust=="mulim"){
      ##
      ## Limit radial plot, slope through origin (mu-lim)
      ##
      TE.adjust   <- m.lim$TE.fixed
      seTE.adjust <- m.lim$seTE.fixed
    }
    else if (method.adjust=="betalim"){
      ##
      ## Limit radial plot, slope best fit (beta-lim)
      ##
      reg.l <- radialregression(m.lim$TE, m.lim$seTE, k)
      ##
      TE.adjust   <- as.vector(reg.l$slope)
      seTE.adjust <- as.vector(reg.l$se.slope)
    }
  }
  ##
  ci.adjust <- meta::ci(TE.adjust, seTE.adjust, level=level.comb)
  ##
  lower.adjust <- ci.adjust$lower
  upper.adjust <- ci.adjust$upper
  zval.adjust <- ci.adjust$z
  pval.adjust <- ci.adjust$p
  
  
  ##
  ## Only recalculate RE confidence interval if argument 'level.comb'
  ## is not missing
  ##
  if (!missing(level.comb)){
    ci.r <- meta::ci(TE.random, seTE.random, level=level.comb)
    ##
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }
  
  
  ## cat("\n****************\n")
  ## cat("*** G-square ***\n")
  ## cat("****************\n")
  ## GSquare <- 1-reg.l$r.squared
  ## print(GSquare)
  

  res <- list(TE=TE,
              seTE=seTE,
              ##
              TE.limit=TE.limit,
              seTE.limit=seTE.limit,
              ##
              studlab=x$studlab,
              ##
              TE.random=TE.random,
              seTE.random=seTE.random,
              lower.random=lower.random,
              upper.random=upper.random,
              zval.random=zval.random,
              pval.random=pval.random,
              w.random=w.random,
              tau=tau,
              ##
              TE.adjust=TE.adjust,
              seTE.adjust=seTE.adjust,
              lower.adjust=lower.adjust,
              upper.adjust=upper.adjust,
              zval.adjust=zval.adjust,
              pval.adjust=pval.adjust,
              ##
              alpha.r=alpha.r,
              beta.r=beta.r,
              ##
              level=level,
              level.comb=level.comb,
              ##
              k=k,
              sm=sm,
              method.adjust=method.adjust,
              ##
              title=title,
              complab=complab,
              outclab=outclab,
              ##
              call=match.call(),
              x=x)
  
  class(res) <- c("limitmeta")
  
  res
}
