funnel.limitmeta <- function(x,
                             from,
                             to,
                             ##
                             pch=21,
                             cex=1,
                             col="black",
                             lwd=1,
                             ##
                             pch.adjust=18,
                             cex.adjust=1.5,
                             col.adjust="gray",
                             lty.line=1,
                             col.line="gray",
                             lwd.line=2,
                             ...){
  
  
  if (!inherits(x, "limitmeta"))
    stop("Argument 'x' must be an object of class \"limitmeta\"")
  
  
  if (x$method.adjust!="beta0"){
    warning(paste("TODO - Funnel plot for adjustment method ''",
                  x$method.adjust, "''?", sep=""))
    return(NULL)
  }
  
  
  TE <- x$TE
  TE.limit <- x$TE.limit
  TE.adjust <- x$TE.adjust
  seTE.tau <- sqrt(1/x$w.random)
  tau <- x$tau
  alpha.r <- x$alpha.r
  beta.r <- x$beta.r
  sm <- x$sm
  
  
  if (missing(from)){
    if (beta.r < 0)
      from <- x$TE.adjust+0.01
    if (beta.r > 0)
      from <- x$TE.adjust-0.01
    if (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR")
      from <- exp(from)
  }
  ##
  if (missing(to)){
    if (beta.r < 0)
      to <- max(TE)
    if (beta.r > 0)
      to <- min(TE)
    if (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR")
      to <- exp(to)
  }
  
  
  ##
  ## Generate funnel plot
  ##
  funnel(x$x, cex=cex, col=col, lwd=lwd, pch=pch, ...)
  
  
  ##
  ## Add line for adjustment method beta0
  ##
  if (sm=="RR" | sm=="OR" | sm=="HR" | sm=="IRR"){
    curve(sqrt((log(x)-beta.r)^2 / alpha.r^2 - tau^2),
          from=from, to=to,
          lty=lty.line, col=col.line, lwd=lwd.line, add=TRUE)
    points(exp(TE.adjust), 0, pch=pch.adjust, cex=cex.adjust, col=col.adjust)
  }
  else{
    curve(sqrt((x-beta.r)^2 / alpha.r^2 - tau^2),
          from=from, to=to,
          lty=lty.line, col=col.line, lwd=lwd.line, add=TRUE)
    points(TE.adjust, 0, pch=pch.adjust, cex=cex.adjust, col=col.adjust)
  }
  
  
  invisible(NULL)
}
