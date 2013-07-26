orbbound <- function(x, k.suspect=1, tau){
  
  ## Copas J, Jackson D. A bound for publication bias based on the
  ## fraction of unpublished studies. Biometrics 2004,
  ## Mar;60(1):146-53.
  
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")
  
  if (missing(tau))
    tau <- x$tau

  if (!(is.numeric(k.suspect)))
    stop("Argument 'k.suspect' must be a numeric vector")
  ##
  if (any(k.suspect < 0))
    stop("Negative values not allowed for argument 'k.suspect'")
  ##
  if (!is.numeric(tau) || length(tau)!=1 || tau < 0)
    stop("Argument 'tau' must be a positive numeric of length 1")
  
  sel <- !is.na(x$seTE)
  
  if (x$hakn)
    warning("Hartung-Knapp adjustment not considered to evaluate outcome reporting bias")

  if (min(k.suspect)!=0)
    k.suspect <- c(0, k.suspect)
  ##
  k.suspect <- sort(k.suspect)
  
  maxbias <- ((x$k + k.suspect)/x$k *
              dnorm(qnorm(x$k/(x$k+k.suspect))) *
              sum(sqrt(x$seTE[sel]^2 + tau^2)^(-1)) /
              sum(sqrt(x$seTE[sel]^2 + tau^2)^(-2))
              )

  if (is.na(x$TE.fixed)){
    sign.f <- NA
    warning("Estimate from fixed effect model is 'NA'")
  }
  else if (x$TE.fixed<0)
    sign.f <- -1
  else if (x$TE.fixed>0)
    sign.f <-  1
  else if (x$TE.fixed==0){
    sign.f <- NA
    warning("Direction of outcome reporting bias is unclear as estimate from fixed effect model is equal to null effect")
  }
  ##
  if (is.na(x$TE.random)){
    sign.r <- NA
    warning("Estimate from random effects model is 'NA'")
  }
  else if (x$TE.random<0)
    sign.r <- -1
  else if (x$TE.random>0)
    sign.r <-  1
  else if (x$TE.random==0){
    sign.r <- NA
    warning("Direction of outcome reporting bias is unclear as estimate from random effects model is equal to null effect")
  }

  ci.f <- meta:::ci(x$TE.fixed-sign.f*maxbias, x$seTE.fixed)
  ci.r <- meta:::ci(x$TE.random-sign.r*maxbias, x$seTE.random)
  
  res <- list(maxbias=maxbias,
              k.suspect=k.suspect,
              tau=tau,
              fixed=ci.f,
              random=ci.r,
              meta=x,
              call=match.call(),
              version=packageDescription("copas")$Version)
  
  class(res) <- "orbbound"
  
  res
}
