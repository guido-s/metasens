radialregression <- function(TE, seTE, tau=0, k){

  seTE.tau <- sqrt(seTE^2 + tau^2)
  
  x <-  1/seTE.tau
  y <- TE/seTE.tau

  reg <- lm(y ~ x)

  intercept <- coefficients(reg)[1]
  slope     <- coefficients(reg)[2]

  se.slope     <- 1/sd(sqrt(1/seTE.tau^2))/sqrt(k-1)
  se.intercept <- se.slope*mean(1/seTE.tau^2)

  res <- list(intercept=intercept,
              se.intercept=se.intercept,
              slope=slope,
              se.slope=se.slope,
              k=k,
              r.squared=summary(reg)$r.squared)
  res
}
