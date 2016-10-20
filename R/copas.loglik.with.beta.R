copas.loglik.with.beta <- function(x, gamma = c(-1.5, 0.08),
                                   TE, seTE) {
  
  mu   <- x[1]
  rho  <- x[2]
  tau  <- x[3]
  beta <- x[4]
  ##
  ## TE   <=> estimated treatment effect
  ## seTE <=> standard error from trials, conditional on publication
  
  
  ## Copas, Shi (2000), Biostatistics, p. 250:
  ##
  u <- gamma[1] + gamma[2] / seTE
  ##
  sigma <- sqrt(seTE^2 / (1 - rho^2 * lambda(u) * (u + lambda(u))))
  ##
  s2t2 <- sigma^2 + tau^2
  ##
  rho.tilde <- rho * sigma / sqrt(s2t2)
  ##
  v <- (u + rho.tilde * (TE - mu - beta * seTE) / (sqrt(s2t2))) /
       sqrt(1 - rho.tilde^2)
  ##
  ## Avoid numerical problems by replacing 0's in pnorm(v):
  ## qnorm(1e-320) = -38.26913
  ## this is towards the smallest value for log
  ##
  v[v < -37] <- -37
  ##
  ## Take minus log-likelihood and minimise it;
  ## leave out log(pnorm(u)) as this is a constant
  ##
  ell <- -(-0.5 * log(s2t2) -
         (TE - mu - beta * seTE)^2 / (2 * s2t2) + log(pnorm(v)))
  
  res <- sum(ell)
  ##
  res
}
