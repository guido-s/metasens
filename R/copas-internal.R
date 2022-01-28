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


copas.loglik.without.beta <- function(x, gamma = c(-1.5, 0.08),
                                      TE, seTE) {
  
  
  mu  <- x[1]
  rho <- x[2]
  tau <- x[3]
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
  v <- (u + rho.tilde * (TE - mu) / (sqrt(s2t2))) /
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
  ell <- -(-0.5 * log(s2t2) - (TE - mu)^2 / (2 * s2t2) + log(pnorm(v)))
  
  res <- sum(ell)
  ##
  res
}
copas.gradient.without.beta <- function(x, gamma, TE, seTE) {
  
  
  mu  <- x[1]
  rho <- x[2]
  tau <- x[3]
  ##
  ## TE   <=> estimated treatment effect
  ## seTE <=> standard error from trials, conditional on publication
  ##
  TE.mu <- TE - mu
  rho2 <- rho^2
  tau2 <- tau^2
  varTE <- seTE^2
  
  
  ## Copas, Shi (2000), Biostatistics, p. 250:
  ##
  u <- gamma[1] + gamma[2] / seTE
  ##
  sigma <- sqrt(varTE / (1 - rho2 * lambda(u) * (u + lambda(u))))
  sigma2 <- sigma^2
  ##
  s2t2 <- sigma2 + tau2
  rho.tilde <- rho * sigma / sqrt(s2t2)
  rho2.tilde <- rho.tilde^2
  ##
  v <- (u + rho.tilde * TE.mu / (sqrt(s2t2))) / sqrt(1 - rho2.tilde)
  ##
  ## avoid numerical problems by replacing 0's in pnorm(v):
  ## qnorm(1e-320) = -38.26913
  ## this is towards the smallest value for log
  ##
  v[v < -37] <- -37
  ##
  ##
  ci2 <- lambda(u) * (u + lambda(u))
  
  
  ##
  ## Derivatives for mu:
  ##
  ## term 1 is zero
  ##
  ## term 2:
  ##
  grad.mu <- TE.mu / s2t2
  ##
  ## term 3 is always 0
  ##
  ## term 4:
  ##
  grad.mu <- grad.mu - (rho.tilde /
                        (sqrt((s2t2) *
                              (1 - rho2.tilde)))) * lambda(v)
  
  
  ##
  ## Derivatives for rho:
  ##
  ## term 1:
  ##
  grad.rho <- -ci2 * rho * sigma2 * sigma2 / (varTE * s2t2)
  ##
  ## term 2:
  ##
  grad.rho <- grad.rho +
    TE.mu^2 * ci2 * rho * sigma2 * sigma2 / (varTE * (s2t2^2))
  ##
  ## term 4:
  ##
  top <- u + rho.tilde * TE.mu / sqrt(s2t2)
  bottom <- sqrt((1 - rho2.tilde))
  ##
  diff.top <- (top - u) / rho - (top - u) * rho * ci2 / (1 - ci2 * rho2) +
               2 * (top - u) * rho * tau2 * ci2 / (varTE + tau2 * (1 - ci2 * rho2))
  ##
  eta <- varTE / (varTE + tau2 * (1 - ci2 * rho2))
  ##
  diff.bottom <- ((1 - rho2.tilde)^(-1.5)) * rho * eta *
                 (1 + (rho2 * tau2 * ci2 * eta) / varTE)
  ##
  grad.rho <- grad.rho + (top * diff.bottom + diff.top / bottom) * lambda(v)
  
  
  ##
  ## gradient for square-root of variance (tau)
  ##
  ## term 1:
  ##
  grad.tau <- -0.5 / s2t2
  ##
  ## term 2:
  ##
  grad.tau <- grad.tau + 0.5 * TE.mu^2 / s2t2^2
  ##
  ## term 4:
  ##
  grad.tau <- grad.tau +
              (-sigma * TE.mu *
                rho / (bottom * s2t2^2) -
                0.5 * top * ((1 - rho2.tilde)^(-1.5)) *
                sigma2 * rho2 / s2t2^2) *
              lambda(v)
  ##
  grad.tau <- 2 * tau * grad.tau
  
  
  ##  
  ## negative gradient as seek to minimise - log likelihood
  ##
  res <- c(sum(-grad.mu), sum(-grad.rho), sum(-grad.tau))
  ##
  res
}


lambda <- function(x)
  dnorm(x) / pnorm(x)

