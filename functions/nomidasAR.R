###mcmc ar function
nomidas.ar <- function(horizon, y, z, iterations){
  for(l in 1:horizon){
    y <- as.matrix(y[-1,])
    z <- z[-nrow(z),]
  }
  
  #gamma prior and initialization - pre-MCMC
  sigma <- list()
  phi <- list()
  
  sigsq <- 3
  
  for(iter in 1:iterations){ 
    #geting big phi
    zz <- t(z) %*% z
    mean0 <- rep(0, ncol(z))
    cov0 <- diag(10, ncol(z))
    sd.star <- solve(solve(cov0) + 1/sigsq * zz) #kxk matrix, prior precision + data precision
    
    phi.star <- sd.star %*% ((solve(cov0) %*% mean0)+ (1/sigsq * t(z) %*% y)) #1xk matrix, means weighted by prior n data precision
    
    ph <- mvrnorm(1, mu=phi.star, Sigma=sd.star) #draw from phi cond conj
    phi[[iter]] <- ph
    
    ###getting sigma (error) 
    a0 <- 1
    b0 <- 1
    
    a1 <- a0 + length(y) 
    b1 <- b0 + t(y - z%*%ph) %*% (y - z%*%ph)
    
    sigsq <- rinvgamma(1, shape=a1/2, rate=b1/2)
    sigma[[iter]] <- sigsq
    
  }
  return(list(phi, sigma))
}
