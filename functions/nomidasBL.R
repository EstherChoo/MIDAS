nomidas.bl <- function(horizon, y, z, iterations){
  for(l in 1:horizon){
    y <- as.matrix(y[-1,])
    z <- z[-nrow(z),]
  }
  
  #gamma prior and initialization - pre-MCMC
  sigma <- list()
  phi <- list()
  
  #initialise first values
  sigsq <- 3
  tau <- diag(1, ncol(z)) #tbc
  lam <- 1 #tbc
  
  for(iter in 1:iterations){
    #geting big phi
    zz <- t(z) %*% z
    sd.star <- sigsq * solve(solve(tau) + zz)
    phi.star <- solve(solve(tau) + zz) %*% t(z) %*% y
    
    ph <- mvrnorm(1, mu=phi.star, Sigma=sd.star) #draw from phi cond conj
    phi[[iter]] <- ph
    
    ###getting sigma (error) 
    a1 <- ncol(z) + length(y) - 1
    b1 <- t(y - z%*%ph) %*% (y - z%*%ph) + (t(ph) %*% solve(tau) %*% ph)
    
    sigsq <- rinvgamma(1, shape=a1/2, rate=b1/2)
    sigma[[iter]] <- sigsq
    
    #getting tau and lambda
    for(t in 1:ncol(z)){
      ig <- as.numeric(sqrt(lam * sigsq / ph[t]^2))
      tau[t,t] <- 1 / rinvgauss(1, ig, lam)
    }
    #adjust hyperparam delta and r such that prior density -> 0 as lam -> inf, 
    delta <- 1 #tbc
    r <- 1 #tbc
    lam <- rgamma(1, ncol(z) + r, sum(diag(tau)/2) + delta )
  }
  return(list(phi, sigma))
}
