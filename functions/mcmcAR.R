###mcmc ar function
#1: free 2: down 3:hump
mcmc.ar <- function(horizon, y, z.low, z.high, iterations, shape){
  for(l in 1:horizon){
    y <- as.matrix(y[-1,])
    z.low <- z.low[-nrow(z.low),]
    for(i in 1:length(z.high)){
      z.high[[i]] <- z.high[[i]][-nrow(z.high[[i]]),]
    }
  }
  
  #gamma prior and initialization - pre-MCMC
  theta <- list()
  sigma <- list()
  phi <- list()
  
  th0 <- list()
  for(i in 1:length(z.high)){
    if(shape == 2){
      th0[[i]] <- c(1, 0.5)
    } else {
    th0[[i]] <- c(0.5,0.5)
    }
  }

  z <- compress(th0, z.high, z.low, shape)
  
  sigsq <- 0.1
  
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
    
    ###metropolis for high freq
    c <- 4 #tuning param
    
    for(i in 1:length(z.high)){
      for(k in 1:2){
        if(shape == 2 & k==1){next}
        
        gam0 <- th0[[i]][k] #th0 is most updated version of theta
        gam1 <- rgamma(1, c*(gam0)^2, c*gam0)
        
        #probability ratio (denom)
        l0 <- sum(dnorm(y, c(z%*%ph), sqrt(sigsq), log=T)) #conditional likelihood w prev gamma
        
        prior0 <- dgamma(gam0, 1, 1, log=T) #prior 
        jump0 <- dgamma(gam1, c*(gam0^2), c*gam0, log=T) #new gamma given old
        
        #probability ratio (numer)
        oldz <- z
        newz <- z
        if(k==1){
          newg <- c(gam1, th0[[i]][2])
        } else {
          newg <- c(th0[[i]][1], gam1)
        }
        newz.high <- midaspoly(newg, z.high[[i]], shape)
        newz[,(ncol(newz)-4+i)] <- newz.high
        
        l1 <- sum(dnorm(y, c(newz%*%ph), sqrt(sigsq), log=T))
        
        prior1 <- dgamma(gam1, 1, 1, log=T)
        jump1 <- dgamma(gam0, c*(gam1^2), c*gam1, log=T) #old gamma given new
        
        prob <- exp((l1 + prior1 + jump1) - (l0 + prior0 + jump0))
        
        if(runif(1) <= prob){
          th0[[i]][k] <- gam1
          z <- newz
        }
        theta[[iter]] <- th0
        names(theta[[iter]]) <- c("Drt", "meanT", "minT", "maxT")
      }
    }
  }
  return(list(theta, phi, sigma))
}
