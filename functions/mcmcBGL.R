###mcmc ar function
mcmc.bgl <- function(horizon, y, z.low, z.high, iterations, shape){
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
  z <- z[,-1]
  sigsq <- 0.1
  
  v <- diag(1, ncol(z))
  lam1 <- 1
  lam2 <- 1
  grp <- c(rep(1, ncol(z.low)), seq(from=2, length.out=length(z.high)))
  
  for(iter in 1:iterations){ 
    #geting big phi - coefficients
    zz <- t(z) %*% z
    sd.star <- solve(solve(v) + zz) #kxk matrix, prior precision + data precision
    
    phi.star <- sd.star %*% t(z) %*% y #1xk matrix, means weighted by prior n data precision
    
    ph <- mvrnorm(1, mu=phi.star, Sigma=sd.star) #draw from phi cond conj
    phi[[iter]] <- ph
    
    ###getting sigma (error)
    a1 <- ncol(z) + length(y)
    b1 <- t(y - z%*%ph) %*% (y - z%*%ph) + (t(ph) %*% solve(v) %*% ph)
    
    sigsq <- rinvgamma(1, shape=a1/2, rate=b1/2)
    sigma[[iter]] <- sigsq
    
    ###indiv lasso terms
    tau2 <- c()
    k <- 1
    for (j in ph){
      tau2[k] <- rinvgauss(n=1,
                           mean=sqrt((lam1*sigsq)/(j^2)) ,
                           shape=lam1)
      k<- k + 1
    }
    
    #group lasso terms
    gamma <- c()
    k<-1
    for (j in unique(grp)){
      ph_ind <- ph[which(grp==j)] #index the grouping, G is index of groups for beta
      gamma[which(grp==j)] <- rinvgauss(n=1,
                                      mean=sqrt(lam2*sigsq)/sum(ph_ind^2),
                                      shape=lam2)
      k<- k + 1
    }
    
    V <- diag(1/(1/tau2+1/gamma))

    ###
    
    lam1 <- rgamma(1, shape=length(tau2) + 1, sum(tau2)/2 + 0.1) #hyperparam for grp
    lam2 <- rgamma(1, shape=length(unique(grp))/2 + 1, sum(gamma)/2 + 0.1) #hyperparam for indiv
    
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
