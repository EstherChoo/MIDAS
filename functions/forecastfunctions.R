###summarise mcmc results into posterior mean and quantiles
summariseMCMC <- function(mcmc.results, burnin){
  fullsummary <- list()
  
  for(h in 1:length(mcmc.results)){
    #summarise phi
    summary <- list()
    phi <- mcmc.results[[h]][[2]]
    phi <- phi[-c(1:burnin)]
    for(col in 1:length(phi[[1]])){
      lst <- list()
      x <- unlist(lapply(phi, "[[", col))
      lst$quantile <- quantile(x, probs=c(0.025, 0.5, 0.975))
      lst$mean <- mean(x)
      lst$samples <- x
      summary$phi[[col]] <- lst
      names(summary$phi)[col] <- names(phi[[1]])[col]
    }
    
    #summarise sigma
    lst <- list()
    sigma <- mcmc.results[[h]][[3]][-c(1:burnin)]
    x <- unlist(sigma)
    lst$quantile <- quantile(x, probs=c(0.025, 0.5, 0.975))
    lst$mean <- mean(x)
    lst$samples <- x
    summary$sigma <- lst
    
    #summarise theta
    theta <- mcmc.results[[h]][[1]][-c(1:(burnin))]
    for(col in 1:length(theta[[1]])){
      x12 <- lapply(theta, "[[", col)
      lst <- list()
      for(param in 1:2){
        lst[[param]] <- list()
        x <- unlist(lapply(x12, "[", param))
        lst[[param]]$quantile <- quantile(x, probs=c(0.025, 0.5, 0.975))
        lst[[param]]$mean <- mean(x)
        lst[[param]]$samples <- x
      }
      summary$theta[[col]] <- lst
      names(summary$theta)[col] <- names(theta[[1]])[col]
    }
  fullsummary[[h]] <- summary
  }
  return(fullsummary)
}

###summarise mcmc results into posterior mean and quantiles
summarisenmidas <- function(mcmc.results, burnin){
  fullsummary <- list()
  
  for(h in 1:length(mcmc.results)){
    #summarise phi
    summary <- list()
    phi <- mcmc.results[[h]][[1]]
    phi <- phi[-c(1:burnin)]
    for(col in 1:length(phi[[1]])){
      lst <- list()
      x <- unlist(lapply(phi, "[[", col))
      lst$quantile <- quantile(x, probs=c(0.025, 0.5, 0.975))
      lst$mean <- mean(x)
      lst$samples <- x
      summary$phi[[col]] <- lst
      names(summary$phi)[col] <- names(phi[[1]])[col]
    }
    
    #summarise sigma
    lst <- list()
    sigma <- mcmc.results[[h]][[2]][-c(1:burnin)]
    x <- unlist(sigma)
    lst$quantile <- quantile(x, probs=c(0.025, 0.5, 0.975))
    lst$mean <- mean(x)
    lst$samples <- x
    summary$sigma <- lst
    
    fullsummary[[h]] <- summary
  }
    
  return(fullsummary)
}

#use mcmc summary to get h-step ahead point forecast
predictMCMC <- function(summary, z.high, z.low, shape){
  preds <- list()
    
  for(h in 1:length(summary)){
    phi <- lapply(summary[[h]]$phi, "[[", "mean")
    phi <- as.matrix(unlist(phi))

    theta <- summary[[h]]$theta
    for(i in 1:length(theta)){
      theta[[i]] <- unlist(lapply(theta[[i]], "[[", "mean"))
    }

    z <- compress(theta, z.high, z.low, shape)
    if(nrow(phi) == 16){
       z <- z[,-1]
    }
    predy <- z %*% phi
    preds[[h]] <- predy
  }
  return(preds)
}

predictnmidas <- function(summary, z){
  preds <- list()
  for(h in 1:length(summary)){
    phi <- sapply(summary[[h]]$phi, "[[", "mean")
    phi <- as.matrix(phi)
    
    predy <- z %*% phi
    preds[[h]] <- predy
  }
  return(preds)
}

#use mcmc summary to get h-step ahead density forecast
predictdensMCMC <- function(summary, z.high, z.low, shape){
  preds <- list()
  iter <- length(summary[[1]][["phi"]][[1]][["samples"]])
  
  for(h in 1:length(summary)){
    prediter <- matrix(nrow=nrow(z.low), ncol=4000)
    for(i in 1:iter){
      phi <- lapply(lapply(summary[[h]]$phi, "[[", "samples"), "[[", i)
      phi <- as.matrix(unlist(phi))
      theta <- summary[[h]]$theta
      for(k in 1:length(theta)){
        theta[[k]] <- unlist(lapply(lapply(theta[[k]], "[[", "samples"), "[", k))
      }
      z <- compress(theta, z.high, z.low, shape)
      if(nrow(phi) == 16){
        z <- z[,-1]
      }
      predy <- z %*% phi
      prediter[,i] <- predy
    }
    
    preds[[h]] <- prediter #each row is 1 timepoint, each column is the density

  }
  return(preds)
}

predictdensnmidas <- function(summary, z){
  preds <- list()
  iter <- length(summary[[1]][["phi"]][[1]][["samples"]])
  
  for(h in 1:length(summary)){
    prediter <- matrix(nrow=nrow(z), ncol=iter)
    for(i in 1:iter){
      phi <- lapply(lapply(summary[[h]]$phi, "[[", "samples"), "[[", i)
      phi <- as.matrix(unlist(phi))
      predy <- z %*% phi
      prediter[,i] <- t(predy)
    }
    
    preds[[h]] <- prediter #each row is 1 timepoint, each column is the density
    
  }
  return(preds)
}

#use mcmc summary to get h-step ahead point forecast for rolling CV
predicttestMCMC <- function(summary, z.high, z.low, shape){
  preds <- list()
  
  for(h in 1:length(summary)){
    phi <- lapply(summary[[h]]$phi, "[[", "mean")
    phi <- as.matrix(unlist(phi))
    
    theta <- summary[[h]]$theta
    for(i in 1:length(theta)){
      theta[[i]] <- unlist(lapply(theta[[i]], "[[", "mean"))
    }
    
    z <- compress(theta, z.high, z.low, shape)
    if(length(phi)==16){
      z <- z[,2:17]
    }
    predy <- z %*% phi
    
    ##density
    prediter <- c()
    iter <- length(summary[[1]][["phi"]][[1]][["samples"]])
    
    for(i in 1:iter){
      phi <- lapply(lapply(summary[[h]]$phi, "[[", "samples"), "[[", i)
      phi <- as.matrix(unlist(phi))
      theta <- summary[[h]]$theta
      
      for(k in 1:length(theta)){
        theta[[k]] <- unlist(lapply(lapply(theta[[k]], "[[", "samples"), "[", i))
      }
        
      z <- compress(theta, z.high, z.low, shape)
      if(length(phi)==16){
        z <- z[,2:17]
      }
      prediter <- c(prediter, z %*% phi)

    }
    
    preds[[h]] <- prediter #each row is 1 timepoint, each column is the density
    
    
    preds[[h]] <- list(predy, prediter)
  }
  return(preds)
}
