midaspoly <- function(theta, x, shape){ #x is a matrix of high freq data for 1 var
  #k is the index of the x
  totalcompressed <- c()
  
  num <- c()
  bigk <- ncol(x)
  
  # t1 <- theta[1] 
  # t2 <- theta[2]

  if(shape == 1){ #free weighting
    t1 <- theta[1] + 1
    t2 <- theta[2] + 1
  } else if(shape == 2){ #downward sloping
    t1 <- 1
    t2 <- theta[2] + 1
  } else { #hump shape
    t1 <- theta[1] + 1
    t2 <- theta[1] + theta[2] + 1
  }

  for(k in 1:bigk){
    xk <- (k-1)/(bigk-1)
    num[k] <- xk^(t1-1) * (1-xk)^(t2-1)
  }
    
  denom <- sum(num)
  coef <- num/denom
    
  for(i in 1:nrow(x)){
    row <- x[i,]
    compressed <- sum(coef*row)
    totalcompressed <- append(totalcompressed, compressed)
  }
  
  return(totalcompressed)
}

#compress high freq vars w midas poly and theta  
compress <- function(th, z.high, z.low, shape){
  z <- as.matrix(z.low)
  for(i in 1:length(z.high)){
    z <- cbind(z, midaspoly(th[[i]], z.high[[i]], shape))
  }
  colnames(z)[13:16] <- c("ah", "rh", "tp", "at")
  z <- cbind(1, z)  
}
  
midasWeights <- function(theta, shape){ #takes list of thetas 
  weights <- list()
  
  for(m in 1:length(climate_high_lst)){
    num <- c()
    bigk <- ncol(climate_high_lst[[m]])
    
    if(shape == 1){ #free weighting
      t1 <- theta[[m]][1] + 1 
      t2 <- theta[[m]][2] + 1
    } else if(shape == 2){ #downward sloping
      t1 <- 1
      t2 <- theta[[m]][2] + 1
    } else { #hump shape
      t1 <- theta[[m]][1] + 1
      t2 <- theta[[m]][1] + theta[[m]][2] + 1
    }
    
    for(k in 1:bigk){
      xk <- (k-1)/(bigk-1)
        num[k] <- xk^(t1-1) * (1-xk)^(t2-1)
    }
    
    denom <- sum(num)
    coef <- num/denom
    weights[[m]] <- coef
  }
  
  return(weights)
}




  