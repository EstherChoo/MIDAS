fulltrace <- function(mcmc.results){
  full <- list()
  len <- length(mcmc.results[[1]][[1]])
  #len <- 10000
  
  for(h in 1:length(mcmc.results)){
    phi <- mcmc.results[[h]][[2]]
    mat <- matrix(NA, nrow=len, ncol=length(phi[[1]]))
    for(col in 1:length(phi[[1]])){
      mat[,col] <- unlist(lapply(phi, "[[", col))
    }
    
    #summarise sigma
    mat <- cbind(mat, unlist(mcmc.results[[h]][[3]]))
    
    #summarise theta
    theta <- mcmc.results[[h]][[1]]
    for(col in 1:length(theta[[1]])){
      x12 <- lapply(theta, "[[", col)
      for(param in 1:2){
        mat <- cbind(mat, unlist(lapply(x12, "[", param))[1:len]) #can remove index later
      }
    }
    
    colnames(mat) <- c(names(mcmc.results[[1]][[2]][[1]]), "sigma", "ar1", "ar2", "rh1", "rh2",
                       "tp1", "tp2", "at1", "at2")
    full[[h]] <- mat
  }
  
  return(full)
}

fulltrace.nm <- function(mcmc.results){
  full <- list()
  len <- length(mcmc.results[[1]][[1]])
  
  for(h in 1:length(mcmc.results)){
    phi <- mcmc.results[[h]][[1]]
    mat <- matrix(NA, nrow=len, ncol=length(phi[[1]]))
    for(col in 1:length(phi[[1]])){
      mat[,col] <- unlist(lapply(phi, "[[", col))
    }
    
    #summarise sigma
    mat <- cbind(mat, unlist(mcmc.results[[h]][[2]]))
    
    colnames(mat) <- c(names(mcmc.results[[1]][[1]][[1]]), "sigma")
    full[[h]] <- mat
  }
  
  return(full)
}

# fulltrace.bgl <- function(mcmc.results){
#   full <- list()
#   
#   for(h in 1:length(mcmc.results)){
#     phi <- mcmc.results[[h]][[2]]
#     mat <- matrix(NA, nrow=5000, ncol=16) #hard-coded!
#     for(col in 1:length(phi[[1]])){
#       mat[,col] <- unlist(lapply(phi, "[[", col))
#     }
#     
#     #summarise sigma
#     mat <- cbind(mat, unlist(mcmc.results[[h]][[3]]))
#     
#     #summarise theta
#     theta <- mcmc.results[[h]][[1]]
#     for(col in 1:length(theta[[1]])){
#       x12 <- lapply(theta, "[[", col)
#       for(param in 1:2){
#         mat <- cbind(mat, unlist(lapply(x12, "[", param))[1:5000]) #can remove index later
#       }
#     }
#     
#     colnames(mat) <- c(names(mcmc.results[[1]][[2]][[1]]), "sigma", "ar1", "ar2", "rh1", "rh2",
#                        "tp1", "tp2", "at1", "at2")
#     full[[h]] <- mat
#   }
#   
#   return(full)
# }
