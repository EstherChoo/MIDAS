library(tidyverse)
library(doParallel)
library(foreach)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)
library(EnvStats)

#folder <- "C:/Users/esthe/OneDrive - Nanyang Technological University/FYP/MyCode"
folder <- "/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/Code"
setwd(folder)
sapply(paste0("./functions/", list.files("./functions")), source)

###Data generating process
load(file="data.Rdata")
nlags <- 30

#coefficients - true value
save_mae <- list()
save_params <- list()

cl <- makeCluster(30, type="FORK")
registerDoParallel(cl)

foreach(nsim=1:100, .packages=c('invgamma', 'statmod', 'MASS', 'mvtnorm', 'tidyverse')) %dopar% {
  cat(nsim)
  #draw theta and phi from gamma and normal dist respectively
  real_theta1 <- c(rgamma(1, 2, 3), rgamma(1, 3, 2.5)) #rain
  real_theta2 <- c(rgamma(1, 1, 1.5), rgamma(1, 4.2, 3.4)) #meanT
  real_theta3 <- c(rgamma(1, 3.6, 0.7), rgamma(1, 3.2, 2.3)) #minT
  real_theta4 <- c(rgamma(1, 1.3, 1.6), rgamma(1, 4, 1)) #maxT
  real_theta <- list(real_theta1, real_theta2, real_theta3, real_theta4)
  
  real_phi <- c(rnorm(1, 3, 0.5), rnorm(1, 7.8, 2), rnorm(1, 4.8, 1.5), rnorm(1, 5.8, 2.8))
  real_weights <- midasWeights(real_theta, 1)
  
  for(i in 1:4){
    x <- real_weights[[i]] * real_phi[i]
    assign(paste0("coef", i), x)
    print(sum(get(paste0("coef", i))))
  }
  
  sim_cases <- coef1 %*% t(climate_high_lst[[1]]) + coef2 %*% t(climate_high_lst[[2]]) + 
    coef3 %*% t(climate_high_lst[[3]]) + coef4 %*% t(climate_high_lst[[4]])
  sim_cases <- t(sim_cases)
  #plot(sim_cases)
  
  
  ###sanity check - glm is able to retrieve coefficients###
  # nmidas <- lm(sim_cases ~ z.agg[,13:36] + sim_lag)
  # tp <- 0
  # for(i in 2:7){
  #   tp <- tp + nmidas[["coefficients"]][[i]]
  # }
  # tp <- tp/6
  
  # sum(coef1_d) * 6
  # testcoef <- coef1_d
  # testz <- rnorm(167, 30, 20)
  # 
  # sum(testcoef) + testcoef %*% testz
  # testcoef %*% (testz + 1)
  ###
  
  #train MIDAS and non-MIDAS models on simulated data
  out.mcmcar1 <- list()
  out.mcmcar2 <- list()
  out.mcmcar3 <- list()
  
  out.mcmcbl1 <- list()
  out.mcmcbl2 <- list()
  out.mcmcbl3 <- list()
  
  out.mcmcbal1 <- list()
  out.mcmcbal2 <- list()
  out.mcmcbal3 <- list()
  
  out.mcmcbgl1 <- list()
  out.mcmcbgl2 <- list()
  out.mcmcbgl3 <- list()
  
  out.nmidasagg <- list()
  
  h<-1
  z.high <- climate_high_lst
  sim_lag <- matrix(0, nrow=666, ncol=12)
  out.mcmcar1[[h]] <- mcmc.ar(h, sim_cases, sim_lag, z.high, iterations=5000, shape=1)
  out.mcmcar2[[h]] <- mcmc.ar(h, sim_cases, sim_lag, z.high, iterations=5000, shape=2)
  out.mcmcar3[[h]] <- mcmc.ar(h, sim_cases, sim_lag, z.high, iterations=5000, shape=3)
  
  out.mcmcbl1[[h]] <- mcmc.bl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=1)
  out.mcmcbl2[[h]] <- mcmc.bl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=2)
  out.mcmcbl3[[h]] <- mcmc.bl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=3)
  
  out.mcmcbal1[[h]] <- mcmc.bal(h, sim_cases, sim_lag, z.high, iterations=5000, shape=1)
  out.mcmcbal2[[h]] <- mcmc.bal(h, sim_cases, sim_lag, z.high, iterations=5000, shape=2)
  out.mcmcbal3[[h]] <- mcmc.bal(h, sim_cases, sim_lag, z.high, iterations=5000, shape=3)
  
  out.mcmcbgl1[[h]] <- mcmc.bgl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=1)
  out.mcmcbgl2[[h]] <- mcmc.bgl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=2)
  out.mcmcbgl3[[h]] <- mcmc.bgl(h, sim_cases, sim_lag, z.high, iterations=5000, shape=3)
  
  #aggdata <- cbind(z.agg[1:1035,13:36], sim_lag)
  out.nmidasagg[[h]] <- nomidas.ar(h, sim_cases, as.matrix(climate_low), 5000)

  #summarise the results of trained model
  in.summary.ar1 <- summariseMCMC(out.mcmcar1, 1000)
  in.summary.ar2 <- summariseMCMC(out.mcmcar2, 1000)
  in.summary.ar3 <- summariseMCMC(out.mcmcar3, 1000)
  
  in.summary.bl1 <- summariseMCMC(out.mcmcbl1, 1000)
  in.summary.bl2 <- summariseMCMC(out.mcmcbl2, 1000)
  in.summary.bl3 <- summariseMCMC(out.mcmcbl3, 1000)
  
  in.summary.bal1 <- summariseMCMC(out.mcmcbal1, 1000)
  in.summary.bal2 <- summariseMCMC(out.mcmcbal2, 1000)
  in.summary.bal3 <- summariseMCMC(out.mcmcbal3, 1000)
  
  in.summary.bgl1 <- summariseMCMC(out.mcmcbgl1, 1000)
  in.summary.bgl2 <- summariseMCMC(out.mcmcbgl2, 1000)
  in.summary.bgl3 <- summariseMCMC(out.mcmcbgl3, 1000)
  
  in.summary.nmagg <- summarisenmidas(out.nmidasagg, 1000)
  
  rm("in.summary.nmar", "in.summary.nmbl", "in.summary.nmbal")
  
  #get coef
  sumfiles <- ls(pattern="in.summary")
  params <- list()
  for(sum in sumfiles){
    sum1 <- get(sum)
    params[[sum]] <- sapply(sum1[[1]][["phi"]], "[[", "mean") 
  }
  
  #evaluate accuracy of coefficients - calc accuracy of coefficients (MAE)
  #get weighted coef over time
  for(sum in sumfiles){
    params[[sum]][["the"]] <- list(NA)
    for(i in 1:4){
      sum1 <- get(sum)
      params[[sum]][["the"]][[i]] <- sapply(sum1[[1]][["theta"]][[i]], "[[", "mean")
    }
  }
  
  mae <- data.frame("Model"=NA, "Error"=NA)
  for(sum in sumfiles){
    w <- list()
    if(sum=="in.summary.nmagg"){
      w[["Drt"]] <- rep(mean(unlist(params[[sum]][1:4])) / nlags, nlags)
      w[["Mt"]] <- rep(mean(unlist(params[[sum]][5:8])) / nlags, nlags)
      w[["MinT"]] <- rep(mean(unlist(params[[sum]][9:12])) / nlags, nlags)
      w[["MaxT"]] <- rep(mean(unlist(params[[sum]][13:16])) / nlags, nlags)
      params[[sum]][["weight"]] <- w
    }else{
      the <- params[[sum]][["the"]]
      w[["Drt"]] <- midasWeights(the, 1)[[1]] * params[[sum]][["ah"]]
      w[["Mt"]] <- midasWeights(the, 1)[[2]] * params[[sum]][["rh"]]
      w[["MinT"]] <- midasWeights(the, 1)[[3]] * params[[sum]][["at"]]
      w[["MaxT"]] <- midasWeights(the, 1)[[4]] * params[[sum]][["tp"]]
      params[[sum]][["weight"]] <- w
    }
    mae <- rbind(mae, c(sum, sum(abs(w[["Drt"]] - coef1)) +
                          sum(abs(w[["Mt"]] - coef2)) +
                          sum(abs(w[["MinT"]] - coef3)) +
                          sum(abs(w[["MaxT"]] - coef4))))
  }
  mae$Error <- as.numeric(mae$Error)
  mae <- mae[-1,]
  
  save("mae", "params", file=paste0("./sim/sim_", nsim, ".Rdata"))

}

###compare distribution of KL test
#get samples of params
#t1=Drt=ah, t2=MeanT=rh, t3=MinT=tp, t4=MaxT=at
load(paste0("./sim/sim_", 1, ".Rdata"))
p <- c("t1_1", "t1_2", "t2_1", "t2_2", "t3_1", "t3_2", "t4_1", "t4_2", "p1", "p2", "p3", "p4")
models <- names(params)

param_df <- expand.grid(p, models)
colnames(param_df) <- c("Parameter", "Model")
for(i in 1:100){
  new <- rep(0, nrow(param_df))
  param_df <- cbind(param_df, new)
  colnames(param_df)[2+i] <- paste0("Sample_", i)
}

for(i in 1:100){
  load(paste0("./sim/sim_", i, ".Rdata"))
  n <- 2+i
  print(i)
  for(model in names(params)){
    print(model)
    if(model == "in.summary.nmagg"){
      param_df[which(param_df$Parameter=="p1" & param_df$Model==model), n] <- mean(unlist(params[[model]][1:4]))
      param_df[which(param_df$Parameter=="p2" & param_df$Model==model), n] <- mean(unlist(params[[model]][5:8]))
      param_df[which(param_df$Parameter=="p3" & param_df$Model==model), n] <- mean(unlist(params[[model]][9:12]))
      param_df[which(param_df$Parameter=="p4" & param_df$Model==model), n] <- mean(unlist(params[[model]][13:16]))
      
    } else {
      param_df[which(param_df$Parameter=="p1" & param_df$Model==model), n] <- params[[model]][["ah"]]
      param_df[which(param_df$Parameter=="p2" & param_df$Model==model), n] <- params[[model]][["rh"]]
      param_df[which(param_df$Parameter=="p3" & param_df$Model==model), n] <- params[[model]][["tp"]] 
      param_df[which(param_df$Parameter=="p4" & param_df$Model==model), n] <- params[[model]][["at"]] 
      
      param_df[which(param_df$Parameter=="t1_1" & param_df$Model==model), n] <- params[[model]][["the"]][[1]][1]
      param_df[which(param_df$Parameter=="t2_1" & param_df$Model==model), n] <- params[[model]][["the"]][[2]][1]
      param_df[which(param_df$Parameter=="t3_1" & param_df$Model==model), n] <- params[[model]][["the"]][[3]][1]
      param_df[which(param_df$Parameter=="t4_1" & param_df$Model==model), n] <- params[[model]][["the"]][[4]][1]
      
      param_df[which(param_df$Parameter=="t1_2" & param_df$Model==model), n] <- params[[model]][["the"]][[1]][2]
      param_df[which(param_df$Parameter=="t2_2" & param_df$Model==model), n] <- params[[model]][["the"]][[2]][2]
      param_df[which(param_df$Parameter=="t3_2" & param_df$Model==model), n] <- params[[model]][["the"]][[3]][2]
      param_df[which(param_df$Parameter=="t4_2" & param_df$Model==model), n] <- params[[model]][["the"]][[4]][2]
    }
  }
}

# sample <- unlist(param_df[which(param_df$Parameter=="p1" & param_df$Model=="in.summary.ar1"), 3:102])
x <- seq(1,10, 0.1)

real_dist <- function(param, x) { #output diff dist based on param_df$Model
  if(param=="t1_1")
    return(dgamma(x, 2, 3))
  if(param=="t3_1") 
    return(dgamma(x, 1, 1.5))
  if(param=="t2_1") 
    return(dgamma(x, 3.6, 0.7))
  if(param=="t4_1") 
    return(dgamma(x, 1.3, 1.6))
  
  if(param=="t1_2") 
    return(dgamma(x, 3, 2.5))
  if(param=="t3_2") 
    return(dgamma(x, 4.2, 3.4))
  if(param=="t2_2") 
    return(dgamma(x, 3.2, 2.3))
  if(param=="t4_2") 
    return(dgamma(x, 4, 1))
  
  if(param=="p1") 
    return(dnorm(x, 3, 0.5))
  if(param=="p2") 
    return(dnorm(x, 7.8, 2))
  if(param=="p3") 
    return(dnorm(x, 4.8, 1.5))
  if(param=="p4") 
    return(dnorm(x, 5.8, 2.8))
  
}


kl_divergence_empirical_vs_theoretical <- function(sample,param) {
  #q function is the empirical sample from sim, p is the actual theoretical distribution
  
  # Estimate the density q(x) from samples
  q_density <- function(val) demp(val, obs=sample)
  
  # Theoretical density p(x)
  p_density <- function(val) real_dist(param, val)
  
  f <- function(x) {p_density(x) * (log(p_density(x)) - log(q_density(x)+1e-10))}
  if(grepl("t", param)){
    x <- seq(0.1, 50, 0.1)
  }else{
    x <- seq(-10, 20, 0.1)
  }
  if(anyNA(f(x))){print(f(x))}
  kl_div <- sum(f(x))
  
  return(kl_div)
}

param_df$kl <- 0
for(n in 1:nrow(param_df)){
  param <- param_df[n, "Parameter"]
  s <- unlist(param_df[n, 3:102])
  param_df[n, "kl"] <- kl_divergence_empirical_vs_theoretical(s, param)
  print(param)
  print(param_df[n, "kl"])
  
}

t <- dplyr::select(param_df, Model, Parameter, kl)
t$Model <- toupper(unlist(lapply(strsplit(as.character(t$Model), "\\."), "[[", 3)))
t$Model <- gsub("1", "-(N)", t$Model)
t$Model <-gsub("2", "-(D)", t$Model)
t$Model <-gsub("3", "-(H)", t$Model)
t$Model <-gsub("NMAGG", "NMIDAS", t$Model)

colours <- c("royalblue1", "royalblue2", "royalblue3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "seagreen2",
             "seagreen3", "seagreen4", "orange1", "orange2", "orange3", "firebrick")
ggplot(t) +
  geom_col(aes(x=Model, y=kl, group=Parameter,fill=Model)) +
  facet_wrap(~Parameter, scales="free") +
  scale_fill_manual(name="Models", values=colours) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

write.csv(param_df, file="./tables/sim_results.csv", row.names=F)



