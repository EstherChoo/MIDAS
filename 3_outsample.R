library(tidyverse)
library(doParallel)
library(foreach)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")
sapply(paste0("./functions/", list.files("./functions/")), source)
load("data.Rdata")

rollingcv <- function(r){
  rowstart <- round(0.7 * nrow(diseases)) #to put for loop
  rowseq <- seq(rowstart, nrow(diseases), 1)
  row <- rowseq[r]
  
  for(d in colnames(diseases)){
    #split the data
    y.train <- as.matrix(diseases[1:row-1,d])
    y.test <- as.matrix(diseases[row:nrow(dis_lag),d])
    
    z.low.train <- dplyr::select(dis_lag[1:row-1,], starts_with(d))
    z.low.test <- dplyr::select(dis_lag[row:nrow(dis_lag),], starts_with(d))
    
    z.high.train <- list()
    z.high.test <- list()
    
    for(i in 1:length(climate_high_lst)){
      z.high.train[[i]] <- climate_high_lst[[i]][1:row-1,]
      z.high.test[[i]] <- climate_high_lst[[i]][row:nrow(dis_lag),]
    }
    
    climate_low.train <- climate_low[1:row-1,]
    climate_low.test <- climate_low[row:nrow(dis_lag),]
    
    z.agg.train <- cbind(climate_low.train, z.low.train)
    z.agg.test <- cbind(climate_low.test, z.low.test)
    
    ####Out-of-sample forecast####
    train.mcmcar1 <- list()
    train.mcmcar2 <- list()
    train.mcmcar3 <- list()

    train.mcmcbl1 <- list()
    train.mcmcbl2 <- list()
    train.mcmcbl3 <- list()

    train.mcmcbgl1 <- list()
    train.mcmcbgl2 <- list()
    train.mcmcbgl3 <- list()

    train.mcmcbal1 <- list()
    train.mcmcbal2 <- list()
    train.mcmcbal3 <- list()

     train.nmidasagg <- list()

    #train each step model on respective train set
    for(h in 1:12){
      train.mcmcar1[[h]] <- mcmc.ar(h, y.train, z.low.train, z.high.train, iterations=5000, shape=1)
      train.mcmcar2[[h]] <- mcmc.ar(h, y.train, z.low.train, z.high.train, iterations=5000, shape=2)
      train.mcmcar3[[h]] <- mcmc.ar(h, y.train, z.low.train, z.high.train, iterations=5000, shape=3)
    }

    cat("timepoint:", row, "progress: midas trained")
    cat("\n")

    for(h in 1:12){
      train.mcmcbl1[[h]] <- mcmc.bl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=1)
      train.mcmcbl2[[h]] <- mcmc.bl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=2)
      train.mcmcbl3[[h]] <- mcmc.bl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=3)
    }
    cat("timepoint:", row, "progress: midas bl trained")
    cat("\n")

    for(h in 1:12){
      train.mcmcbal1[[h]] <- mcmc.bal(h, y.train, z.low.train, z.high.train, iterations=5000, shape=1)
      train.mcmcbal2[[h]] <- mcmc.bal(h, y.train, z.low.train, z.high.train, iterations=5000, shape=2)
      train.mcmcbal3[[h]] <- mcmc.bal(h, y.train, z.low.train, z.high.train, iterations=5000, shape=3)
    }
    cat("timepoint:", row, "progress: midas bal trained")
    cat("\n")

    for(h in 1:12){
      train.mcmcbgl1[[h]] <- mcmc.bgl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=1)
      train.mcmcbgl2[[h]] <- mcmc.bgl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=2)
      train.mcmcbgl3[[h]] <- mcmc.bgl(h, y.train, z.low.train, z.high.train, iterations=5000, shape=3)
    }
    cat("timepoint:", row, "progress: midas bgl trained")
    cat("\n")

    #summarise the results of trained model
    out.summary.ar1 <- summariseMCMC(train.mcmcar1, 1000)
    out.summary.ar2 <- summariseMCMC(train.mcmcar2, 1000)
    out.summary.ar3 <- summariseMCMC(train.mcmcar3, 1000)

    out.summary.bl1 <- summariseMCMC(train.mcmcbl1, 1000)
    out.summary.bl2 <- summariseMCMC(train.mcmcbl2, 1000)
    out.summary.bl3 <- summariseMCMC(train.mcmcbl3, 1000)

    out.summary.bal3 <- summariseMCMC(train.mcmcbal3, 1000)
    out.summary.bal1 <- summariseMCMC(train.mcmcbal1, 1000)
    out.summary.bal2 <- summariseMCMC(train.mcmcbal2, 1000)

    out.summary.bgl1 <- summariseMCMC(train.mcmcbgl1, 1000)
    out.summary.bgl2 <- summariseMCMC(train.mcmcbgl2, 1000)
    out.summary.bgl3 <- summariseMCMC(train.mcmcbgl3, 1000)
    cat("timepoint", row, "summarised")
    cat("\n")

    #save each timepoint
    pred.ar1 <- predicttestMCMC(out.summary.ar1, z.high.test, z.low.test, 1)
    pred.ar2 <- predicttestMCMC(out.summary.ar2, z.high.test, z.low.test, 2)
    pred.ar3 <- predicttestMCMC(out.summary.ar3, z.high.test, z.low.test, 3)

    pred.bl1 <- predicttestMCMC(out.summary.bl1, z.high.test, z.low.test, 1)
    pred.bl2 <- predicttestMCMC(out.summary.bl2, z.high.test, z.low.test, 2)
    pred.bl3 <- predicttestMCMC(out.summary.bl3, z.high.test, z.low.test, 3)

    pred.bal1 <- predicttestMCMC(out.summary.bal1, z.high.test, z.low.test, 1)
    pred.bal2 <- predicttestMCMC(out.summary.bal2, z.high.test, z.low.test, 2)
    pred.bal3 <- predicttestMCMC(out.summary.bal3, z.high.test, z.low.test, 3)

    pred.bgl1 <- predicttestMCMC(out.summary.bgl1, z.high.test, z.low.test, 1)
    pred.bgl2 <- predicttestMCMC(out.summary.bgl2, z.high.test, z.low.test, 2)
    pred.bgl3 <- predicttestMCMC(out.summary.bgl3, z.high.test, z.low.test, 3)

    savepred.ar <- list(pred.ar1, pred.ar2, pred.ar3)
    savepred.bl <- list(pred.bl1, pred.bl2, pred.bl3)
    savepred.bal <-list(pred.bal1, pred.bal2, pred.bal3)
    savepred.bgl <- list(pred.bgl1, pred.bgl2, pred.bgl3)

    names(savepred.ar) <- c("ar1", "ar2", "ar3")
    names(savepred.bl) <- c("bl1", "bl2", "bl3")
    names(savepred.bal) <- c("bal1", "bal2", "bal3")
    names(savepred.bgl) <- c("bgl1", "bgl2", "bgl3")
    
    for(h in 1:12){
      train.nmidasagg[[h]] <- nomidas.ar(h, y.train, as.matrix(z.agg.train), iterations=5000)
    }
    cat("nmidasagg trained")
    cat("\n")

    summary.nmidasagg <- summarisenmidas(train.nmidasagg, 1000)
    
    pred.nmidasagg <- predictnmidas(summary.nmidasagg, as.matrix(z.agg.test))
    
    dens.nmidasagg <- predictdensnmidas(summary.nmidasagg, as.matrix(z.agg.test))
    
    nmidasagg <- list()
    for(h in 1:12){
      nmidasagg[[h]] <- list(pred.nmidasagg[[h]], dens.nmidasagg[[h]])
    }
    
    #save each model into rdata file
    assign(paste0("arpred", row), savepred.ar)
    paste0("arpred", row) %>%
      save(list=., file = paste0("./rollingcv/",d, "_ar", row, ".Rdata" ))
    rm(list=ls(pattern="(arpred)"))
    
    assign(paste0("blpred", row), savepred.bl)
    paste0("blpred", row) %>%
      save(list=., file = paste0("./rollingcv/", d, "_bl", row, ".Rdata" ))
    rm(list=ls(pattern="(blpred)"))
    
    assign(paste0("balpred", row), savepred.bal)
    paste0("balpred", row) %>%
      save(list=., file = paste0("./rollingcv/",d, "_bal", row, ".Rdata" ))
    rm(list=ls(pattern="(balpred)"))
    
    assign(paste0("bglpred", row), savepred.bgl)
    paste0("bglpred", row) %>%
      save(list=., file = paste0("./rollingcv/",d, "_bgl", row, ".Rdata" ))
    rm(list=ls(pattern="(bglpred)"))
    
    save(nmidasagg, file=paste0("./rollingcv/", d, "_nmidasagg", row, ".Rdata"))
    
    cat("timepoint:", row, "completed")
    cat("/n")
  }

}

cl <- makeCluster(21)
registerDoParallel(cl)

foreach(r=1:21, .packages=c('invgamma', 'statmod', 'MASS', 'mvtnorm', 'tidyverse')) %dopar% {
  sink("rollingcv_log.txt", append=TRUE) #saves the printed stuff into log.txt
  rollingcv(r=r) #rollingcv iterates across "diseases" specified in foreach function
}
  

#assess
for(i in 1:12){
  ypred <- pred.ar1[[i]][[1]]
  yactual <- y.test
  mape <- mean(abs((yactual-ypred)/yactual))
  print(mape)
}

