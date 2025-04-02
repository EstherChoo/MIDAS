library(coda)
library(fixest)
library(doParallel)
library(foreach)
library(stringr)

#setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")

load("data.Rdata")
#gewecke, traceplot 
#qqplot, mse, r2, adj-r2, fitted vs actual plot
#residual autocorrelation, dic, waic
#choose best model using dic/waic
#1-3: ar, 4-6: bal, 7-9: bgl, 10-12:bl, 13: nmidasagg, 14-16:nmidasar/bl/bal

sink("./logs/in_assess_log.txt", append=TRUE)
for(d in colnames(diseases)){
  
  sapply(paste0("./functions/", list.files("./functions/")), source)
  cat(d)
  
  insamplemeasures <- list()
  for(i in 1:13){
    insamplemeasures[[i]] <- list()
  }
  
  #1. Trace-plot
  load(paste0("./out/2a_", d, "insampleresults.Rdata"))
  ar1trace <- fulltrace(out.mcmcar1)
  ar2trace <- fulltrace(out.mcmcar2)
  ar3trace <- fulltrace(out.mcmcar3)
  
  bl1trace <- fulltrace(out.mcmcbl1)
  bl2trace <- fulltrace(out.mcmcbl2)
  bl3trace <- fulltrace(out.mcmcbl3)
  
  bal1trace <- fulltrace(out.mcmcbal1)
  bal2trace <- fulltrace(out.mcmcbal2)
  bal3trace <- fulltrace(out.mcmcbal3)
  
  bgl1trace <- fulltrace(out.mcmcbgl1)
  bgl2trace <- fulltrace(out.mcmcbgl2)
  bgl3trace <- fulltrace(out.mcmcbgl3)
  
  nmagg.trace <- fulltrace.nm(out.nmidasagg)
  
  save(list=ls(pattern="(trace)"), file=paste0("./out/4a_", d, "_traces.Rdata"))
  rm(list=ls(pattern="(trace)"))
  print("traceplot done")
  
  # 2. geweke diagnostics - equality of means btw first 10% and second half of chain
  # majority of points fall within 2sd
  # load("./out/4a_traces.Rdata")
  # geweke.plot(as.mcmc(bal1trace[[1]][1001:5000,])) #manually check plots
  
  
  #3. residual autocorrelation - independence of residuals -> all information are captured in the model structure
  load(file="./data.Rdata")
  load(file=paste0("./out/2c_", d, "insamplepoint.Rdata"))
  i <- 1
  y <- as.numeric(c(diseases[[d]]))
  for(pred in ls(pattern="(in.pred)")){
    lst <- list()
    lst.res <- list()
    pred <- get(pred)
    for(h in 1:12){
      n <- nrow(pred[[h]]) - h
      res <- y - c(pred[[h]])
      #res <- y[(h+1):nrow(y),] - pred[[h]][1:n]
      lst[[h]] <- acf(res, type="correlation", lag.max=12)
      lst.res[[h]] <- res
    }
    insamplemeasures[[i]][["res.acf"]] <- lst
    insamplemeasures[[i]][["residuals"]] <- lst.res
    i <- i + 1
  }
  
  
  #res autocorr tables
  modelnames <- c("AR1", "AR2", "AR3", "BL1", "BL2", "BL3", "BAL1", "BAL2", "BAL3",
                  "BGL1", "BGL2", "BGL3", "NMIDAS-AGG")
  
  racf <- lapply(insamplemeasures, "[[", "res.acf")
  resacf <- data.frame(`Lags`=rep(c(1:13), 12), `Forecast Horizon`=rep(c(1:12), each=13))
  acfval <- matrix(ncol=0, nrow=156)
  for(i in 1:13){
    acfval <- cbind(acfval, unlist(lapply(racf[[i]], "[[", "acf")))
  }
  resacf <- cbind(resacf, acfval)
  colnames(resacf)[3:15] <- modelnames
  
  write.csv(resacf, paste0("./tables/", d, "_resacf.csv"))
  print("residual done")
  
  #4. qqplot - normality of residuals - manual check
  #fat tails
  # model <- 10
  # horizon <- 1
  # r <- residuals[[model]][[horizon]]
  # sd.r <- sqrt(sum(y[(horizon+1):nrow(y)]-residuals[[model]][[horizon]])^2/(n-2))
  # qqnorm(r/sd.r)
  # qqline(r/sd.r, col="red")
  
  #5. r2, adj-r2, mse
  r2.mse <- function(pred, y){
    if(str_detect(pred, "ar|bal|bgl|bl")){
      p <- 16
    } else if (str_detect(pred, "nmagg")){
      p <- 36
    } else { p <- 12}
    r2list <- list()
    ar2list <- list()
    mselist <- list()
    pred <- get(pred)
    for(h in 1:12){
      #n <- nrow(pred[[h]]) - h
      mse <- mean((y - pred[[h]])^2)
      ssr <- sum((pred[[h]] - mean(y))^2)
      sst <- sum((y - mean(y))^2)
      r2 <- ssr/sst
      r2list[[h]] <- r2
      ar2list[[h]] <- 1 - ((1-r2)*(n-1))/(n-p-1)
      mselist[[h]] <- mse
    }
    return(list(r2list,ar2list, mselist))
  }
  
  
  i <- 1
  for(pred in ls(pattern="(in.pred)")){
    r <- r2.mse(pred, y)
    insamplemeasures[[i]][["r2"]] <- r[[1]]
    insamplemeasures[[i]][["adj-r2"]] <- r[[2]]
    insamplemeasures[[i]][["mse"]] <- r[[3]]
    i <- i + 1
  }
  
  print("r2, ar2, mse done")
  
  #5. AIC
  i <- 1
  for(pred in ls(pattern="(in.pred)")){
    if(str_detect(pred, "ar|bal|bgl|bl")){
      p <- 16
    } else if (str_detect(pred, "nmagg")){
      p <- 28
    } else { p <- 12}
    n <- nrow(get(pred)[[1]])
    rss <- c()
    aic <- list()
    for(h in 1:12){ 
      aic[[h]] <- n * log(sum(insamplemeasures[[i]][["residuals"]][[h]]^2) / n) + 2*p
    }
    insamplemeasures[[i]][["aic"]] <- aic
    i <- i + 1
  }
  
  #save as dataframe for easy copying into appendix
  models <- c("AR1", "AR2", "AR3", "BAL1", "BAL2", "BAL3", "BGL1", "BGL2", "BGL3", "BL1", "BL2", "BL3", 
              "NMIDASAGG")
  
  msedf <- data.frame(horizon=seq(1,12))
  r2df <- data.frame(horizon=seq(1,12))
  adjr2df <- data.frame(horizon=seq(1,12))
  aicdf <- data.frame(horizon=seq(1,12))
  
  for(i in  1:13){
    tempdf1 <- data.frame(mse=unlist(insamplemeasures[[i]][["mse"]]))
    tempdf2 <- data.frame(r2=unlist(insamplemeasures[[i]][["r2"]]))
    tempdf3 <- data.frame(adjr2=unlist(insamplemeasures[[i]][["adj-r2"]]))
    tempdf4 <- data.frame(aic=unlist(insamplemeasures[[i]][["aic"]]))
    colnames(tempdf1)[1] <- models[i]; colnames(tempdf2)[1] <- models[i]; colnames(tempdf3)[1] <- models[i]; colnames(tempdf4)[1] <- models[i]
    msedf <- cbind(msedf, tempdf1)
    r2df <- cbind(r2df, tempdf2)
    adjr2df <- cbind(adjr2df, tempdf3)
    aicdf <- cbind(aicdf, tempdf4)
  }
  
  write.csv(msedf, file=paste0("./tables/", d, "_mse.csv"), row.names=F)
  write.csv(r2df, file=paste0("./tables/", d, "_r2.csv"), row.names=F)
  write.csv(adjr2df, file=paste0("./tables/", d, "_adjr2.csv"), row.names=F)
  write.csv(aicdf, file=paste0("./tables/", d, "_aic.csv"), row.names=F)
  #save(insamplemeasures, file=paste0("./out/4a_", d, "_insamplemeasures.Rdata"))
  
  print("aic done")
  
  # #6. DIC
  # cl <- makeCluster(12, type="FORK")
  # registerDoParallel(cl)
  # 
  # load(paste0("./out/2d_", d, "insampledensity.Rdata"))
  # load(paste0("./out/2b_", d, "mcmcsummary.Rdata"))
  # predfiles <- ls(pattern="(in.dens)")
  # mcmcfiles <- ls(pattern="(in.summary)")
  # #bayes mean and sd is expectation of samples
  # llh <- function(lst){
  #   mu <- lst[[1]]
  #   sd <- lst[[2]] 
  #   return(sum(dnorm(dens, mu, sd, log=T)))}
  # 
  # llh.helper <- function(dens, sig){
  #   new <- list()
  #   for(i in 1:4000){
  #     new[[i]] <- list(dens[i], sqrt(sig[i]))
  #   }
  #   return(new)
  # }
  # 
  # dic.list <- list()
  # print("DIC running")
  # for(model in 1:13) {
  #   cat(model)
  #   cat("\n")
  #   pred <- get(predfiles[model])
  #   mcmc <- get(mcmcfiles[model])
  #   lst <- list()
  #   for(h in c(1:12)) {
  #     dic <- 0
  #     sig <- mcmc[[h]][["sigma"]][["samples"]]
  #     for(i in 1:666){
  #       dens <- pred[[h]][i,]
  #       samples <- llh.helper(dens, sig)
  #       l <- llh(list(mean(dens), mean(sqrt(sig))))
  #       pdic <- 2*(l - mean(sapply(samples, llh)))
  #       dic <- dic + -2*l + 2*pdic
  #     }
  #     lst[[h]] <- dic/666
  #   }
  #   dic.list[[model]] <- lst
  # }
  # 
  # stopCluster(cl)
  # names(dic.list) <- c("AR1", "AR2", "AR3", "BAL1", "BAL2", "BAL3", "BGL1", "BGL2", "BGL3", "BL1", "BL2", "BL3", 
  #                                "NMIDASAGG")
  # save(dic.list, file=paste0("./out/", d, "_dic.Rdata"))
  # 
  # dic.df <- as.data.frame(sapply(dic.list, c))
  # dic.df <- apply(dic.df,2,as.numeric)
  # write.csv(dic.df, paste0("./tables/", d, "_dic.csv"))
  # 
  # # model <- 12 #bgl3 - h8 h10
  # # pred <- get(predfiles[model])
  # # mcmc <- get(mcmcfiles[model])
  # # h <- 10
  # # dic <- 0
  # # sig <- mcmc[[h]][["sigma"]][["samples"]]
  # # for(i in 1:1047){
  # #   dens <- pred[[h]][i,]
  # #   samples <- llh.helper(dens, sig)
  # #   l <- llh(list(mean(dens), mean(sqrt(sig))))
  # #   pdic <- 2*(l - mean(sapply(samples, llh)))
  # #   dic <- dic + -2*l + 2*pdic
  # # }
  # # print(dic/1047)
  # 
  # 
  # #7. WAIC
  # load(paste0("./out/2d_", d, "insampledensity.Rdata"))
  # load(paste0("./out/2b_", d, "mcmcsummary.Rdata"))
  # predfiles <- ls(pattern="(in.dens)")
  # mcmcfiles <- ls(pattern="(in.summary)")
  # 
  # llh <- function(lst){
  #   mu <- lst[[1]]
  #   sd <- lst[[2]] 
  #   return(sum(dnorm(dens, mu, sd, log=T)))}
  # 
  # lkh <- function(lst){
  #   mu <- lst[[1]]
  #   sd <- lst[[2]] 
  #   return(sum(dnorm(dens, mu, sd, log=F)))}
  # 
  # llh.helper <- function(dens, sig){
  #   new <- list()
  #   for(i in 1:4000){
  #     new[[i]] <- list(dens[i], sqrt(sig[i]))
  #   }
  #   return(new)
  # }
  # 
  # cl <- makeCluster(12, type="FORK")
  # registerDoParallel(cl)
  # 
  # print("waic running")
  # waic.list <- list()
  # for(model in 1:13){
  #   cat(model)
  #   cat("\n")
  #   pred <- get(predfiles[model])
  #   mcmc <- get(mcmcfiles[model])
  #   waic.list[[model]] <- foreach(h=1:12) %dopar% {
  #     lppd <- 0
  #     pwaic <- 0
  #     sig <- mcmc[[h]][["sigma"]][["samples"]]
  #     for(i in 1:666){
  #       print(i)
  #       dens <- pred[[h]][i,]
  #       samples <- llh.helper(dens, sig)
  #       lppd <- log(mean(sapply(samples, lkh)))
  #       pwaic <- (lppd - mean(sapply(samples, llh)))
  #     }
  #     return(-2*lppd + 2*2*pwaic)
  #   }
  # }
  # 
  # waic.list <- lapply(waic.list, unlist)
  # 
  # stopCluster(cl)
  # names(waic.list) <- modelnames
  # save(waic.list, file=paste0("./out/", d, "_waic.Rdata"))
  # 
  # waic.df <- as.data.frame(sapply(waic.list, c))
  # waic.df <- apply(waic.df,2,as.numeric)
  # write.csv(waic.df, paste0("./tables/", d, "_waic.csv"))
}
sink()
 