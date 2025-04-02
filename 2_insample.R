library(tidyverse)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)

##set-up
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")
sapply(paste0("./functions/", list.files("./functions/")), source)
load(file="data.Rdata")

for(d in colnames(diseases)){
  y <- diseases[d]
  print(d)
  z.low <- dplyr::select(dis_lag, starts_with(d))
  z.high <- climate_high_lst
  z.agg <- as.matrix(cbind(z.low, climate_low))
  
  ##create empty lists that will hold model results
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
  
  ##train models and save lists
  for(h in 1:12){
    print(h)
    out.mcmcar1[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=1)
    out.mcmcar2[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=2)
    out.mcmcar3[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=3)
  }
  
  
  for(h in 1:12){
    print(h)
    out.mcmcbl1[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=1)
    out.mcmcbl2[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=2)
    out.mcmcbl3[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=3)
  }
  
  
  for(h in 1:12){
    print(h)
    out.mcmcbal1[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=1)
    out.mcmcbal2[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=2)
    out.mcmcbal3[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=3)
  }
  
  for(h in 1:12){
    print(h)
    out.mcmcbgl1[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=1)
    out.mcmcbgl2[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=2)
    out.mcmcbgl3[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=3)
  }
  
  for(h in 1:12){
    out.nmidasagg[[h]] <- nomidas.ar(h, y, z.agg, 5000)
  }
  
  save(list=c("out.mcmcar1", "out.mcmcar2", "out.mcmcar3", "out.mcmcbl1", "out.mcmcbl2", "out.mcmcbl3", 
              "out.mcmcbal1", "out.mcmcbal2", "out.mcmcbal3", "out.mcmcbgl1", "out.mcmcbgl2", "out.mcmcbgl3",
              "out.nmidasagg"),
       file=paste0("./out/2a_", d, "insampleresults.Rdata"))
  
  ##summarise results from mcmc
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
  
  #needed for traceplots
  save(list=grep("in.summary", ls(), value=T),
       file=paste0("./out/2b_", d, "mcmcsummary.Rdata"))
  
  ##use summary to get point forecasts
  in.pred.ar1 <- predictMCMC(in.summary.ar1, z.high, z.low, 1)
  in.pred.ar2 <- predictMCMC(in.summary.ar2, z.high, z.low, 2)
  in.pred.ar3 <- predictMCMC(in.summary.ar3, z.high, z.low, 3)
  
  in.pred.bl1 <- predictMCMC(in.summary.bl1, z.high, z.low, 1)
  in.pred.bl2 <- predictMCMC(in.summary.bl2, z.high, z.low, 2)
  in.pred.bl3 <- predictMCMC(in.summary.bl3, z.high, z.low, 3)
  
  in.pred.bal1 <- predictMCMC(in.summary.bal1, z.high, z.low, 1)
  in.pred.bal2 <- predictMCMC(in.summary.bal2, z.high, z.low, 2)
  in.pred.bal3 <- predictMCMC(in.summary.bal3, z.high, z.low, 3)
  
  in.pred.bgl1 <- predictMCMC(in.summary.bgl1, z.high, z.low, 1)
  in.pred.bgl2 <- predictMCMC(in.summary.bgl2, z.high, z.low, 2)
  in.pred.bgl3 <- predictMCMC(in.summary.bgl3, z.high, z.low, 3)
  
  in.pred.nmagg <- predictnmidas(in.summary.nmagg, z.agg)
  
  save(list=grep("in.pred", ls(), value=T),
       file=paste0("./out/2c_", d, "insamplepoint.Rdata"))
  
  ##use summary to get density forecasts
  in.dens.ar1 <- predictdensMCMC(in.summary.ar1, z.high, z.low, 1)
  in.dens.ar2 <- predictdensMCMC(in.summary.ar2, z.high, z.low, 2)
  in.dens.ar3 <- predictdensMCMC(in.summary.ar3, z.high, z.low, 3)
  
  in.dens.bl1 <- predictdensMCMC(in.summary.bl1, z.high, z.low, 1)
  in.dens.bl2 <- predictdensMCMC(in.summary.bl2, z.high, z.low, 2)
  in.dens.bl3 <- predictdensMCMC(in.summary.bl3, z.high, z.low, 3)
  
  in.dens.bal1 <- predictdensMCMC(in.summary.bal1, z.high, z.low, 1)
  in.dens.bal2 <- predictdensMCMC(in.summary.bal2, z.high, z.low, 2)
  in.dens.bal3 <- predictdensMCMC(in.summary.bal3, z.high, z.low, 3)
  
  in.dens.bgl1 <- predictdensMCMC(in.summary.bgl1, z.high, z.low, 1)
  in.dens.bgl2 <- predictdensMCMC(in.summary.bgl2, z.high, z.low, 2)
  in.dens.bgl3 <- predictdensMCMC(in.summary.bgl3, z.high, z.low, 3)
  
  in.dens.nmagg <- predictdensnmidas(in.summary.nmagg, z.agg)
  
  save(list=grep("in.dens", ls(), value=T),
       file=paste0("./out/2d_", d, "insampledensity.Rdata"))
}

###DIC/WAIC table
all_dic <- data.frame(matrix(nrow=0, ncol=15))
all_waic <- data.frame(matrix(nrow=0, ncol=15))
for(d in colnames(diseases)){
  load(paste0("./out/", d, "_dic.Rdata"))
  load(paste0("./out/", d, "_waic.Rdata"))
  waic <- as.data.frame(waic.list)
  waic$`Best Model` <- colnames(waic)[apply(waic, 1, which.min)]
  waic$Disease <- d
  dic.list <- lapply(dic.list, unlist)
  dic <- as.data.frame(dic.list)
  dic$`Best Model` <- colnames(dic)[apply(dic, 1, which.min)]
  dic$Disease <- d
  all_dic <- rbind(all_dic, dic)
  all_waic <- rbind(all_waic, waic)
}

table(all_waic$`Best Model`) #BAL2
table(all_dic$`Best Model`) #BGL2
write.csv(all_dic, file="./tables/dic.csv", row.names=F)
write.csv(all_waic, file="./tables/waic.csv", row.names=F)

