library(tidyverse)
library(doParallel)
library(foreach)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)
library(forecast)
library(EnvStats)
library(scoringRules)
library(ggplot2)
library(cowplot)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")
sapply(paste0("./functions/", list.files("./functions/")), source)
load("data.Rdata")

###run out of sample for pre-covid and during covid period

sa <- function(){ #truncate to week 52, 2019
  for(d in colnames(diseases)){
    if(d!="Cpox") next
    
    y <- as.matrix(diseases[1:406, d])
    print(d)
    z.low <- dplyr::select(dis_lag[1:406,], starts_with(d))
    
    for(i in 1:4){
      climate_high_lst[[i]] <- climate_high_lst[[i]][1:406,]
    }
    z.high <- climate_high_lst
    climate_low <- climate_low[1:406,]
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
         file=paste0("./out/sa_", d, "insampleresults.Rdata"))
    
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
         file=paste0("./out/sa_", d, "mcmcsummary.Rdata"))
  }
}


sa()

cl <- makeCluster(4)
registerDoParallel(cl)

foreach(d=colnames(diseases), .packages=c('invgamma', 'statmod', 'MASS', 'mvtnorm', 'tidyverse')) %dopar% {
  sapply(paste0("./functions/", list.files("./functions/")), source)
  load("data.Rdata")
  diseases<- diseases[-2]
  sa_cov(d) #rollingcv iterates across "diseases" specified in foreach function
}

sa_cov <- function(d){ #truncate to week 52, 2019
  
    y <- as.matrix(diseases[407:666, d])
    print(d)
    z.low <- dplyr::select(dis_lag[407:666,], starts_with(d))
    for(i in 1:4){
      climate_high_lst[[i]] <- climate_high_lst[[i]][407:666,]
    }
    z.high <- climate_high_lst
    climate_low <- climate_low[407:666,]
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
         file=paste0("./out/sacov_", d, "insampleresults.Rdata"))
    
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
         file=paste0("./out/sacov_", d, "mcmcsummary.Rdata"))
}


load("data.Rdata")

coef.df <- function(model){
  lst <- model[[1]][["phi"]]
  if(length(lst)==17){lst <- lst[2:17]}
  if(length(model[[1]])==2){
    names(lst) <- c("Cases Lag 1", "Cases Lag 2", "Cases Lag 3", "Cases Lag 4", "Cases Lag 5",
                    "Cases Lag 6", "Cases Lag 7", "Cases Lag 8", "Cases Lag 9", "Cases Lag 10", "Cases Lag 11",
                    "Cases Lag 12", "Total Precipitation (mm)", "Mean Temp (\u00B0C)", "Min Temp (\u00B0C)", "Max Temp (\u00B0C)")
  } else {
    names(lst) <- c("Cases Lag 1", "Cases Lag 2", "Cases Lag 3", "Cases Lag 4", "Cases Lag 5",
                    "Cases Lag 6", "Cases Lag 7", "Cases Lag 8", "Cases Lag 9", "Cases Lag 10", "Cases Lag 11",
                    "Cases Lag 12", "Total Precipitation (mm)", "Max Temp (\u00B0C)", "Min Temp (\u00B0C)", "Mean Temp (\u00B0C)")
  }
  
  conf <- lapply(lst, "[[", "quantile")
  df <- data.frame(conf.low = sapply(conf, "[[", "2.5%"),
                   conf.high = sapply(conf, "[[", "97.5%"),
                   estimate = sapply(lst, "[[", "mean"),
                   term = names(lst))
  return(df)
}

modelnames <- c("AR-(N)", "AR-(D)", "AR-(H)", "BAL-(N)", "BAL-(D)", "BAL-(H)", "BGL-(N)", "BGL-(D)", "BGL-(H)",
                "BL-(N)", "BL-(D)", "BL-(H)", "NMIDAS")
dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox", "HFMD"="HFMD")


coefplot <- list()
coefplotc <- list()

for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/sa_", d, "mcmcsummary.Rdata"))
  
  full_clim <- data.frame()
  i <- 1
  for(model in ls(pattern="in.summary")){
    print(model)
    df <- coef.df(get(model))
    tp_ind <- which(df$term=="Total Precipitation (mm)")
    df[tp_ind, 1:3] <- df[tp_ind, 1:3]
    df <- cbind(df, "model" = modelnames[i])
    i <- i + 1
    
    full_clim <- rbind(full_clim, df[13:16,])
  }
  
  full_clim <- full_clim %>%
    mutate(model = factor(model, levels=unique(model))) %>%
    mutate(term = factor(term, levels=unique(term)))
  
  load(paste0("./out/sacov_", d, "mcmcsummary.Rdata"))
  
  full_clim2 <- data.frame()
  i <- 1
  for(model in ls(pattern="in.summary")){
    print(model)
    df <- coef.df(get(model))
    tp_ind <- which(df$term=="Total Precipitation (mm)")
    df[tp_ind, 1:3] <- df[tp_ind, 1:3]
    df <- cbind(df, "model" = modelnames[i])
    i <- i + 1
    
    full_clim2 <- rbind(full_clim2, df[13:16,])
  }
  
  full_clim2 <- full_clim2 %>%
    mutate(model = factor(model, levels=unique(model))) %>%
    mutate(term = factor(term, levels=unique(term)))
  
  # full_clim2$covid <- "post"
  # full_clim$covid <- "pre"
  # 
  # full_clim <- rbind(full_clim, full_clim2)  
  
  colours <- c("royalblue1", "royalblue2", "royalblue3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "seagreen2",
               "seagreen3", "seagreen4", "orange1", "orange2", "orange3", "firebrick")
  
  ymin <- c()
  ymax <- c()
  full <- full_clim
  n <- 4
  
  for(i in seq(1, n, 2)){
    if(i==1){
      ymin <- c(ymin, as.numeric(full$term[[i]])-0.6)
    } else {
      ymin <- c(ymin, as.numeric(full$term[[i]])-0.5)
    }
    
    if(i==n){
      ymax <- c(ymax, as.numeric(full$term[[i]])+0.6)
    } else {
      ymax <- c(ymax, as.numeric(full$term[[i]])+0.5)
    }
  }
  
  stripe_midas <- data.frame(ymin=ymin, ymax=ymax, xmin=-Inf, xmax=Inf)
  
  coefplot[[d]] <- local({ggplot(full_clim) +
      aes(y=term, x=estimate, color=model, linetype=model, group=model) + 
      geom_rect(xmin=stripe_midas$xmin[1], ymin=stripe_midas$ymin[1], xmax=stripe_midas$xmax[1], ymax=stripe_midas$ymax[1], fill="grey90", inherit.aes=F)+
      geom_rect(xmin=stripe_midas$xmin[2], ymin=stripe_midas$ymin[2], xmax=stripe_midas$xmax[2], ymax=stripe_midas$ymax[2], fill="grey90", inherit.aes=F)+
      geom_point(aes(fill="Posterior Mean"), position = position_dodge(width=.75)) + 
      geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),position=position_dodge(width=.75), height=0) + 
      labs(x="Change in 1-week-ahead cases", y="", colour="Models", title=dict[d], fill=NULL, linetype="Models") + 
      geom_vline(xintercept=0, linetype="dashed") +
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +      scale_color_manual(name="Models", values=colours) +
      scale_linetype_manual(name="Models", values=c(rep("dotdash",7), "solid", rep("dotdash",4), "dashed")) +
      scale_fill_manual(values=NA)})
  
  coefplotc[[d]] <- local({ggplot(full_clim2) +
      aes(y=term, x=estimate, color=model, linetype=model, group=model)
      geom_rect(xmin=stripe_midas$xmin[1], ymin=stripe_midas$ymin[1], xmax=stripe_midas$xmax[1], ymax=stripe_midas$ymax[1], fill="grey90", inherit.aes=F)+
      geom_rect(xmin=stripe_midas$xmin[2], ymin=stripe_midas$ymin[2], xmax=stripe_midas$xmax[2], ymax=stripe_midas$ymax[2], fill="grey90", inherit.aes=F)+
      geom_point(aes(fill="Posterior Mean"), position = position_dodge(width=.75)) + 
      geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),position=position_dodge(width=.75), height=0) + 
      labs(x="Change in 1-week-ahead cases", y="", colour="Models", title=dict[d], fill=NULL, linetype="Models") + 
      geom_vline(xintercept=0, linetype="dashed") +
      theme_bw() +
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +      scale_color_manual(name="Models", values=colours) +
      scale_linetype_manual(name="Models", values=c(rep("dotdash",7), "solid", rep("dotdash",4), "dashed")) +
      scale_fill_manual(values=NA)})
  
}

p1 <- ggpubr::ggarrange(plotlist=coefplot, common.legend=T, nrow=4, ncol=1, legend="none")
p2 <- ggpubr::ggarrange(plotlist=coefplotc, common.legend=T, nrow=4, ncol=1, legend="bottom")   
finalplot <- ggarrange(p1, p2)
load("data.Rdata")
sapply(paste0("./functions/", list.files("./functions/")), source)

bestmodel <- c()
for(d in colnames(diseases)){
  load(paste0("./out/", d, "_waic.Rdata"))
  waic <- lapply(waic.list, "[[", 1)
  bestmodel[d] <- names(which.min(waic))
}

#panel plot
#df: dis, day, weight, weightlow, weighthigh, model, variable

var <- c("Total Precipitation", "Mean Temperature", "Min Temperature", "Max Temperature")

coefdf <- data.frame(matrix(ncol=8, nrow=0))
colnames(coefdf) <- c("Disease", "Variable", "Model", "Weight", "Weightlow", "Weighthigh", "Day", "Max")

dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox")

for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/sa_", d, "mcmcsummary.Rdata"))
  
  the <- list()
  the.low <- list()
  the.high <- list()
  
  summary <- get(paste0("in.summary.", tolower(bestmodel[d])))
  
  for(i in 1:4){
    the[[i]] <- sapply(summary[[1]][["theta"]][[i]], "[[", "mean")
    
    the.low[[i]] <- sapply(lapply(summary[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "2.5%")
    the.high[[i]] <- sapply(lapply(summary[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "97.5%")
  }
  
  for(i in 1:4){
    shape <- as.numeric(substr(bestmodel[d], nchar(bestmodel[d]), nchar(bestmodel[d])))
    w <- midasWeights(the, shape)[[i]] * 100
    wl <- midasWeights(the.low, shape)[[i]] * 100
    wh <- midasWeights(the.high, shape)[[i]] * 100
    
    
    df <- data.frame(Disease=dict[d], Variable=var[i], Model=bestmodel[d], Weight=w, Weightlow=wl, Weighthigh=wh, Day=1:30, Max=F)
    df[which.max(df$Weight), "Max"] = T
    coefdf <- rbind(coefdf, df)
  }
  
}

df_max <- filter(coefdf, Max==T)
df_max_precov <- df_max
coefdf_precov <- coefdf

for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/sacov_", d, "mcmcsummary.Rdata"))
  
  the <- list()
  the.low <- list()
  the.high <- list()
  
  summary <- get(paste0("in.summary.", tolower(bestmodel[d])))
  
  for(i in 1:4){
    the[[i]] <- sapply(summary[[1]][["theta"]][[i]], "[[", "mean")
    
    the.low[[i]] <- sapply(lapply(summary[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "2.5%")
    the.high[[i]] <- sapply(lapply(summary[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "97.5%")
  }
  
  for(i in 1:4){
    shape <- as.numeric(substr(bestmodel[d], nchar(bestmodel[d]), nchar(bestmodel[d])))
    w <- midasWeights(the, shape)[[i]] * 100
    wl <- midasWeights(the.low, shape)[[i]] * 100
    wh <- midasWeights(the.high, shape)[[i]] * 100
    
    
    df <- data.frame(Disease=dict[d], Variable=var[i], Model=bestmodel[d], Weight=w, Weightlow=wl, Weighthigh=wh, Day=1:30, Max=F)
    df[which.max(df$Weight), "Max"] = T
    coefdf <- rbind(coefdf, df)
  }
  
}

df_max <- filter(coefdf, Max==T)

ggplot(coefdf_precov) +
  geom_line(aes(x = Day, y = Weight, colour = "Parameter Posterior Mean", linetype = "Parameter Posterior Mean")) +
  geom_line(aes(x = Day, y = Weightlow, colour = "2.5% Parameter LB", linetype = "2.5% Parameter LB")) +
  geom_line(aes(x = Day, y = Weighthigh, colour = "97.5% Parameter UB", linetype = "97.5% Parameter UB")) +
  geom_point(data = df_max, aes(x = Day, y = Weight)) +
  geom_label_repel(
    data = df_max_precov, 
    aes(x = Day, y = Weight, 
        label = paste0("Day ", Day, ", ", round(Weight, 1), "% (", round(Weightlow, 1), 
                       "%, ", round(Weighthigh, 1), "%)")), 
    nudge_y = 2
  ) +
  facet_grid(Disease ~ Variable, scales = "free") +
  scale_color_manual(
    values = c("2.5% Parameter LB" = "blue", "97.5% Parameter UB" = "red", "Parameter Posterior Mean" = "black")
  ) +
  scale_linetype_manual(values = c(2, 2, 1)) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +  # Set x-axis breaks from 0 to 30 in increments of 5
  labs(
    x = "Days", 
    y = "Contribution of days to change in 1-week ahead cases (%)", 
    colour = NULL, 
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14), # Increase legend title size
    axis.title = element_text(size = 16),   # Increase axis title size
    axis.text = element_text(size = 14)     # Increase axis tick labels size
  )

ggplot(coefdf) +
  geom_line(aes(x = Day, y = Weight, colour = "Parameter Posterior Mean", linetype = "Parameter Posterior Mean")) +
  geom_line(aes(x = Day, y = Weightlow, colour = "2.5% Parameter LB", linetype = "2.5% Parameter LB")) +
  geom_line(aes(x = Day, y = Weighthigh, colour = "97.5% Parameter UB", linetype = "97.5% Parameter UB")) +
  geom_point(data = df_max, aes(x = Day, y = Weight)) +
  geom_label_repel(
    data = df_max, 
    aes(x = Day, y = Weight, 
        label = paste0("Day ", Day, ", ", round(Weight, 1), "% (", round(Weightlow, 1), 
                       "%, ", round(Weighthigh, 1), "%)")), 
    nudge_y = 2
  ) +
  facet_grid(Disease ~ Variable, scales = "free") +
  scale_color_manual(
    values = c("2.5% Parameter LB" = "blue", "97.5% Parameter UB" = "red", "Parameter Posterior Mean" = "black")
  ) +
  scale_linetype_manual(values = c(2, 2, 1)) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +  # Set x-axis breaks from 0 to 30 in increments of 5
  labs(
    x = "Days", 
    y = "Contribution of days to change in 1-week ahead cases (%)", 
    colour = NULL, 
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    axis.line = element_line(color = 'black'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 14), # Increase legend title size
    axis.title = element_text(size = 16),   # Increase axis title size
    axis.text = element_text(size = 14)     # Increase axis tick labels size
  )




