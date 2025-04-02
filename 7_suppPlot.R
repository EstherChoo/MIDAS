library(ggplot2)
library(tidyverse)
library(ggpubr)

######################################################
#################### coef plot #######################
######################################################
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

for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/2b_", d, "mcmcsummary.Rdata"))
  rm("in.summary.nmar", "in.summary.nmbal", "in.summary.nmbl")
  
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
  
}

p2 <- ggpubr::ggarrange(plotlist=coefplot, common.legend=T, nrow=2, ncol=2, legend="bottom") 

######################################################
#################### fit plot ########################
######################################################
#load data
load("data.Rdata")
raw <- read.csv("climatedisease.csv")
ind <- dplyr::select(raw, "Start", "End")
diseases <- cbind(ind[13:678,], diseases)
diseases <- diseases[-4]

dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox")

#get best model for each disease
bestmodel <- c()
for(d in colnames(diseases)[3:6]){
  load(paste0("./out/", d, "_waic.Rdata"))
  waic <- lapply(waic.list, "[[", 1)
  bestmodel[d] <- names(which.min(waic))
}

date <- ind[["Start"]]
date <- as.Date(date, tryFormats = c("%m/%d/%Y"))
yearbreaks <- date[match(unique(year(date)), year(date))]
monthbreaks <- as.Date(c())
for(year in unique(year(date))){
  subsetdates <- date[year(date)==year]
  monthbreaks <- c(monthbreaks, subsetdates[match(unique(month(subsetdates)), month(subsetdates))])
}

date <- as.Date(date[13:678])

# diseases$date <- date[13:678]
# p1_df <- pivot_longer(diseases, names_to="Disease", values_to="Cases", cols=3:6) #hard coded
# p1_df$Disease <- dict[p1_df$Disease]

allplots <- list()
for(d in colnames(diseases)[3:6]){
  load(paste0("./out/2c_", d, "insamplepoint.Rdata"))
  pred <- get(paste0("in.pred.", tolower(bestmodel[d])))
  diseases[d] <- as.numeric(diseases[[d]])
  
  plot1.fitted1 <- ggplot() +
    geom_point(aes(x=date, y=diseases[[d]], colour="Actual")) + 
    geom_line(aes(x=date, y=pred[[1]][1:666], colour="Fitted")) +
    scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
    ylab("Case Counts") +
    xlab("Year") +
    labs(title=paste0("Fitted ", dict[d], " Case Counts (1 Week Ahead)")) +
    theme(legend.position = "top") +
    theme_bw() + 
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          axis.ticks.length = unit(5, "pt")) +
    theme(axis.text.x=element_text(hjust=-0.7)) +
    theme(
      legend.text = element_text(size = 14),   # Increase legend text size
      legend.title = element_text(size = 16),  # Increase legend title size
      legend.key.size = unit(1.5, "lines")     # Increase legend key size (symbols)
    ) +
    scale_x_date(minor_breaks=monthbreaks, breaks=yearbreaks,labels=sprintf("%02d", seq(12, 24, 1)) )
  
  plot1.fitted2 <- ggplot() +
    geom_point(aes(x=date[6:666], y=diseases[[d]][6:666], colour="Actual")) + 
    geom_line(aes(x=date[6:666], y=pred[[5]][1:661], colour="Fitted")) +
    scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
    ylab("Case Counts") +
    xlab("Year") +
    labs(title=paste0("Fitted ", dict[d], " Case Counts (5 Weeks Ahead)")) +
    theme(legend.position = "top") +
    theme_bw() + 
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          axis.ticks.length = unit(5, "pt")) +
    theme(axis.text.x=element_text(hjust=-0.7)) +
    theme(
      legend.text = element_text(size = 14),   # Increase legend text size
      legend.title = element_text(size = 16),  # Increase legend title size
      legend.key.size = unit(1.5, "lines")     # Increase legend key size (symbols)
    ) +
    scale_x_date(minor_breaks=monthbreaks, breaks=yearbreaks,labels=sprintf("%02d", seq(12, 24, 1)) )
  
  plot1.fitted3 <-  ggplot() +
    geom_point(aes(x=date[11:666], y=diseases[[d]][11:666], colour="Actual")) + 
    geom_line(aes(x=date[11:666], y=pred[[10]][1:656], colour="Fitted")) +
    scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
    ylab("Case Counts") +
    xlab("Year") +
    labs(title=paste0("Fitted ", dict[d], " Case Counts (10 Weeks Ahead)")) +
    theme(legend.position = "top") +
    theme_bw() + 
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          axis.ticks.length = unit(5, "pt")) +
    theme(axis.text.x=element_text(hjust=-0.7)) +
    theme(
      legend.text = element_text(size = 14),   # Increase legend text size
      legend.title = element_text(size = 16),  # Increase legend title size
      legend.key.size = unit(1.5, "lines")     # Increase legend key size (symbols)
    ) +
    scale_x_date(minor_breaks=monthbreaks, breaks=yearbreaks,labels=sprintf("%02d", seq(12, 24, 1)) )
  
  if(d=="Dengue") { 
    allplots[[d]] <- ggarrange(plot1.fitted1, plot1.fitted2, plot1.fitted3, ncol=3, common.legend=T)
  } else {
    allplots[[d]] <- ggarrange(plot1.fitted1, plot1.fitted2, plot1.fitted3, ncol=3, legend=F)
    }
  

}

png(file="../plots/s1-FittedvsActual.png", width=16, height=12, units="in", res=300)
ggarrange(plotlist=allplots, nrow=4)
dev.off()

################################################################
#################### overfit/underfit plot #####################
################################################################
load("./data.Rdata")
dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox")

bestmodel <- c()
for(d in colnames(diseases)){
  load(paste0("./out/", d, "_waic.Rdata"))
  waic <- lapply(waic.list, "[[", 1)
  bestmodel[d] <- names(which.min(waic))
}

plotlist <- list()

#fix the plots - why is it giving all the same plot 

for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/2c_", d, "insamplepoint.Rdata"))
  pred <- get(paste0("in.pred.", tolower(bestmodel[d])))
  diseases[d] <- as.numeric(diseases[[d]])
  
  plot1.fit1 <- ggplot() +
    geom_point(aes(x=diseases[2:666,d],y=pred[[1]][1:665], colour="red"), show.legend=F) +
    geom_abline()+
    xlab("Actual Case Counts") +
    ylab("Fitted Case Counts") +
    labs(title=paste0("1 Week Ahead (", dict[d], ")")) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
          aspect.ratio=1) +
    coord_equal(ratio=1)

  
  plot1.fit2 <- ggplot() +
    geom_point(aes(x=diseases[6:666,d],y=in.pred.bal2[[5]][1:661], colour="red"), show.legend=F) +
    geom_abline()+
    xlab("Actual Case Counts") +
    ylab("Fitted Case Counts") +
    labs(title=paste0("5 Weeks Ahead (", dict[d], ")")) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
          aspect.ratio=1) +
    coord_equal(ratio=1)
  
  plot1.fit3 <- ggplot() +
    geom_point(aes(x=diseases[11:666,d],y=pred[[10]][1:656], colour="red"), show.legend=F) +
    geom_abline()+
    xlab("Actual Case Counts") +
    ylab("Fitted Case Counts") +
    labs(title=paste0("10 Weeks Ahead (", dict[d], ")")) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title=element_text(face="bold"),
          plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"),
          aspect.ratio=1) +
    coord_equal(ratio=1)

  plotlist[[d]] <- plot_grid(plot1.fit1, plot1.fit2, plot1.fit3, ncol=1)
  
}




png(file="../plots/s2-underoverpredict.png", width=14, height=10, units="in", res=300)
plot_grid(plotlist=plotlist, ncol=4)
dev.off()
