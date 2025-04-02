library(tidyverse)
library(ggplot2)
library(ggprism)
library(cowplot)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")

######################################################
################### plot 1: cases ####################
######################################################
load("data.Rdata")
raw <- read.csv("climatedisease.csv")
ind <- dplyr::select(raw, "Start", "End")
diseases <- cbind(ind[13:678,], diseases)
diseases <- diseases[-4]

dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox")

date <- ind[["Start"]]
date <- as.Date(date, tryFormats = c("%m/%d/%Y"))
yearbreaks <- date[match(unique(year(date)), year(date))]
monthbreaks <- as.Date(c())
for(year in unique(year(date))){
  subsetdates <- date[year(date)==year]
  monthbreaks <- c(monthbreaks, subsetdates[match(unique(month(subsetdates)), month(subsetdates))])
}

diseases$date <- date[13:678]
p1_df <- pivot_longer(diseases, names_to="Disease", values_to="Cases", cols=3:6) #hard coded
p1_df$Disease <- dict[p1_df$Disease]

cases <- ggplot(p1_df) +
  geom_line(aes(x=date, y=Cases, group="Disease")) +
  theme_bw() + 
  facet_wrap(vars(Disease), scales="free") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.length = unit(5, "pt"),
        axis.text.x=element_text(hjust=-0.4),
        strip.background=element_rect(fill="lightcyan")) +
  labs(x="Year", y="Case Counts") +
  scale_x_date(breaks=yearbreaks, labels=sprintf("%02d", seq(12, 24, 1))) 
  #guides(x=guide_axis(minor.ticks=T))


pdf(file="../plots/fig1-cases.pdf", width=10, height=6, paper="a4r")
cases
dev.off()

png(file="../plots/fig1-cases.png", width=10, height=6, units="in", res=300)
cases
dev.off()

######################################################
########### plot 2: weights of best model ###########
######################################################
#1 get best model for each disease
#plot weights of each model
#plot coef beside and join

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
  if(d=="Salmonella")next
  load(paste0("./out/2b_", d, "mcmcsummary.Rdata"))
  
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




######################################################
########## plot 3: out of sample  ############
######################################################

load("data.Rdata")
plotlist <- list()

dict <- c("Dengue" = "Dengue", "Salmonella"="Salmonella", "URTI"="URTI", "Conjunc"= "Conjunctivitis",
          "Cpox" = "Chickenpox")
diseases <- diseases[-2]

for(d in colnames(diseases)){
  load(file=paste0("./out/", d, "outsamplemeasures.Rdata"))
  
  #MAPE
  mae1 <- gather(data.frame(mapedf), key=Model, value=MAPE, 1:13)
  mae1$MAPE <- as.numeric(mae1$MAPE) *100
  mae1 <- cbind(mae1, c(rep("Midas", 144), rep("No Midas", 12)))
  mae1 <- cbind(V1=c(1:12), mae1)
  colnames(mae1)[4] <- "type"
  
  colours <- c("royalblue1", "royalblue2", "royalblue3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "seagreen2",
               "seagreen3", "seagreen4", "orange1", "orange2", "orange3", "#E50000")
  
  plot2.mae <- ggplot(data=mae1, aes(x=V1, y=MAPE, group=Model, colour=factor(Model))) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="none") +
    geom_line(linewidth=1) +
    geom_point(aes(shape=type), size=3) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 4)) +
    labs(colour="Model", shape="") +
    xlab("Forecast Horizon (Weeks)") +
    ylab("MAPE (%)")  +
    scale_x_continuous(breaks=c(1,4,8,12), minor_breaks = 1:12)
  
  #PIT
  pit.ks <- pit.ks[1:12,]
  pit.ks.p <- pit.ks.p[1:12,]
  pit1 <- gather(data.frame(pit.ks), key=Model, value=Value, 1:13)
  pit1 <- cbind(pit1, c(rep("Midas", 144),rep("No Midas", 12)))
  pit1 <- cbind(c(1:12), pit1)
  colnames(pit1)[4] <- "type"
  colnames(pit1)[1] <- "V1"
  
  plot2.pit <- ggplot(data=pit1, aes(x=V1, y=Value, group=Model, colour=factor(Model))) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="none") +
    geom_line(linewidth=1) +
    geom_point(aes(shape=type), size=3) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 4)) +
    labs(colour="Model", shape="") +
    xlab("Forecast Horizon (Weeks)") +
    ylab("KS Test") +
    scale_x_continuous(breaks=c(1,4,8,12), minor_breaks = 1:12)
  
  
  #Logscore
  ls.mat <- ls.mat[1:12,]
  ls1 <- gather(data.frame(ls.mat), key=Model, value=`Logarithmic Score`, 1:13)
  ls1 <- cbind(ls1, c(rep("Midas", 144), rep("No Midas", 12)))
  ls1 <- cbind(c(1:12), ls1)
  colnames(ls1)[4] <- "type"
  colnames(ls1)[1] <- "V1"
  
  plot2.ls <- ggplot(data=ls1, aes(x=V1, y=`Logarithmic Score`, group=Model, colour=factor(Model))) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="none") +
    geom_line(linewidth=1) +
    geom_point(aes(shape=type), size=3) +
    scale_colour_manual(values=colours) +
    scale_shape_manual(values=c(1, 4)) +
    labs(colour="Model", shape="") +
    xlab("Forecast Horizon (Weeks)") +
    ylab("Logarithmic Score") +
    scale_x_continuous(breaks=c(1,4,8,12), minor_breaks = 1:12)
  
  #CRPS
  crps1 <- gather(data.frame(crps.mat), key=Model, value=`Continuous Ranked Probability Score`, 1:13)
  crps1 <- cbind(crps1, c(rep("Midas", 144), rep("No Midas", 12)))
  crps1 <- cbind(c(1:12), crps1)
  colnames(crps1)[4] <- "type"
  colnames(crps1)[1] <- "V1"
  #crpsmodels <- c("AR-(1)", "AR-(2)", "AR-(3)", "BAL-(1)", "BAL-(2)", "BAL-(3)", "BGL-(1)", "BGL-(2)", "BGL-(3)",
  #                "BL-(1)", "BL-(2)", "BL-(3)", "DLM")
  modelnames <- c("AR-(N)", "AR-(D)", "AR-(H)", "BAL-(N)", "BAL-(D)", "BAL-(H)", "BSGL-(N)", "BSGL-(D)", "BSGL-(H)",
                  "BL-(N)", "BL-(D)", "BL-(H)", "NMIDAS")
  
  
  plot2.crps <- ggplot(data=crps1, aes(x=V1, y=`Continuous Ranked Probability Score`, group=Model, colour=factor(Model))) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 12)) +
    geom_line(linewidth=1) +
    geom_point(aes(shape=type), size=3) +
    scale_colour_manual(values=colours, labels=modelnames) +
    #scale_colour_discrete(labels=crpsmodels) +
    scale_shape_manual(values=c(1, 4)) +
    labs(colour="Model", shape="") +
    xlab("Forecast Horizon (Weeks)") +
    ylab("CRPS") +
    scale_x_continuous(breaks=c(1,4,8,12), minor_breaks = 1:12)
  
  k <- which(colnames(diseases)==d) - 1
  
  finalplot <- plot_grid(plot2.mae, plot2.pit, plot2.ls, plot2.crps +
                          theme(legend.position="none"), labels=LETTERS[(1+k*4):(4+k*4)], nrow=1)
  label <- ggdraw() + draw_label(dict[d], size=20)
  plotlist[[d]] <- plot_grid(label, finalplot, ncol=1, rel_heights=c(0.1,1))
  savelegend <- ggpubr::get_legend(plot2.crps) 
  }

alld <- plot_grid(plotlist=plotlist, ncol=1)


pdf(file="../plots/fig3-outsampleperf.pdf", width=15, height=15)
plot_grid(alld, legend=savelegend, ncol=2, rel_widths=c(3,0.3))
dev.off()

png(file="../plots/fig3-outsampleperf.png", width=15, height=15, units="in", res=300)
plot_grid(alld, legend=savelegend, ncol=2, rel_widths=c(3,0.3))
dev.off()

######################################################
################ plot 4: DM test  ###################
######################################################


##diebold mariano test##
#DM test (1 step ahead)
#darker colour -> lower pval -> x is better than y
dmplot <- list()
for(d in colnames(diseases)){
  if(d=="Salmonella") next
  load(paste0("./out/", d, "outsamplemeasures.Rdata"))
  colnames(dmpdf)[1] <- "AR-(N)"
  colnames(dmpdf) <- gsub(colnames(dmpdf), pattern="BGL", replacement="BSGL")
  dmdatp <- dmpdf[2:14,]
  diag(dmdatp) <- NA
  dmdatlongp <- gather(dmdatp, key=x, value=pvalue)
  dmdatlongp <- cbind(dmdatlongp, y=rep(colnames(dmdatp), 13))
  dmdatlongp[which(dmdatlongp$pvalue < 0.05), "signif"] <- "*"
  
  colnames(dmdf)[1] <- "AR-(N)"
  colnames(dmdf) <- gsub(colnames(dmdf), pattern="BGL", replacement="BSGL")
  dmdat <- dmdf[2:14,]
  diag(dmdat) <- NA
  dmdatlong <- gather(dmdat, key=x, value="Test Statistic")
  dmdatlong <- cbind(dmdatlong, y=rep(colnames(dmdat), 13))
  dmdatlong$signif <- dmdatlongp$signif
  dmdatlong$`Test Statistic` <- as.numeric(dmdatlong$`Test Statistic`)
  
  
  dmplot[[d]] <-local({
    ggplot(data=dmdatlong, aes(x=x, y=y, fill=`Test Statistic`)) +
      geom_tile() +
      geom_text(aes(label=signif, colour="p-value < 0.05"), size=10) +
      labs(x="Base Model", y="Competing Model", title=dict[d]) +
      scale_colour_manual(name="", values="black") +
      scale_fill_gradient2(low = "navy", mid="lavender", high = "firebrick") +
      annotate("rect", xmin=c(12.5, 0.5), xmax=c(13.5, 13.5), ymin=c(0.5,12.5), ymax=c(13.5,13.5), colour="coral", fill="transparent", linewidth=1) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.title=element_text(face="bold"),
            axis.text.x = element_text(angle = 45, hjust=1)) +
      guides(color= guide_legend(override.aes = list(label = "*", size = 10)))
    
    
  }) 
}

pdf(file=paste0("../plots/fig4-", "dieboldmariano.pdf"), width=11, height=11)
ggarrange(plotlist=dmplot, common.legend=T, labels="AUTO")
dev.off()

png(file=paste0("../plots/fig4-", "dieboldmariano.png"), width=11, height=11, units="in", res=300)
ggarrange(plotlist=dmplot, common.legend=T, labels="AUTO")
dev.off()

######################################################
######### plot 5: weights of simulation coefs ########
######################################################
load("data.Rdata")
sapply(paste0("./functions/", list.files("./functions/")), source)

get_weights <- function(m) {
  simres <- read.csv("./tables/sim_results.csv")

  simres <- simres[grepl("t", simres$Parameter),]
  simres <- dplyr::filter(simres, Model==m)
  simres$Mean <- rowMeans(simres[3:102])
  simres$Lo <- apply(simres[3:102],1,quantile,probs=0.025)
  simres$Hi <- apply(simres[3:102],1,quantile,probs=0.975)
  
  simresp <- read.csv("./tables/sim_results.csv")
  simresp <- simresp[grepl("p", simresp$Parameter),]
  simresp <- dplyr::filter(simresp, Model==m)
  simresp$Mean <- rowMeans(simresp[3:102])
  
  sim_coef <- list()
  
  sim_theta <- list()
  sim_theta[[1]] <- c(simres[1, "Mean"], simres[2, "Mean"])
  sim_theta[[2]] <- c(simres[3, "Mean"], simres[4, "Mean"])
  sim_theta[[3]] <- c(simres[5, "Mean"], simres[6, "Mean"])
  sim_theta[[4]] <- c(simres[7, "Mean"], simres[8, "Mean"])
  sim_coef[["mean"]] <- midasWeights(sim_theta, 1)
  # 
  # for(i in 1:4){
  #   sim_coef[["mean"]][[i]] <- simresp$Mean[i] * sim_coef[["mean"]][[i]]
  # }
  # 
  # real_theta <- list()
  # real_theta[[1]] <- c(2/3, 3/2.5) #rain
  # real_theta[[2]] <- c(1/1.5, 4.2/3.4) #meanT
  # real_theta[[3]] <- c(3.6/0.7, 3.2/2.3) #minT
  # real_theta[[4]] <- c(1.3/1.6, 4) #maxT
  # 
  # real_coef <- midasWeights(real_theta, 1)
  # 
  # real_p <- c(3, 7.8, 4.8, 5.8)
  # for(i in 1:4){
  #   real_coef[[i]] <- real_coef[[i]] * real_p[i]
  # }
  # par(mfrow=c(2,2))
  # for(i in 1:4){
  #   plot(x=1:30, y=real_coef[[i]], col="red", type="l")
  #   lines(x=1:30, y=sim_coef[["mean"]][[1]], type="l")
  #   lines(x=1:30, y=rep(simresp$Mean[1], 30))
  # }
  # 
  
  sim_theta_lo <- list()
  sim_theta_lo[[1]] <- c(simres[1, "Lo"], simres[2, "Lo"])
  sim_theta_lo[[2]] <- c(simres[3, "Lo"], simres[4, "Lo"])
  sim_theta_lo[[3]] <- c(simres[5, "Lo"], simres[6, "Lo"])
  sim_theta_lo[[4]] <- c(simres[7, "Lo"], simres[8, "Lo"])
  sim_coef[["lo"]] <- midasWeights(sim_theta_lo, 1)
  
  sim_theta_hi <- list()
  sim_theta_hi[[1]] <- c(simres[1, "Hi"], simres[2, "Hi"])
  sim_theta_hi[[2]] <- c(simres[3, "Hi"], simres[4, "Hi"])
  sim_theta_hi[[3]] <- c(simres[5, "Hi"], simres[6, "Hi"])
  sim_theta_hi[[4]] <- c(simres[7, "Hi"], simres[8, "Hi"])
  sim_coef[["hi"]] <- midasWeights(sim_theta_hi, 1)
  
  return(sim_coef)
}

simweightplotter <- function() {
  nlags <-30

  real_theta <- list()
  real_theta[[1]] <- c(2/3, 3/2.5) #rain
  real_theta[[2]] <- c(1/1.5, 4.2/3.4) #meanT
  real_theta[[3]] <- c(3.6/0.7, 3.2/2.3) #minT
  real_theta[[4]] <- c(1.3/1.6, 4) #maxT
  
  real_coef <- midasWeights(real_theta, 1)
  
  #best model coefs (BGL-N)
  vars <- c("Total Precipitation", "Max Temp", "Min Temp", "Mean Temp")
  coefdf <- matrix(nrow=600, ncol=4)
  colnames(coefdf) <- c("Days", "Variable", "Value", "Model")
  coefdf[,"Days"] <- rep(1:30, 4)
  coefdf[,"Variable"] <- rep(vars, each=30)
  coefdf[, "Model"] <-  rep(c("97.5% UB", "BSGL-(N)", "2.5% LB", "Real", "NMIDAS"), each=120)
  coefdf <- as.data.frame(coefdf)
  
  coefdf[which(coefdf$Model=="NMIDAS"), "Value"] <- 3.33
  coefdf[which(coefdf$Model=="Real"), "Value"] <- unlist(real_coef) * 100
  coefdf[which(coefdf$Model=="BSGL-(N)"), "Value"] <- unlist(get_weights("in.summary.bgl1")[["mean"]]) * 100
  coefdf[which(coefdf$Model=="2.5% LB"), "Value"] <- unlist(get_weights("in.summary.bgl1")[["lo"]]) * 100
  coefdf[which(coefdf$Model=="97.5% UB"), "Value"] <- unlist(get_weights("in.summary.bgl1")[["hi"]]) * 100
  
  coefdf$Value <- as.numeric(coefdf$Value)
  coefdf$Days <- as.numeric(coefdf$Days)
  
  maxcoef <- coefdf %>%
    group_by(Model, Variable) %>%
    summarise(max_val = max(Value), max_day = Days[which.max(Value)]) %>%
    filter(Model != "NMIDAS")

  ggplot(coefdf) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_blank()) +
    geom_line(aes(x=Days, y=Value, color=Model, group=Model, linetype=Model), linewidth=1) +
    geom_label_repel(data=maxcoef, 
                     show.legend=F,
                     aes(x=max_day, y=max_val, group= Model, colour = Model,
                                       label = paste0("Day ", max_day, ", ", round(max_val, 1), "%"))) +
    facet_wrap(~Variable) +
    labs(y="Contribution of days to 1-week ahead cases (%)") +
    scale_color_manual(values=c("Real"="brown","97.5% UB"="lightseagreen", "BSGL-(N)"="steelblue", 
                                "2.5% LB"="darkturquoise", "NMIDAS"="darkseagreen")) +
    scale_linetype_manual(values=c("Real"="solid","97.5% UB"="dashed", "BSGL-(N)"="solid", 
                                   "2.5% LB"="dashed", "NMIDAS"="solid"))
}
simweightplotter()

