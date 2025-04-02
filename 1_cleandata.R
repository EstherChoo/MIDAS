library(tidyverse)
library(zoo)
library(tseries)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")
sapply(paste0("./functions/", list.files("./functions/")), source)

####Load data####
data <- read.csv(file="./climatedisease.csv")

diseases <- dplyr::select(data, Dengue, 
                          Salmonellosis_non_enteric_fevers, Acute_Upper_Respiratory_Tract_infections, Acute_Conjunctivitis,
                          Chickenpox, HFMD)
colnames(diseases) <- c("Dengue", "Salmonella", "URTI", "Conjunc", "Cpox")
climate <- dplyr::select(data, starts_with("Drt"), starts_with("MT"), starts_with("MinT"), starts_with("MaxT"))
climate_low <- dplyr::select(climate, ends_with("epiweek"))
climate_high <- dplyr::select(climate, !ends_with("epiweek"))
climate_low$Drt_epiweek <- na.approx(climate_low$Drt_epiweek)

for(n in 1:ncol(climate_high)){ #interpolate using na.approx, if at start of ts then use avg of future 7 days to back-interpolate
  app <- na.approx(climate_high[n])
  if(length(app)<678){
    avg <- mean(app[1:7])
    pad <- 678-length(app) 
    app <- c(rep(avg,pad) , app)
  }
  climate_high[n] <- app
}

climate_high <- climate_high[-c(1:12),]
climate_high_lst <- list()
for(i in 1:4){
  climate_high_lst[[i]] <- as.matrix(climate_high[,(1+(i-1)*30):(30+(i-1)*30)])
}

#lag diseases
dis_lag <- data.frame(rep(NA, nrow(diseases)))
dises <- colnames(diseases)
c <- 1
for(dis in dises){
  for(i in 1:12){
    c <- c + 1
    dis_lag <- cbind(dis_lag, dplyr::lag(diseases[dis],i))
    colnames(dis_lag)[c] <- paste0(dis, "_lag", i)
  }
}
dis_lag <- na.omit(dis_lag[-1])
for(c in 1:4){
  for(i in 1:4){
    climate_low <- cbind(climate_low, dplyr::lag(climate_low[c],i))
    colnames(climate_low)[ncol(climate_low)] <- paste0(colnames(climate_low)[c], "_lag", i)
  }
}
climate_low <- climate_low[-c(1:12),]
climate_low <- climate_low[-c(1:4)]

diseases <- diseases[-c(1:12),]

#stationary test
for(i in 1:6){
  print(adf.test(diseases[[i]]))
}

adf.test(diff(diseases[[3]]))


save(list=c("dis_lag", "climate_low", "climate_high_lst", "diseases"), file="data.Rdata")





