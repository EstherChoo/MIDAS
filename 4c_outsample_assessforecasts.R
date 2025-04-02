library(tidyverse)
library(forecast)
library(EnvStats)
library(scoringRules)

setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/code")
sapply(paste0("./functions/", list.files("./functions/")), source)

load("./data.Rdata")

for(d in colnames(diseases)){
  load(paste0("./out/4b_", d, "outptpred.Rdata"))
  load(paste0("./out/4b_", d, "outdenspred.Rdata"))

  #1. mape 2. rmaebbb  3. dm test 4. pit+ks 5. logscore 6. crps
  #1-3: ar, 4-6: bal, 7-9: bgl, 10-12:bl, 13: nmidas
  
  model.names <-c("A  R-(N)", "AR-(D)", "AR-(H)", "BAL-(N)", "BAL-(D)", "BAL-(H)", "BGL-(N)", "BGL-(D)", "BGL-(H)", "BL-(N)",
                  "BL-(D)", "BL-(H)", "NMIDAS")
  
  calc.err <- function(ypred, yactual){ #func to calculate residuals
    h <- length(ypred)
    res <- list()
    mape <- list()
    for(i in 1:h){
        realyactual <- as.numeric(yactual[(i+1):length(yactual)])
        realypred <- as.numeric(ypred[[i]][1:(length(ypred[[i]])-i)])
        #realyactual <- descale(og.y, yactual[(1+i):length(yactual)])
        #realypred <- descale(og.y, ypred[[i]][1:(length(ypred[[i]])-i)])
        res[[i]] <- abs(ypred[[i]][1:(length(ypred[[i]])-i)] - yactual[(1+i):length(yactual)])
        mape[[i]] <- mean(abs((realyactual-realypred)/realyactual))
        mae[[i]] <- mean(abs(realyactual-realypred))
    }
    return(list(res, mape, mae))
  }
  
  testy <- as.numeric(diseases[475:666, d])
  
  #residuals and mean abs error
  i <- 1
  error <- list()
  mape <- list()
  mae <- list()
  for(model in out.ptpred){
    result <- calc.err(model, testy)
    error[[i]] <- result[[1]]
    mape[[i]] <- result[[2]]
    mae[[i]] <- result[[3]]
    
    i <- i+1
  }
  mapedf <- sapply(mape, c)
  colnames(mapedf) <- model.names
  write.csv(mapedf, paste0("./tables/", d, "_mape.csv"))
  
  #relative mean abs error
  rmae <- list()
  for(h in 1:12){
    mat <- matrix(0, nrow=13, ncol=13) #total 16 models
    for(i in 1:length(out.ptpred)){
      for(j in i:length(out.ptpred)){
        mat[i,j] <- mae[[i]][[h]]/mae[[j]][[h]] #val < 1 => model i (row) better than model j (col)
        mat[j,i] <- mae[[j]][[h]]/mae[[i]][[h]]
      }
    }
    rmae[[h]] <- mat
  }
  
  newdf <- data.frame(t(model.names))
  
  for(i in 1:length(rmae)){
    newdf <- rbind(newdf, rep(i,13))
    colnames(rmae[[i]]) <- names(newdf)
    newdf <- rbind(newdf, rmae[[i]])
  }
  
  colnames(newdf) <- model.names
  newdf <- newdf[-1,]
  write.csv(newdf, paste0("./tables/", d, "_rmae.csv"))
  
  
  #diebold mariano test
  dm <- list()
  dmp <- list()
  for(h in 1:12){
    dmmat <- matrix(0, nrow=13, ncol=13)
    dmmatp <- matrix(0, nrow=13, ncol=13)
    for(i in 1:length(out.ptpred)){
      for(j in 1:length(out.ptpred)){
        if(i==j){
          next
        }
        dmmat[i,j] <- dm.test(error[[i]][[h]][1:192], error[[j]][[h]][1:192], h=1, alternative="g")$statistic #is model in col better than model in row
        dmmatp[i,j] <- dm.test(error[[i]][[h]][1:192], error[[j]][[h]][1:192], h=1, alternative="g")$p.value
      }
    }
    dm[[h]] <- dmmat
    dmp[[h]] <- dmmatp
  }
  
  #dm test statistic
  dmdf <- data.frame(t(model.names))
  
  for(i in 1:length(dm)){
    dmdf <- rbind(dmdf, rep(i,13))
    colnames(dm[[i]]) <- names(dmdf)
    dmdf <- rbind(dmdf, dm[[i]])
  }
  
  colnames(dmdf) <- model.names
  dmdf <- dmdf[-1,]
  write.csv(dmdf, paste0("./tables/", d, "_dieboldmariano.csv"))
  
  #dm test p value
  dmpdf <- data.frame(t(model.names))
  
  for(i in 1:length(dmp)){
    dmpdf <- rbind(dmpdf, rep(i,13))
    colnames(dmp[[i]]) <- names(dmpdf)
    dmpdf <- rbind(dmpdf, dmp[[i]])
  }
  
  colnames(dmpdf) <- model.names
  dmpdf <- dmpdf[-1,]
  write.csv(dmpdf, paste0("./tables/", d, "dieboldmarianopvalues.csv"))
  
    
  ##DENSITY PERFORMANCE
  #pit
  pit <- function(mat, actualy, h, unif){ #takes in matrix of timepts x density - each col is a density for 1 timepoint
    pitval <- c()
    for(n in 1:ncol(mat)){
      pitval <- c(pitval, pemp(actualy[n+h-1], mat[,n]))
    }
    return(ks.test(pitval, unif))
  }
  
  pit.ks <- matrix(0, nrow=12, ncol=13)
  pit.ks.p <- matrix(0, nrow=12, ncol=13)
  i <- 1
  unif <- runif(1000)
  for(file in out.denspred){
    for(h in 1:12){
      ks <- pit(file[[h]], testy, h=1, unif) #get pit values
      pit.ks[h,i] <- ks$statistic #check pit values against uniform using ks
      pit.ks.p[h,i] <- ks$p.value
    }
    i <- i + 1
  }
  
  colnames(pit.ks) <- model.names
  colnames(pit.ks.p) <- model.names
  
  #add signif to pit.ks.p
  symp <- symnum(pit.ks.p, corr = FALSE,
                 cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c("***","**","*over","."," "))
  
  for(row in 1:12){
    for(col in 1:13){
      pit.ks.p[row, col] <- paste0("$", pit.ks.p[row,col], "^{", symp[row, (col-1)], "}$")
    }
  }
  
  write.csv(pit.ks, paste0("./tables/", d, "pit.csv"))
  write.csv(pit.ks.p, paste0("./tables/", d, "pitpvalues.csv"))
  
  #log score and crps
  ls.mat <- matrix(0, nrow=12, ncol=13)
  crps.mat <- matrix(0, nrow=12, ncol=13)
  i <- 1
  for(i in 1:13){
    file <- out.denspred[[i]]
    for(h in 1:12){
      mat <- file[[h]]
      mat <- mat[-c((nrow(mat)-h+1):nrow(mat)),]
      yactual <- testy[(h+1):length(testy)]
      ls <- logs_sample(yactual, mat)
      ls <- mean(ls[which(is.finite(ls))])
      crps <- mean(crps_sample(yactual, mat))
      ls.mat[h,i] <- ls
      crps.mat[h,i] <- crps
    }
  }
  
  colnames(ls.mat) <- model.names
  colnames(crps.mat) <- model.names
  write.csv(ls.mat, paste0("./tables/",d, "_logscore.csv"))
  write.csv(crps.mat, paste0("./tables/", d, "_crps.csv"))
   
  save(list=c("mapedf", "rmae", "dmdf", "dmpdf", "ls.mat", "crps.mat", "pit.ks", "pit.ks.p"), file=paste0("./out/", d, "outsamplemeasures.Rdata"))
}
