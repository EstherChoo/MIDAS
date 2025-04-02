#combine all separate timept files
#[ar/bal/bgl/bl][timept].Rdata
#"ar1", "ar2", "ar3", "bal1", "bal2", "bal3", "bgl1", "bgl2", "bgl3", "bl1", "bl2", "bl3", "nmidasar", "nmidasbal", "nmidasbl", "nmidasagg"
#obj 1: pt pred - list of 15 model, each model list of 20 horizons, each horizon is 524x1
#obj 2: dens pred - list of 15 model, each model list of 20 horizons, each horizon is 524x4000

#setwd("C:/Users/esthe/Downloads/Acads/FYP/MyCode")
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/MIDAS/Code")

#loading files
files <- list.files("./rollingcv")
load("data.Rdata")
rowstart <- round(0.7 * nrow(diseases)) #to put for loop
rowseq <- seq(rowstart, nrow(diseases), 1)

#set up structure of output data
out.ptpred <- list()
out.denspred <- list()
for(i in 1:13){
  out.ptpred[[i]] <- list()
  out.denspred[[i]] <- list()
}
for(i in 1:13){
  for(j in 1:12){
    out.ptpred[[i]][[j]] <- matrix(nrow=192, ncol=1)
    out.denspred[[i]][[j]] <- matrix(nrow=192, ncol=4000)
  }
}
names(out.ptpred) <- c("AR-(N)", "AR-(D)", "AR-(H)", "BAL-(N)", "BAL-(D)", "BAL-(H)", "BGL-(N)", "BGL-(D)", "BGL-(H)", "BL-(N)",
                       "BL-(D)", "BL-(H)", "NMIDAS")
names(out.denspred) <- names(out.ptpred)

#TO CHANGE (ar-0, bal-3, bgl-6, bl-9, nmidasagg-13)
models <- c("ar", "bal", "bgl", "bl")
inds <- c(0, 3, 6, 9)

for(d in colnames(diseases)){
  for(m in 1:4){
    model <- models[m]
    ind <- inds[m]
    for(i in rowseq){
      load(file=paste0("./rollingcv/",d, "_", model, i, ".Rdata"))
      assign("lst", get(paste0(model, "pred", i)))
      
      n <- i - 474
      for(mod in 1:3){
        if(mod==2 & m==5) break
        for(hor in 1:12){
          #out.ptpred[[(mod+ind)]][[hor]][n:(n+9), 1] <- lst[[mod]][[hor]][[1]][1:10]
          #out.denspred[[(mod+ind)]][[hor]][n:(n+9),] <- matrix(lst[[mod]][[hor]][[2]])
          out.ptpred[[(mod+ind)]][[hor]] <- lst[[mod]][[hor]][[1]]
          out.denspred[[(mod+ind)]][[hor]] <- matrix(lst[[mod]][[hor]][[2]], ncol=4000, byrow=F)
        }
      }
      rm(list=paste0(model, "pred", i))
    }
    
  }
  
  load(file=paste0("./rollingcv/",d, "_nmidasagg", i, ".Rdata"))
  for(hor in 1:12){
    out.ptpred[["NMIDAS"]][[hor]] <- nmidasagg[[hor]][[1]]
    out.denspred[["NMIDAS"]][[hor]] <- nmidasagg[[hor]][[2]]
  }
  
  
  
  save(out.ptpred, file=paste0("./out/4b_", d, "outptpred.Rdata"))
  save(out.denspred, file=paste0("./out/4b_", d, "outdenspred.Rdata"))
}
