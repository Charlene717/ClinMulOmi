#### Dataframe with cutoff setting ####
## Test Function
mayo3 <- DFCutoffSet(mayo, cutoffKM = ROCResultSeq_mayo5[["cutoff"]][["5_years"]],
                     OSTimeSetting = 5,
                     Tar="mayoscore5", time = "time", censor="censor")

DFCutoffSet <- function(mayo, cutoffKM, OSTimeSetting,
                        Tar="mayoscore5", time = "time", censor="censor") {

  ##### Set 2 Groups for KM curve from Index #####
  mayo$ROCGroup <- ""
  for (i in c(1:length(mayo[,1]))) {
    if(mayo[i,Tar] >= cutoffKM){
      mayo$ROCGroup[i] <- "High"
    }else if(mayo[i,Tar] < cutoffKM){
      mayo$ROCGroup[i] <- "Low"
    }else{
      mayo$ROCGroup[i] <- 3
    }
  }
  mayo <- mayo[!mayo$ROCGroup==3,]

  ##### Revise the survival time and status #####
  ## !!!Please check the survfit package whether the setting of N-years survival has already included in this package
  mayo$ReTime <- ""
  mayo$Status <- ""

  for (i in c(1:length(mayo[,1]))) {
    if(mayo[i,time] >= OSTimeSetting*365){
      mayo$ReTime[i] <- OSTimeSetting*365
      if(mayo[i,censor]==0){
        mayo$Status[i] <- 0
      }else{
        mayo$Status[i] <- NA
      }

    }else if(mayo[i,time] < OSTimeSetting*365){
      mayo$ReTime[i] <- mayo[i,time]
      mayo$Status[i] <- mayo[i,censor]
    }else{
      mayo$ReTime[i] <- c("error")
      mayo$Status[i] <- c("error")
    }
  }

  mayo$ReTime <- mayo$ReTime %>% as.numeric()
  mayo$Status <- mayo$Status %>% as.numeric()
  # mayo <- mayo[!is.na(mayo$Status),]


  return(mayo)
}




##### KM plot ####
fit<- survfit(Surv(ReTime, Status) ~ 1, data = mayo)
ggsurvplot(fit)

fit2<- survfit(Surv(ReTime, Status) ~ ROCGroup, data = mayo)
ggsurvplot(fit2)

