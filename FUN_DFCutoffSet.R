DFCutoffSet <- function(variables) {
  ## Set cutoff of N years survival
  cutoffKM <- cutoff3
  OSTimeSetting <- year3


  return(df)
}

## Presetting
Tar="mayoscore5"
time = "time"
censor="censor"

#### Dataframe with cutoff setting ####
## Set cutoff of N years survival
cutoffKM <- ROCResultSeq_mayo5[["cutoff"]][["5_years"]]
OSTimeSetting <- 5

## Set 2 Groups for KM curve from Index
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
mayo2 <- mayo
mayo2$ReTime <- ""
mayo2$Status <- ""

for (i in c(1:length(mayo2[,1]))) {
  if(mayo2[i,time] >= OSTimeSetting*365){
    mayo2$ReTime[i] <- OSTimeSetting*365
    if(mayo2[i,censor]==0){
      mayo2$Status[i] <- 0
    }else{
      mayo2$Status[i] <- NA
    }

  }else if(mayo2[i,time] < OSTimeSetting*365){
    mayo2$ReTime[i] <- mayo2[i,time]
    mayo2$Status[i] <- mayo2[i,censor]
  }else{
    mayo2$ReTime[i] <- c("error")
    mayo2$Status[i] <- c("error")
  }
}

mayo2$ReTime <- mayo2$ReTime %>% as.numeric()
mayo2$Status <- mayo2$Status %>% as.numeric()
# mayo2 <- mayo2[!is.na(mayo2$Status),]

##### KM plot ####
fit<- survfit(Surv(ReTime, Status) ~ 1, data = mayo2)
ggsurvplot(fit)

fit2<- survfit(Surv(ReTime, Status) ~ ROCGroup, data = mayo2)
ggsurvplot(fit2)

