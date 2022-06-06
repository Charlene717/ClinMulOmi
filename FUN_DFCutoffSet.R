
DFCutoffSet <- function(mayo, cutoffKM, OSTimeSetting,
                        Tar="mayoscore5", time = "time", censor="censor") {

  ##### Set 2 Groups for KM curve from Index #####
  # ROCGroupName <- paste0(Tar,"_",OSTimeSetting,"yr_","ROCGrp")
  ROCGroupName <- paste0("ROCGrp")

  mayo[,ROCGroupName]<- ""
  for (i in c(1:length(mayo[,1]))) {
    if(mayo[i,Tar] >= cutoffKM){
      mayo[,ROCGroupName][i] <- "High"
    }else if(mayo[i,Tar] < cutoffKM){
      mayo[,ROCGroupName][i] <- "Low"
    }else{
      mayo[,ROCGroupName][i] <- 3
    }
  }
  mayo <- mayo[!mayo[,ROCGroupName]==3,]

  ##### Revise the survival time and status #####
  ## !!!Please check the survfit package whether the setting of N-years survival has already included in this package
  # ReTimeName <- paste0(Tar,"_",OSTimeSetting,"yr_","time")
  # StatusName <- paste0(Tar,"_",OSTimeSetting,"yr_","Status")
  ReTimeName <- paste0("ReTime")
  StatusName <- paste0("Status")

  mayo[,ReTimeName]<- ""
  mayo[,StatusName]<- ""

  for (i in c(1:length(mayo[,1]))) {
    if(mayo[i,time] >= OSTimeSetting*365){
      mayo[,ReTimeName][i] <- OSTimeSetting*365
      if(mayo[i,censor]==0){
        mayo[,StatusName][i] <- 0
      }else{
        mayo[,StatusName][i] <- NA
      }

    }else if(mayo[i,time] < OSTimeSetting*365){
      mayo[,ReTimeName][i] <- mayo[i,time]
      mayo[,StatusName][i] <- mayo[i,censor]
    }else{
      mayo[,ReTimeName][i] <- c("error")
      mayo[,StatusName][i] <- c("error")
    }
  }

  mayo[,ReTimeName] <- mayo[,ReTimeName] %>% as.numeric()
  mayo[,StatusName] <- mayo[,StatusName] %>% as.numeric()
  # mayo <- mayo[!is.na(mayo[,StatusName]),]


  return(mayo)
}

