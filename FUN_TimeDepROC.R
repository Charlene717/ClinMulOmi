## Ref: https://www.bioinfo-scrounger.com/archives/Time-dependent-ROC/
## Paper: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0332-6

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "Example"
  Sampletype = "mayo"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  timesseq.set <- seq(from=1, to=10, by=2)

  data(mayo)

  ## Test fuction
  ROCResult <- TimeDepROC(mayo,timesseq.set,Tar="mayoscore5",time = "time", censor="censor")

TimeDepROC <- function(mayo,timesseq.set,Tar="mayoscore5",time = "time", censor="censor") {
  ##### Load Packages #####
    #### Basic installation ####
    Package.set <- c("tidyverse","timeROC","survival","survivalROC")
    ## Check whether the installation of those packages is required from basic
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        install.packages(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)


  ##### timeROC  #####
  time_roc_res <- timeROC(
    T = mayo[,time],
    delta = mayo[,censor],
    marker = mayo[,Tar],
    cause = 1,
    weighting="marginal",
    times = timesseq.set*365,
    ROC = TRUE,
    iid = TRUE
  )

  time_roc_res$AUC
  confint(time_roc_res, level = 0.95)$CI_AUC

  ## Plot
  ## Color ## Ref: https://www.jianshu.com/p/21971df8e2e4
  # library(wesanderson)
  # names(wes_palettes)

  ## Color ## Ref: https://bookdown.org/wangminjie/R4DS/tidyverse-ggplot2-colors.html
  library(colorspace)


  for (i in 1:length(timesseq.set)) {
    if(i==1){
      plot(time_roc_res, time=timesseq.set[i]* 365, col = rainbow(length(timesseq.set))[i], title = FALSE)
    }else{
      plot(time_roc_res, time=timesseq.set[i] * 365, add=TRUE, col=rainbow(length(timesseq.set))[i])
    }

  }
  legend("bottomright",paste0(timesseq.set," Years"),
         col=rainbow(length(timesseq.set)), lty=1, lwd=2)
  rm(i)


  ##### Beautify Figures #####
  for (i in 1:length(timesseq.set)) {
    if(i==1){
      time_ROC_df <- data.frame( time_roc_res$TP[, i], time_roc_res$FP[, i])
      colnames(time_ROC_df) <- c(paste0("TP_",timesseq.set[i],"years"),paste0("FP_",timesseq.set[i],"years"))

    }else{
      time_ROC_df_Temp <- data.frame( time_roc_res$TP[, i], time_roc_res$FP[, i])
      colnames(time_ROC_df_Temp) <- c(paste0("TP_",timesseq.set[i],"years"),paste0("FP_",timesseq.set[i],"years"))
      time_ROC_df <- cbind(time_ROC_df,time_ROC_df_Temp)
    }

  }
  rm(i,time_ROC_df_Temp)

  library(ggplot2)
  ## Using-geom-line-in-a-for-loop
  ## Ref: https://stackoverflow.com/questions/27852620/using-geom-line-in-a-for-loop

  P.ROC <-  ggplot(data = time_ROC_df)
  for (i in  1:length(timesseq.set)) {
    P.ROC <- P.ROC + geom_line(aes_string(x = time_ROC_df[,timesseq.set[i]+1], y = time_ROC_df[,timesseq.set[i]]), size = 1.5, color = qualitative_hcl(length(timesseq.set))[i])+
      annotate("text", x = 0.75, y = 0.55 -i*0.05, size = 7,
               label = paste0("AUC at ",timesseq.set[i]," years = ", sprintf("%.3f", time_roc_res$AUC[[i]])), color = qualitative_hcl(length(timesseq.set))[i]
      )

  }
  P.ROC

  P.ROC2 <- P.ROC + geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2)+
    theme_classic() + # White background
    labs(x = "False positive rate", y = "True positive rate") +
    theme(axis.line = element_line(colour = "black",
                                   size = 0, linetype = "solid")) + # Change the line type and color of axis lines
    theme(axis.text.x = element_text(face="bold", color="black",
                                     size=17, angle=0),
          axis.text.y = element_text(face="bold", color="black",
                                     size=17, angle=0)) +  # Change the appearance and the orientation angle
    ggtitle("ROC Curve")+ # Change the main title and axis labels

    theme(
      plot.title = element_text(color="black", size=20, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="black", size=20, face="bold"),
      axis.title.y = element_text(color="black", size=20, face="bold"), # Change the color, the size and the face of  the main title, x and y axis labels
      aspect.ratio=1)+
    theme(panel.background = element_rect(fill = "white", colour = "black",  size = 2.5))

  print(P.ROC2)

  plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
  legend("bottomright",colnames(mayo)[3], col = c("red"), lty=1, lwd=2)


  pdf(file = paste0(Save.Path,"/",ProjectName,"_ROC.pdf"),
      width = 7,  height = 7
  )
    P.ROC2
  dev.off()

  ##### Optimal threshold for ROC（cutoff）#####
  cutoff <- mayo[,3][which.max(time_ROC_df$TP_3year - time_ROC_df$FP_3year)]

  Output <- list(time_roc_res = time_roc_res,
                 time_ROC_df = time_ROC_df,
                 cutoff = cutoff)

  return(Output)
}




####################################################################################################


## Compare two time-dependent AUC
  time_roc_res2 <- timeROC(
    T = mayo$time,
    delta = mayo$censor,
    marker = mayo$mayoscore4,
    cause = 1,
    weighting="marginal",
    times = timesseq.set*365,
    ROC = TRUE,
    iid = TRUE
  )

  time_roc_res2$AUC
  compare(time_roc_res, time_roc_res2, adjusted = TRUE)$p_values_AUC

  plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
  plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue", add=TRUE)
  legend("bottomright",c("mayoscore5", "mayoscore4"), col = c("red","blue"), lty=1, lwd=2)
  dev.off()


  pdf(file = paste0(Save.Path,"/",ProjectName,"_AUC.pdf"),
      width = 7,  height = 7
  )

    plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
    plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue", add=TRUE)
    legend("bottomright",c("mayoscore5", "mayoscore4"), col = c("red","blue"), lty=1, lwd=2)
  dev.off()

  result.confint <- confint(object = time_roc_res, level = 0.95,
                            n.sim = 3000)

  pdf(file = paste0(Save.Path,"/",ProjectName,"_ROC_AUC.pdf"),
      width = 7,  height = 7
  )
    P.ROC

    plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
    plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue", add=TRUE)
    legend("bottomright",c("mayoscore5", "mayoscore4"), col = c("red","blue"), lty=1, lwd=2)
  dev.off()

## Optimal threshold for ROC（cutoff）
  mayo$mayoscore5[which.max(time_ROC_df$TP_3year - time_ROC_df$FP_3year)]



