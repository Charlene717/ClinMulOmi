#### Time-dependent-ROC
## Ref: https://www.bioinfo-scrounger.com/archives/Time-dependent-ROC/
## Paper: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0332-6

## To-do list ##
# - [] ggplot: Time-dependent AUC plot

TimeDepROC <- function(mayo,timesseq.set,Tar="mayoscore5",time = "time", censor="censor",
                       timeROC.lt = list(cause = 1, weighting="marginal", ROC = TRUE, iid = TRUE),
                       save.path = "", Filename="") {
  ##### Load Packages #####
    #### Basic installation ####
    Package.set <- c("tidyverse","timeROC","survival","survivalROC","colorspace")
    ## Check whether the installation of those packages is required from basic
    for (i in 1:length(Package.set)) {
      if (!requireNamespace(Package.set[i], quietly = TRUE)){
        install.packages(Package.set[i])
      }
    }
    ## Load Packages
    lapply(Package.set, library, character.only = TRUE)
    rm(Package.set,i)


  ##### timeROC #####
  formals(timeROC)[names(timeROC.lt)] <- timeROC.lt

  time_roc_res <- timeROC(
    T = mayo[,time],
    delta = mayo[,censor],
    marker = mayo[,Tar],
    times = timesseq.set*365
  )

  time_roc_res$AUC
  confint(time_roc_res, level = 0.95)$CI_AUC

  ##### Plot ROC Curve #####
  #### Color setting reference ####
    ## Color ## Ref: https://www.jianshu.com/p/21971df8e2e4
    # library(wesanderson)
    # names(wes_palettes)

    ## Color ## Ref: https://bookdown.org/wangminjie/R4DS/tidyverse-ggplot2-colors.html
    ## Color ## Ref: https://cran.r-project.org/web/packages/colorspace/vignettes/colorspace.html

  #### Basic ROC plot ####
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


  #### Beautify ROC plot ####
    ## Dataframe for ggplot
    if(length(timesseq.set)==1){
      for (i in 1:length(timesseq.set)) {
        time_ROC_df <- data.frame( time_roc_res$TP[, 2], time_roc_res$FP[, 2])
        colnames(time_ROC_df) <- c(paste0("TP_",timesseq.set[i],"years"),paste0("FP_",timesseq.set[i],"years"))
      }
      rm(i)
    }else{
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

    }

    ## ggplot
    library(ggplot2)
    ## Using-geom-line-in-a-for-loop
    ## Ref: https://stackoverflow.com/questions/27852620/using-geom-line-in-a-for-loop

    P.ROC <-  ggplot(data = time_ROC_df)

    if(length(timesseq.set)==1){
      for (i in  1:length(timesseq.set)) {
        P.ROC <- P.ROC + geom_line(aes_string(x = time_ROC_df[,paste0("FP_",timesseq.set[i],"years")], y = time_ROC_df[,paste0("TP_",timesseq.set[i],"years")]), size = 1.5, color = qualitative_hcl(length(timesseq.set))[i])+
          annotate("text", x = 0.75, y = 0.55 -i*0.05, size = 7,
                   label = paste0("AUC at ",timesseq.set[i]," years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = qualitative_hcl(length(timesseq.set))[i]
          )

      }
      P.ROC
    }else{
      for (i in  1:length(timesseq.set)) {
        P.ROC <- P.ROC + geom_line(aes_string(x = time_ROC_df[,paste0("FP_",timesseq.set[i],"years")], y = time_ROC_df[,paste0("TP_",timesseq.set[i],"years")]), size = 1.5, color = qualitative_hcl(length(timesseq.set))[i])+
          annotate("text", x = 0.75, y = 0.55 -i*0.05, size = 7,
                   label = paste0("AUC at ",timesseq.set[i]," years = ", sprintf("%.3f", time_roc_res$AUC[[i]])), color = qualitative_hcl(length(timesseq.set))[i]
          )

      }
      P.ROC
    }

    P.ROC2 <- P.ROC + geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2)+
      theme_classic() + # White background
      labs(x = "False positive rate", y = "True positive rate") +
      theme(axis.line = element_line(colour = "black",
                                     size = 0, linetype = "solid")) + # Change the line type and color of axis lines
      theme(axis.text.x = element_text(face="bold", color="black",
                                       size=17, angle=0),
            axis.text.y = element_text(face="bold", color="black",
                                       size=17, angle=0)) +  # Change the appearance and the orientation angle
      ggtitle(paste0("ROC Curve (",gsub("\\.", " ",Tar),")"))+ # Change the main title and axis labels

      theme(
        plot.title = element_text(color="black", size=20, face="bold", hjust = 0.5),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"), # Change the color, the size and the face of  the main title, x and y axis labels
        aspect.ratio=1)+
      theme(panel.background = element_rect(fill = "white", colour = "black",  size = 2.5))

    print(P.ROC2)

  #### Basic time-dependent AUC plot ####
    ## timesseq setting
    maxandnext=function(x){list(max(x),max(x[-which.max(x)]))} ## Ref: https://bbs.pinggu.org/forum.php?mod=viewthread&action=printable&tid=1419015
    maxandnext.set <- maxandnext(res_quantiseq_Temp$OStime)
    timesseq.set2 <- seq(from=1, to=floor(as.numeric(maxandnext.set[2])/365), by=1) ## timesseq.set2 <- c(1,3,5,7,9)


    plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
    legend("bottomright",Tar, col = c("red"), lty=1, lwd=2)

    time_roc_res2 <- timeROC(
      T = mayo[,time],
      delta = mayo[,censor],
      marker = mayo[,Tar],
      cause = 1,
      weighting="marginal",
      times = timesseq.set2*365,
      ROC = TRUE,
      iid = TRUE
    )
    plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue")
    legend("bottomright",Tar, col = c("blue"), lty=1, lwd=2)

    # ## (Pending) ## df for ggplot
    # result.confint <- confint(object = time_roc_res, level = 0.95,
    #                           n.sim = 3000)

  #### Export PDF ####
    ## Create new folder
    if (!dir.exists(save.path)){
      dir.create(save.path)
    }

    pdf(file = paste0(save.path,"/ROC_",Filename,".pdf"),
        width = 10,  height = 10
    )
    # ROC
      print(P.ROC2)
      plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
      legend("bottomright",Tar, col = c("red"), lty=1, lwd=2)

      plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue")
      legend("bottomright",Tar, col = c("blue"), lty=1, lwd=2)
    dev.off()

  ##### Optimal threshold for ROC（cutoff）#####

    for (i in 1:length(timesseq.set)) {
      if(i==1){
        cutoff.df <- as.data.frame(mayo[,Tar][which.max(time_ROC_df[,paste0("TP_",timesseq.set[i],"years")] - time_ROC_df[,paste0("FP_",timesseq.set[i],"years")])])
        colnames(cutoff.df) <- paste0(timesseq.set[i],"_years")

      }else{
        cutoff_Temp.df <- as.data.frame(mayo[,Tar][which.max(time_ROC_df[,paste0("TP_",timesseq.set[i],"years")] - time_ROC_df[,paste0("FP_",timesseq.set[i],"years")])])
        colnames(cutoff_Temp.df) <- paste0(timesseq.set[i],"_years")
        cutoff.df <- cbind(cutoff.df,cutoff_Temp.df)
      }

    }
    rm(i,cutoff_Temp.df)
    rownames(cutoff.df) <- "Cutoff"

    Output <- list(time_roc_res = time_roc_res,
                   time_ROC_df = time_ROC_df,
                   cutoff = cutoff.df)

    return(Output)
    }




