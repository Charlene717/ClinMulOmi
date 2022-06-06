#### Time-dependent-ROC
## Ref: https://www.bioinfo-scrounger.com/archives/Time-dependent-ROC/
## Paper: https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0332-6

#### To-do list ####
# - [] ggplot: Time-dependent AUC plot

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Current path and new folder setting* #####
  ProjectName = "Example"
  Sampletype = "mayo"
  #ProjSamp.Path = paste0(Sampletype,"_",ProjectName)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",Sampletype)
  Save.Path = paste0(getwd(),"/",Version)

  Save_Path = save.path
  ## Create new folder
  if (!dir.exists(Save_Path)){
    dir.create(Save_Path)
  }

## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Parameter setting* #####
  timesseq.set <- seq(from=1, to=10, by=2)  ## timesseq.set <- c(1,3,5,7,9)

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

  ##### Function setting #####
    ## Call function
    source("FUN_TimeDepROC.R")

##### Load database ####
  data(mayo)

##### TimeDepROC #####

  ROCResultSeq_mayo4 <- TimeDepROC(mayo,timesseq.set,Tar="mayoscore4",time = "time", censor="censor",
                          save.path = Save.Path , Filename="Seq_mayo4")

  ROCResultSeq_mayo5 <- TimeDepROC(mayo,timesseq.set,Tar="mayoscore5",time = "time", censor="censor",
                             save.path = Save.Path , Filename="Seq_mayo5")


  ROCResult_5Y <- TimeDepROC(mayo,5,Tar="mayoscore4",time = "time", censor="censor",
                                save.path = Save.Path , Filename="5years_mayo5")

  ROCResult_10Y <- TimeDepROC(mayo,10,Tar="mayoscore5",time = "time", censor="censor",
                            save.path = Save.Path , Filename="10years_mayo5")

  ## Compare two time-dependent AUC
  plotAUCcurve(ROCResultSeq_mayo4[["time_roc_res"]], conf.int=TRUE, col="red")
  plotAUCcurve(ROCResultSeq_mayo5[["time_roc_res"]], conf.int=TRUE, col="blue", add=TRUE)
  legend("bottomright",c("mayoscore4", "mayoscore5"), col = c("red","blue"), lty=1, lwd=2)
  dev.off()

  ## Export PDF
  pdf(file = paste0(Save.Path,"/",ProjectName,"_CompareAUC.pdf"),
      width = 7,  height = 7
  )
    plotAUCcurve(ROCResultSeq_mayo4[["time_roc_res"]], conf.int=TRUE, col="red")
    plotAUCcurve(ROCResultSeq_mayo5[["time_roc_res"]], conf.int=TRUE, col="blue", add=TRUE)
    legend("bottomright",c("mayoscore4", "mayoscore5"), col = c("red","blue"), lty=1, lwd=2)
  dev.off()

  # ## (Pending) ## df for ggplot
  # ## df for ggplot
  # result.confint <- confint(object = ROCResultSeq_mayo5[["time_roc_res"]], level = 0.95,
  #                           n.sim = 3000)



####################################################################################################

  ##### Plot the KM curve #####
  ## Ref: https://blog.yjtseng.info/post/2020-05-13-survival-curve/
  ## Ref: http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
  ## Ref: https://www.rdocumentation.org/packages/survminer/versions/0.4.9/topics/ggsurvplot
  library(survival)
  library(survminer)



  #### Plot KM curve ####
  fit<- survfit(Surv(time, censor) ~ 1, data = mayo)
  ggsurvplot(fit)

  # OS
  fit_OS <- survfit(Surv(Survival.in.days, Vital.status) ~ NKIndexROC, data = Bulk_GSE13255_4_M1_PhenoROC2)
  fit_OS_P <- surv_pvalue(fit_OS)
  fit_OS_ggplot <-  ggsurvplot(fit_OS, pval = TRUE,
                               risk.table = TRUE,
                               legend.title = paste0(Marker,"(ROC)"),
                               legend.labs = c(paste0(Marker,"_High"), paste0(Marker,"_Low")),
                               title=paste0("OS"," (",OSTimeSetting," years survival)"),
                               log.rank.weights="1",xlim = c(0, OSTimeSetting*365))
  fit_OS_ggplot


