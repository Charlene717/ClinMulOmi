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
    source("FUN_DFCutoffSet.R")

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

  ## Dataframe with cutoff setting
  mayo <- DFCutoffSet(mayo, cutoffKM = ROCResultSeq_mayo5[["cutoff"]][["5_years"]],
                       OSTimeSetting = 5,
                       Tar="mayoscore5", time = "time", censor="censor")

  ## KM plot

  # Draw survival curves without grouping
  fit_all <- survfit(Surv(ReTime, Status) ~ 1, data = mayo)
  ggsurvplot(fit_all)

  # Draw survival curves with grouping
  fit <- survfit(Surv(ReTime, Status) ~ ROCGrp, data = mayo)
  # Basic plots
  ggsurvplot(fit)
  # Add p-value,
  ggsurvplot(fit, pval = TRUE, pval.coord = c(100, 0.03))
  # Add confidence interval
  ggsurvplot(fit, pval = TRUE, pval.coord = c(100, 0.03),
             conf.int = TRUE) # Add confidence interval
  # Add number at risk table
  ggsurvplot(fit, pval = TRUE, pval.coord = c(100, 0.03),
             conf.int = TRUE, risk.table = TRUE)

  ##
  Tar="mayoscore5"
  OSTimeSetting=5
  ggsurvplot(fit,
             ## Setting of main Fig
             # # Change font size, style and color at the same time
             title=paste0(Tar," (",OSTimeSetting," years survival)"),
             # main = "Survival curve", # No function
             # font.main = c(16, "bold", "darkblue"),
             # font.x = c(14, "bold.italic", "darkblue"),
             # font.y = c(14, "bold.italic", "darkblue"),
             # font.tickslab = c(12, "plain", "darkgreen"),
             # legend = c(0.2, 0.2), # Change legend posistion
             legend.title = paste0(Tar," (ROC)"),
             legend.labs = c(paste0(Tar,"_High"), paste0(Tar,"_Low")),

             size = 1,  # change line size
             # linetype = "strata", # change line type by groups
             palette = c("#ef476f", "#0077b6"), # custom color palette
             conf.int = TRUE, # Add confidence interval
             pval = TRUE, # Add p-value,
             pval.coord = c(100, 0.03), # Change p-value posistion
             xlim = c(0, 5*365), # Change x axis limits

             ## Set risk.table
             risk.table = TRUE, # Add risk table
             break.time.by = 250, # break time axis by 250
             risk.table.col = "strata" # Change risk table color by groups
  ) -> P.KM

  P.KM

  ## Export PDF
  pdf(file = paste0(Save.Path,"/",ProjectName,"_KMCurve.pdf"),
      width = 7,  height = 7
  )
  P.KM %>% print()
  dev.off()

################################################################

