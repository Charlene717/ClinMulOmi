## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

#### Load the required libraries ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","timeROC","Seurat","survival","survminer","survivalROC","colorspace")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("basilisk","zellkonverter","SeuratDisk")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  options(stringsAsFactors = FALSE)

##### Function setting #####
  ## Call function
  source("FUN_TimeDepROC.R")
  source("FUN_DFCutoffSet.R")

##### Current path and new folder setting*  #####
  ProjectName = "TOP2A" # Secret, ECM, CC
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Load Data #####
## Load Seurat RData
  load("PRJCA001063_seuratObject.RData")

## Load Bulk RData
  BulkGE.df <- read.delim(file = "TCGA_PAAD_HiSeqV2")
  row.names(BulkGE.df) <- BulkGE.df[,1]
  BulkGE.df <- BulkGE.df[,-1]

## Load survival data
  Survival.df <- read.delim(file = "survival_PAAD_survival.txt")



##### Generation of signature matrices from single cell data #####
  ##### Prepossessing  #####
  seed=1
  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- ScaleData(scRNA.SeuObj, verbose = FALSE)
  # AnnoNames.set <- colnames(scRNA.SeuObj@meta.data)

  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- FindVariableFeatures(object = scRNA.SeuObj)

  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- RunPCA(scRNA.SeuObj, features = VariableFeatures(object = scRNA.SeuObj))

  PCAdims <- 30
  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- RunUMAP(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims)

  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- FindNeighbors(scRNA.SeuObj, reduction = "pca", dims = 1:PCAdims)
  set.seed(seed) # Fix the seed
  scRNA.SeuObj <- FindClusters(scRNA.SeuObj, resolution = 0.5)

  #### Plot UMAP ###
  DimPlot(scRNA.SeuObj, reduction = "umap")

  scRNA.SeuObj2 <- scRNA.SeuObj

  DimPlot(scRNA.SeuObj2, reduction = "umap",group.by ="Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()


  #### Find Markers ###
  # ## Identify conserved cluster markers
  # # find markers for every cluster compared to all remaining cells, report only the positive ones
  # set.seed(1) # Fix the seed
  # scCTMarker <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # #scCTMarker <- FindAllMarkers(scRNA.SeuObj)
  #
  # ## Ref: https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html


  Cell_type.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()
  Idents(scRNA.SeuObj) <- scRNA.SeuObj$Cell_type

  ## About 1 hour
  for (i in 1:length(Cell_type.set)) {
    if(i==1){
      CTMarker.df <- FindMarkers(scRNA.SeuObj, ident.1= Cell_type.set[i],only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25) %>%
                     rownames_to_column(var="Gene") %>%
                     select(Gene,avg_log2FC)
      colnames(CTMarker.df)[2] <- Cell_type.set[i] %>% as.character()

    }else{
      CTMarker_Temp <- FindMarkers(scRNA.SeuObj, ident.1= Cell_type.set[i],only.pos = TRUE,
                                   min.pct = 0.25, logfc.threshold = 0.25) %>%
        rownames_to_column(var="Gene") %>%
        select(Gene,avg_log2FC)
      colnames(CTMarker_Temp)[2] <- Cell_type.set[i] %>% as.character()

      CTMarker.df <- full_join(CTMarker.df,CTMarker_Temp, by = "Gene")
    }
  }
  rm(i,CTMarker_Temp)
  CTMarker.df_Ori <-CTMarker.df
  CTMarker.df[is.na(CTMarker.df)] <-0
  CTMarker.df <- CTMarker.df %>% column_to_rownames(var="Gene")
#************************************************************************************************************************#
  ##### Deconvolution #####
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(immunedeconv)
  library(tibble)

  # res_quantiseq = deconvolute_consensus_tme_custom(as.matrix(BulkGE.df),
  #                                                    VariableFeatures(object = scRNA.SeuObj))
  res_quantiseq = deconvolute_epic_custom(BulkGE.df, CTMarker.df, VariableFeatures(object = scRNA.SeuObj))
  # res_quantiseq = deconvolute_base_custom(BulkGE.df,
  #                                         CTMarker.df,
  #                                         n_permutations = 100,
  #                                         log10 = F)
  res_quantiseq <- res_quantiseq %>% as.data.frame()
  res_quantiseq$cell_type <- rownames(res_quantiseq)

  res_quantiseq %>% gather(sample, fraction, -cell_type) %>%
                    # plot as stacked bar chart
                    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
                    geom_bar(stat='identity') +
                    coord_flip() +
                    scale_fill_brewer(palette="Paired") +
                    scale_x_discrete(limits = rev(levels(res_quantiseq)))

  res_quantiseq %>% gather(sample, score, -cell_type) %>%
                    ggplot(aes(x=sample, y=score, color=cell_type)) +
                    geom_point(size=4) +
                    facet_wrap(~cell_type, scales="free_x", ncol=3) +
                    scale_color_brewer(palette="Paired", guide=FALSE) +
                    coord_flip() +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#************************************************************************************************************************#
  ##### TimeDepROC #####
  ## Data prepocessing
  res_quantiseq_Temp <- res_quantiseq %>% t() %>% as.data.frame() %>%
                        rownames_to_column(var="sample")
  res_quantiseq_Temp$sample <- gsub("\\.", "-", res_quantiseq_Temp$sample)

  ## Keep primary tumor only (01: Primary Solid Tumor) ## Ref: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  res_quantiseq_Temp <- res_quantiseq_Temp[grepl("01$", res_quantiseq_Temp$sample),]


  res_quantiseq_Temp <- left_join(Survival.df,res_quantiseq_Temp)
  res_quantiseq_Temp <- res_quantiseq_Temp[!is.na(res_quantiseq_Temp[,ncol(res_quantiseq_Temp)]),]


  res_quantiseq_Temp <- res_quantiseq_Temp[,-5:-11]
  colnames(res_quantiseq_Temp) <- gsub("\\.", "", colnames(res_quantiseq_Temp))
  colnames(res_quantiseq_Temp) <- gsub(" ", "", colnames(res_quantiseq_Temp))


  ## as.numeric
  # res_quantiseq_Temp <- data.frame(apply(res_quantiseq_Temp, 2, function(x) as.numeric(as.character(x))))
  for (i in 5:ncol(res_quantiseq_Temp)) {
    res_quantiseq_Temp[,i] <- as.numeric(res_quantiseq_Temp[,i])
  }
  rm(i)

  ## How to round a data.frame in R that contains some character variables?
  ## Ref: https://stackoverflow.com/questions/9063889/how-to-round-a-data-frame-in-r-that-contains-some-character-variables
  res_quantiseq_Temp <- data.frame(lapply(res_quantiseq_Temp, function(y) if(is.numeric(y)) round(y, 5) else y))


  ## timesseq setting
  maxandnext=function(x){list(max(x),max(x[-which.max(x)]))} ## Ref: https://bbs.pinggu.org/forum.php?mod=viewthread&action=printable&tid=1419015
  maxandnext.set <- maxandnext(res_quantiseq_Temp$OStime)
  timesseq.set <- seq(from=1, to=floor(as.numeric(maxandnext.set[2])/365), by=2) ## timesseq.set <- c(1,3,5,7,9)

  ROCResultSeq <- TimeDepROC(res_quantiseq_Temp,timesseq.set,
                             Tar = "Ductalcelltype2", #as.character(Cell_type.set[3])
                             time = "OStime", censor="OS",
                             save.path = Save.Path , Filename="PAAD")



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



#************************************************************************************************************************#
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






