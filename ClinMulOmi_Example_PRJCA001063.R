## Ref: https://github.com/sqjin/CellChat
## Ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages #####
  #### Basic and BiocManager installation ####
  source("FUN_Package_InstLoad.R")
  FUN_Basic.set <-c("tidyverse","timeROC","Seurat","survival","survminer","survivalROC","colorspace")
  FUN_BiocManager.set <- c("basilisk","zellkonverter","SeuratDisk")

  FUN_Package_InstLoad(Basic.set = FUN_Basic.set, BiocManager.set = FUN_BiocManager.set)
  rm(FUN_Basic.set, FUN_BiocManager.set)


##### Function setting #####
  ## Call function
  source("FUN_TimeDepROC.R")
  source("FUN_DFCutoffSet.R")

##### Current path and new folder setting*  #####
  ProjectName = "ClinMulOmi"
  Filename = "PADC"
  Version = paste0(Sys.Date(),"_",ProjectName,"_PADC")
  Save.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(Save.Path)){
    dir.create(Save.Path)
  }

##### Load Data #####
## Load Seurat RData
  # load("PRJCA001063_seuratObject.RData")
  load("D:/Dropbox/#_Dataset/Cancer/PDAC/2022-11-24_scRNA_SeuObj_PDAC_ROGUE.RData")

  scRNA.SeuObj$Cell_type <- scRNA.SeuObj$singleR_classic_PredbyscRNA2

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
  set_cibersort_binary("CIBERSORT.R")
  set_cibersort_mat("LM22.txt")
  library(tibble)

  ## Man pages for omnideconv/immunedeconv
  # https://rdrr.io/github/omnideconv/immunedeconv/man/

  # res_quantiseq = deconvolute_consensus_tme_custom(as.matrix(BulkGE.df),
  #                                                  VariableFeatures(object = scRNA.SeuObj))

  res_quantiseq = deconvolute_epic_custom(BulkGE.df, CTMarker.df, VariableFeatures(object = scRNA.SeuObj))

  # ## Cibersort
  # ## Ref: https://omnideconv.org/immunedeconv/articles/immunedeconv.html#special-case-cibersort
  # BulkGE.df[BulkGE.df==0] <- 0.00001
  # CTMarker.df[CTMarker.df==0] <- 0.00001
  # res_quantiseq = deconvolute_cibersort_custom(BulkGE.df, CTMarker.df)

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
                    scale_x_discrete(limits = rev(levels(res_quantiseq))) -> p.Decon1
  p.Decon1

  res_quantiseq %>% gather(sample, score, -cell_type) %>%
                    ggplot(aes(x=sample, y=score, color=cell_type)) +
                    geom_point(size=4) +
                    facet_wrap(~cell_type, scales="free_x", ncol=3) +
                    scale_color_brewer(palette="Paired", guide=FALSE) +
                    coord_flip() +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> p.Decon2
  p.Decon2

  pdf(file = paste0(Save.Path,"/Decon_",Filename,".pdf"),
      width = 10,  height = 10
  )
    p.Decon1
    p.Decon2
  dev.off()
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
  colnames(res_quantiseq_Temp) <- gsub(" ", ".", colnames(res_quantiseq_Temp))


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
  timesseq.set <- c(1,3,5)

  ROCResultSeq <- TimeDepROC(res_quantiseq_Temp,timesseq.set,
                             Tar = "Ductal.cell.type.2", #as.character(Cell_type.set[3])
                             time = "OStime", censor="OS",
                             save.path = Save.Path , Filename="PAAD")

  ## All cell type
  Cell_type.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique() %>% as.character()
  Cell_type.set <- gsub(" ", ".", Cell_type.set)

  ROCResultSeq.lt <- list()
  ## Export PDF
  pdf(file = paste0(Save.Path,"/ROC_ALL_",Filename,".pdf"),
      width = 10,  height = 10
  )
    for (i in 1:length(Cell_type.set)) {
      ROCResultSeq.lt[[Cell_type.set[i]]]<- TimeDepROC(res_quantiseq_Temp,timesseq.set,
                                                       Tar = Cell_type.set[i], #as.character(Cell_type.set[3])
                                                       time = "OStime", censor="OS",
                                                       save.path = Save.Path ,
                                                       Filename = paste0("PAAD_",
                                                                       gsub("\\.", " ",Cell_type.set[i])))%>% print()
    }

  dev.off()


  ## Compare two time-dependent AUC
  pdf(file = paste0(Save.Path,"/",ProjectName,"_CompareAUC.pdf"),
      width = 17,  height = 7
  )
  for (i in 1:length(Cell_type.set)) {
    if(i==1){
      plotAUCcurve(ROCResultSeq.lt[[i]][["time_roc_res"]], conf.int=TRUE, col=rainbow(length(Cell_type.set))[i])
      #legend("topright",Cell_type.set[i], col = rainbow(length(Cell_type.set))[i], lty=1, lwd=2)

    }else{
      plotAUCcurve(ROCResultSeq.lt[[i]][["time_roc_res"]], conf.int=TRUE, col=rainbow(length(Cell_type.set))[i], add=TRUE)
      #legend("topright",Cell_type.set[i], col = rainbow(length(Cell_type.set))[i], lty=1, lwd=2)
    }
  }
  legend("topleft",gsub("\\.", " ",Cell_type.set), col = rainbow(length(Cell_type.set)), lty=1, lwd=2)

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

  Target = "Ductal.cell.type.2"
  OSTimeSetting = 1

  ## Dataframe with cutoff setting
  res_quantiseq_Temp <- DFCutoffSet(res_quantiseq_Temp, OSTimeSetting = 1,
                                    cutoffKM = ROCResultSeq.lt[[Target]][["cutoff"]][["1_years"]],
                                    Tar = Target, time = "OStime", censor="OS")

  ## KM plot
  # Draw survival curves without grouping
  fit_all <- survfit(Surv(ReTime, Status) ~ 1, data = res_quantiseq_Temp)
  ggsurvplot(fit_all)

  # Draw survival curves with grouping
  fit <- survfit(Surv(ReTime, Status) ~ ROCGrp, data = res_quantiseq_Temp)
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

  ## Beautify Figure
  Target2 <- gsub("\\.", " ",Target)
  ggsurvplot(fit,
             ## Setting of main Fig
             # # Change font size, style and color at the same time
             title=paste0(Target2," (",OSTimeSetting," years survival)"),
             # main = "Survival curve", # No function
             # font.main = c(16, "bold", "darkblue"),
             # font.x = c(14, "bold.italic", "darkblue"),
             # font.y = c(14, "bold.italic", "darkblue"),
             # font.tickslab = c(12, "plain", "darkgreen"),
             # legend = c(0.2, 0.2), # Change legend posistion
             legend.title = paste0(Target2," (ROC)"),
             legend.labs = c(paste0(Target2,"_High"), paste0(Target2,"_Low")),

             size = 1,  # change line size
             # linetype = "strata", # change line type by groups
             palette = c("#ef476f", "#0077b6"), # custom color palette
             conf.int = TRUE, # Add confidence interval
             pval = TRUE, # Add p-value,
             pval.coord = c(100, 0.03), # Change p-value posistion
             xlim = c(0, OSTimeSetting*365), # Change x axis limits

             ## Set risk.table
             risk.table = TRUE, # Add risk table
             break.time.by = 250, # break time axis by 250
             risk.table.col = "strata" # Change risk table color by groups
  ) -> P.KM

  P.KM

  ## Export PDF
  pdf(file = paste0(Save.Path,"/",ProjectName,"_",Target,"_KMCurve.pdf"),
      width = 8,  height = 7
  )
  P.KM %>% print()
  dev.off()

  ## All cell type
  pdf(file = paste0(Save.Path,"/",ProjectName,"_KMCurve.pdf"),
      width = 8,  height = 7
  )
  for (i in 1:length(Cell_type.set)) {
    Target = Cell_type.set[i]
    OSTimeSetting = 1

    ## Dataframe with cutoff setting
    res_quantiseq_Temp <- DFCutoffSet(res_quantiseq_Temp, OSTimeSetting = 1,
                                      cutoffKM = ROCResultSeq.lt[[Target]][["cutoff"]][["1_years"]],
                                      Tar = Target, time = "OStime", censor="OS")
    # Draw survival curves with grouping
    fit <- survfit(Surv(ReTime, Status) ~ ROCGrp, data = res_quantiseq_Temp)

    ## Beautify Figure
    Target2 <- gsub("\\.", " ",Target)
    ggsurvplot(fit,
               ## Setting of main Fig
               # # Change font size, style and color at the same time
               title=paste0(Target2," (",OSTimeSetting," years survival)"),
               legend.title = paste0(Target2," (ROC)"),
               legend.labs = c(paste0(Target2,"_High"), paste0(Target2,"_Low")),

               size = 1,  # change line size
               # linetype = "strata", # change line type by groups
               palette = c("#ef476f", "#0077b6"), # custom color palette
               conf.int = TRUE, # Add confidence interval
               pval = TRUE, # Add p-value,
               pval.coord = c(100, 0.03), # Change p-value posistion
               xlim = c(0, OSTimeSetting*365), # Change x axis limits

               ## Set risk.table
               risk.table = TRUE, # Add risk table
               break.time.by = 250, # break time axis by 250
               risk.table.col = "strata" # Change risk table color by groups
    ) -> P.KM

    print(P.KM)

  }
  dev.off()

  ##### Cox proportional hazards model #####
  ## Ref: https://www.datacamp.com/tutorial/survival-analysis-R
  ## Ref: https://www.reneshbedre.com/blog/survival-analysis.html
  ## Ref: https://cran.r-project.org/web/packages/survminer/vignettes/ggforest-show-interactions-hazard-ratio.html
  ## Ref: https://rpkgs.datanovia.com/survminer/reference/ggforest.html

  ##** Ref: https://wangcc.me/LSHTMlearningnote/cox.html#cox%E5%9B%9E%E6%AD%B8%E6%A8%A1%E5%9E%8B%E4%B8%AD%E5%8C%85%E6%B6%B5%E7%9A%84%E5%81%87%E8%A8%AD
  ##** Ref: https://codingnote.cc/zh-tw/p/13034/

  # Fit a Cox proportional hazards model
  surv_object <- Surv(time = res_quantiseq_Temp$OStime, event = res_quantiseq_Temp$OS)
  fit.coxph <- coxph(surv_object ~  Stellate.cell + Macrophage.cell + Endothelial.cell + T.cell + B.cell+
                                    Endocrine.cell + Fibroblast.celly + Macrophage.cell +
                                    Ductal.cell.type.2 + Endocrine.cell + Ductal.cell.type.1 + Acinar.cell,
                     data = res_quantiseq_Temp)


  # fit.coxph <- coxph(surv_object ~ res_quantiseq_Temp[,5] + res_quantiseq_Temp[,6] + res_quantiseq_Temp[,7] +
  #                    Ductal.cell.type.2,
  #                    data = res_quantiseq_Temp)
  ggforest(fit.coxph, data = res_quantiseq_Temp,
           noDigits = 3)
  pdf(file = paste0(Save.Path,"/",ProjectName,"_HR.pdf"),
      width = 8,  height = 7
  )
  ggforest(fit.coxph, data = res_quantiseq_Temp,
           noDigits = 3)
  dev.off()

  #### Save RData ####
  # save.image("ClinMulOmi_Example_PRJCA001063.RData")
  save.image(paste0("D:/Dropbox/#_Dataset/Cancer/PDAC/",Version,"_ClinMulOmi.RData"))

#######################################################################################
# #### CNV ####
#   ##### Heatmap plotting #####
#   ## Set column annotation
#   column_ha_T = HeatmapAnnotation(
#     Condition = anno_colum.df$Body.weight,
#     Condition2 = anno_colum.df$Area,
#     col = list(Condition = c("Mild"="#9b6ab8", "Severe"="#6e6970"),
#                Condition2= c("Mild"="#9b6ab8", "Severe"="#6e6970","Medium"="#b57545")),
#     show_legend = T
#   )
#
#   ## generate color of top annotation
#   col_exp <-  colorRamp2(
#     c(0, 0.025, 0.05),
#     c("#248a5c", "#52bf8e","#bbedd7")
#   )
#   col_exp2 <-  colorRamp2(
#     c(-17, 0, 17),
#     c("#1a5691", "#96cbff", "#d1e8ff")
#   )
#
#   row_ha = rowAnnotation(
#     p.value = anno_row.df$p.value.Severe.vs..Mild.,
#     LogFC = anno_row.df$Difference.Severe.vs..Mild.,
#     col = list(p.value = col_exp, LogFC = col_exp2),
#     show_legend = T
#   )
#
#
#   Heatmap(
#     matrix.df[-1],
#     # column_title = target_gene,
#     # column_title_side = "top",
#     cluster_rows = T,
#     cluster_columns = T,
#     show_column_names = F,
#     show_row_names = F,
#     name = "M-Value",
#     # set color
#     col = colorRamp2(
#       c(0, 0.5, 1),
#       c("#1c77d9", "#1a2938", "#ffe182")
#     ),
#     show_heatmap_legend = T,
#     use_raster = F,
#     top_annotation = column_ha_T,
#     right_annotation = row_ha
#   ) -> P.Heatmap
#
#   P.Heatmap %>% print
#
#
#   ##### Export PDF #####
#
#   pdf(
#     file = paste0(getwd(), "/",Version,"/", Sys.Date(), "_MValue_Heatmap.pdf"),
#     width = 7, height = 7
#   )
#   P.Heatmap
#
#   graphics.off()
#
#
#   ## From: D:\Dropbox\##_GitHub\##_PHH_Lab\Heatmap_Cachexia_DNAMeth
