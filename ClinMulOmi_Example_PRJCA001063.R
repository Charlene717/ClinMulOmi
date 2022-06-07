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


##### Identify conserved cluster markers  #####
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

  #####
  DimPlot(scRNA.SeuObj, reduction = "umap")

  scRNA.SeuObj2 <- scRNA.SeuObj

  DimPlot(scRNA.SeuObj2, reduction = "umap",group.by ="Cell_type", label = TRUE, pt.size = 0.5) + NoLegend()

  #####

  ## Identify conserved cluster markers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  set.seed(1) # Fix the seed
  scCTMarker <- FindAllMarkers(scRNA.SeuObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #scCTMarker <- FindAllMarkers(scRNA.SeuObj)

  ## Ref: https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html


  Cell_type.set <- scRNA.SeuObj@meta.data[["Cell_type"]] %>% unique()

  Idents(scRNA.SeuObj) <- scRNA.SeuObj$Cell_type
  FibcCTMarker <- FindMarkers(scRNA.SeuObj, ident.1= Cell_type.set[1],only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

########################################################################
  ##### Deconvolution #####
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(immunedeconv)
  library(tibble)

  gene_expression_matrix <- read.delim(file = "TCGA_PAAD_HiSeqV2")
  row.names(gene_expression_matrix) <- gene_expression_matrix[,1]
  gene_expression_matrix <- gene_expression_matrix[,-1]
  res_quantiseq = deconvolute_base_custom(
    gene_expression_matrix,
    signature,
    n_permutations = 100,
    log10 = TRUE)
  res_quantiseq <- res_quantiseq %>% as.data.frame()
  res_quantiseq$cell_type <-  rownames(res_quantiseq)

  res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(res_quantiseq)))
  res_quantiseq %>%
    gather(sample, score, -cell_type) %>%
    ggplot(aes(x=sample, y=score, color=cell_type)) +
    geom_point(size=4) +
    facet_wrap(~cell_type, scales="free_x", ncol=3) +
    scale_color_brewer(palette="Paired", guide=FALSE) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




