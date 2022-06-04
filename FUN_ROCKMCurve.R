## Ref:https://www.bioinfo-scrounger.com/archives/Time-dependent-ROC/

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)


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


##### Load datasets  #####
library(timeROC)
library(survival)
library(survivalROC)

data(mayo)

time_roc_res <- timeROC(
  T = mayo$time,
  delta = mayo$censor,
  marker = mayo$mayoscore5,
  cause = 1,
  weighting="marginal",
  times = c(3 * 365, 5 * 365, 10 * 365),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res$AUC
confint(time_roc_res, level = 0.95)$CI_AUC

plot(time_roc_res, time=3 * 365, col = "red", title = FALSE)
plot(time_roc_res, time=5 * 365, add=TRUE, col="blue")
plot(time_roc_res, time=10 * 365, add=TRUE, col="green")
legend("bottomright",c("3 Years" ,"5 Years", "10 Years"),
       col=c("red", "blue", "green"), lty=1, lwd=2)

# 圖形美化
time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2],
  TP_10year = time_roc_res$TP[, 3],
  FP_10year = time_roc_res$FP[, 3]
)
library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )


## 比較兩個time-dependent AUC
time_roc_res2 <- timeROC(
  T = mayo$time,
  delta = mayo$censor,
  marker = mayo$mayoscore4,
  cause = 1,
  weighting="marginal",
  times = c(3 * 365, 5 * 365, 10 * 365),
  ROC = TRUE,
  iid = TRUE
)

time_roc_res2$AUC
compare(time_roc_res, time_roc_res2, adjusted = TRUE)$p_values_AUC

plotAUCcurve(time_roc_res, conf.int=TRUE, col="red")
plotAUCcurve(time_roc_res2, conf.int=TRUE, col="blue", add=TRUE)
legend("bottomright",c("mayoscore5", "mayoscore4"), col = c("red","blue"), lty=1, lwd=2)


## ROC的最佳閾值（cutoff）
mayo$mayoscore5[which.max(time_ROC_df$TP_3year - time_ROC_df$FP_3year)]

##------------## Ch ##------------##
## 要怎麼標記點?


