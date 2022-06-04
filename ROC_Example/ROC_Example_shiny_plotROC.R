# https://cran.r-project.org/web/packages/plotROC/vignettes/examples.html
rm(list = ls()) # Clean variable

library(plotROC)
shiny_plotROC()


## Loading required package: ggplot2
set.seed(2529)
D.ex <- rbinom(200, size = 1, prob = .5)
M1 <- rnorm(200, mean = D.ex, sd = .65)
M2 <- rnorm(200, mean = D.ex, sd = 1.5)

test <- data.frame(D = D.ex, D.str = c("Healthy", "Ill")[D.ex + 1], 
                   M1 = M1, M2 = M2, stringsAsFactors = FALSE)

basicplot <- ggplot(test, aes(d = D, m = M1)) + geom_roc()
basicplot
