# Code for collecting the results of all previous steps
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# WARNING: Running this script sequentially will take a long time (serveral days).
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports ---------------------------------------------------------------------------
library(matrixStats)
library(dplyr)
library(haven)
library(tidyr)
library(readr)


compute_postior_moments <- function(posterior_draws, mname) {
  moments <- data.frame(id=1:ncol(posterior_draws))
  x <- as.matrix(posterior_draws)
  moments[[paste0(mname,"PostMean")]] <- colMeans(x)
  moments[[paste0(mname,"PostSD")]]   <- colSds(x)
  moments[[paste0(mname,"Post01")]]   <- colQuantiles(x, probs=0.01)
  moments[[paste0(mname,"Post99")]]   <- colQuantiles(x, probs=0.99)
  moments$id <- NULL
  return(moments)
}


# Load data ---------------------------------------------------------------
weight_table <- read_csv("out/model-averaging-weights.csv")
ols_data <- read_csv("out/all-OLS-predictions.csv")



# Model averaging posterior predictions ---------------------------------------------
short_names <- c("IndYr", "IndYrE", "IndYrI",
                 "Firm", "FirmE", "FirmI")

av_TA_co_preds <- matrix(0, nrow=2000, ncol=nrow(ols_data))
for (i in 1:nrow(weight_table)) {
  filename <- weight_table$model[i]
  print(filename)
  TA_co_preds <- readRDS(paste0("out/modelfits/", filename, "-ta-coef-only-preds.rds"))
  TA_co_preds <- TA_co_preds * weight_table$StackWeight[i]
  av_TA_co_preds  <- av_TA_co_preds + TA_co_preds
  remove(TA_co_preds)
}
saveRDS(av_TA_co_preds, "out/modelfits/averaged-ta-coef-only-preds.rds")


av_comparisons <- ols_data %>%
  mutate(AveragedPostMean = colMeans(av_TA_co_preds),
         AveragedPostSD = colSds(av_TA_co_preds),
         AveragedPost01 = colQuantiles(av_TA_co_preds, probs=0.01),
         AveragedPost99 = colQuantiles(av_TA_co_preds, probs=0.99))

write_dta(av_comparisons, "out/all-model-predictions.dta")
