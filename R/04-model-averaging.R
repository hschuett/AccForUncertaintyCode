# Code for averaging previously fit models
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# WARNING: Running this script sequentially will take a long time (serveral days).
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports and Setup -----------------------------------------------------------------
library(loo)
options(mc.cores = parallel::detectCores())



# Loading lok_liks ------------------------------------------------------------------
model_names <- c("industryyear", "industryyear-ext", "industryyear-int",
                 "firm", "firm-ext", "firm-int")
file_names <- sapply(model_names, function(x) paste0("out/modelfits/", x, "-loklik.rds"))
K <- length(file_names)

log_lik_list <- vector("list", length=K)
for (k in 1:K){
  print(file_names[k])
  log_lik_list[[k]] <- readRDS(file_names[k])
}

wts <- loo_model_weights(log_lik_list, method="stacking")
print(wts)
# Method: stacking
# ------
#   weight
# model1 0.042
# model2 0.000
# model3 0.106
# model4 0.017
# model5 0.574
# model6 0.261

moweights <- data.frame(
  model = model_names,
  StackWeight = as.vector(wts)
  )
#              model  StackWeight
# 1     industryyear 4.238248e-02
# 2 industryyear-ext 1.476753e-07
# 3 industryyear-int 1.056742e-01
# 4             firm 1.742363e-02
# 5         firm-ext 5.736472e-01
# 6         firm-int 2.608724e-01

write.csv(moweights, "out/model-averaging-weights.csv", row.names=F)
