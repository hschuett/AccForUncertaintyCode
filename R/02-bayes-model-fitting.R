# Code for Fitting six Bayesian Jones Models
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# WARNING: Running this script sequentially will take a long time (serveral days).
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports and Setup ------------------------------------------------------------
library(haven)
library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = FALSE)

# code that fits the model and saves output. Squeezed into one function:
run_model <- function(modelname,
                      filename,
                      parameters,
                      data) {
  # fit model
  fit <- stan(filename, data=input_data,
              iter=2000, warmup=1000, chains=2, seed=444,
              verbose=FALSE)

  # Save Results
  bayes_predicted <- as.matrix(fit, par="y_fit")
  saveRDS(bayes_predicted, paste0("out/modelfits/", modelname, "-ta-predicions.rds"))
  rm(bayes_predicted)

  bayes_predicted <- as.matrix(fit, par="y_fit_wo")
  saveRDS(bayes_predicted, paste0("out/modelfits/", modelname, "-ta-coef-only-preds.rds"))
  rm(bayes_predicted)

  ind_loklik <- as.matrix(fit, par="log_lik")
  saveRDS(ind_loklik, paste0("out/modelfits/", modelname, "-loklik.rds"))
  rm(ind_loklik)

  bayes_ind_coef <- as.matrix(fit, par="b")
  write.csv(bayes_ind_coef, gzfile(paste0("out/modelfits/", modelname, "-coef.csv.gz")))
  rm(bayes_ind_coef)

  bayes_results <- as.matrix(fit, par=parameters)
  write.csv(bayes_results, gzfile(paste0("out/modelfits/", modelname, "-results.csv.gz")))
  rm(bayes_results)

  sum_table <- summary(fit, pars=parameters)$summary
  write.csv(sum_table, paste0("out/modelfits/", modelname, "-summary.csv"))

  for (para in parameters){
    tp1 <- traceplot(fit, pars=c(para))
    ggsave(filename = paste0("out/modelfits/", modelname, "-trace-", para, ".pdf"), plot = tp1)
  }

  print(fit, par=parameters)
  return(fit)
}



# Load Data --------------------------------------------------------------------
Sample <- read_dta("out/data/final-scaled-sample.dta")
Sample$Int <- 1
summary(Sample)
# gvkey               fyear          sic                  TA
# Length:139264      Min.   :1989   Length:139264      Min.   :-37.3374
# Class :character   1st Qu.:1997   Class :character   1st Qu.: -0.3217
# Mode  :character   Median :2003   Mode  :character   Median :  0.1124
# Mean   :2003                      Mean   :  0.0000
# 3rd Qu.:2010                      3rd Qu.:  0.4349
# Max.   :2017                      Max.   : 78.3030
# ChCRev              PPE               InAt              RoA
# Min.   :-24.0143   Min.   :-1.3428   Min.   :-0.6422   Min.   :-19.8465
# 1st Qu.: -0.3927   1st Qu.:-0.8081   1st Qu.:-0.6096   1st Qu.: -0.2202
# Median : -0.1469   Median :-0.2170   Median :-0.4649   Median :  0.1667
# Mean   :  0.0000   Mean   : 0.0000   Mean   : 0.0000   Mean   :  0.0000
# 3rd Qu.:  0.3026   3rd Qu.: 0.6514   3rd Qu.: 0.1414   3rd Qu.:  0.4683
# Max.   : 39.9658   Max.   :11.1814   Max.   : 4.4276   Max.   :  9.8125
# Lev             StdSale            sic2             IndYearID
# Min.   :-1.1252   Min.   :-0.8791   Length:139264      Min.   :   1.0
# 1st Qu.:-0.9509   1st Qu.:-0.6365   Class :character   1st Qu.: 491.0
# Median :-0.1218   Median :-0.3183   Mode  :character   Median : 760.0
# Mean   : 0.0000   Mean   : 0.0000                      Mean   : 833.6
# 3rd Qu.: 0.6191   3rd Qu.: 0.2731                      3rd Qu.:1141.0
# Max.   : 3.5240   Max.   :35.9346                      Max.   :1738.0
# FirmID           Int
# Min.   :    1   Min.   :1
# 1st Qu.: 2603   1st Qu.:1
# Median : 6454   Median :1
# Mean   : 7026   Mean   :1
# 3rd Qu.:11143   3rd Qu.:1
# Max.   :16753   Max.   :1



# Select which model to run ----------------------------------------------------



# Industry-year-varying coefficients model -------------------------------------
# The
# Bayesian version of the most common approach that fits a regression based on industry
# year subsamples Rather than doing that we use a hierarchical model To pool industry-year
# coefficients and pool accross coefficients
modelname <- "industryyear"
filename <- "stan/industryyear.stan"
model_vars <- c("InAt", "ChCRev", "PPE")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$IndYearID),
                  K = length(model_vars),
                  IndYearID = Sample$IndYearID)
fit <- run_model(modelname, filename, parameters, input_data)



# Extended industry-year-varying coefficients model ----------------------------
modelname <- "industryyear-ext"
filename <- "stan/industryyear-ext.stan"
model_vars <- c("InAt", "ChCRev", "PPE", "RoA", "Lev", "StdSale")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$IndYearID),
                  K = length(model_vars),
                  IndYearID = Sample$IndYearID)
fit <- run_model(modelname, filename, parameters, input_data)



# Industryyear model with intercept instead of InAt ----------------------------
modelname <- "industryyear-int"
filename <- "stan/industryyear-ext.stan"
model_vars <- c("Int", "ChCRev", "PPE", "RoA", "Lev", "StdSale")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$IndYearID),
                  K = length(model_vars),
                  IndYearID = Sample$IndYearID)
fit <- run_model(modelname, filename, parameters, input_data)



# Firm-varying coefficients model ----------------------------------------------
modelname <- "firm"
filename <- "stan/firm.stan"
model_vars <- c("InAt", "ChCRev", "PPE")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$FirmID),
                  K = length(model_vars),
                  FirmID = Sample$FirmID)
fit <- run_model(modelname, filename, parameters, input_data)



# Extended firm-varying coefficients model -------------------------------------
modelname <- "firm-ext"
filename <- "stan/firm-ext.stan"
model_vars <- c("InAt", "ChCRev", "PPE", "RoA", "Lev", "StdSale")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$FirmID),
                  K = length(model_vars),
                  FirmID = Sample$FirmID)
fit <- run_model(modelname, filename, parameters, input_data)



# Firm model with intercept instead of InAt ------------------------------------
modelname <- "firm-int"
filename <- "stan/firm-ext.stan"
model_vars <- c("Int", "ChCRev", "PPE", "RoA", "Lev", "StdSale")
parameters <- c("mu_b", "sigma", "L_Omega", "tau")
input_data = list(TA = Sample$TA,
                  x = Sample[, model_vars],
                  N = nrow(Sample),
                  J = max(Sample$FirmID),
                  K = length(model_vars),
                  FirmID = Sample$FirmID)
fit <- run_model(modelname, filename, parameters, input_data)



# END --------------------------------------------------------------------------
