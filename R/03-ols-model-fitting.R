# Code for computing standard OLS industry-year modified jones models
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# WARNING: Running this script sequentially will take a long time (serveral days).
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports and Setup -----------------------------------------------------------------
library(matrixStats)
library(tidyr)
library(dplyr)
library(stringr)
library(haven)
library(broom)
library(lfe)
library(purrr)



# Data -------------------------------------------------------------------------
Sample <- read_dta("out/data/final-scaled-sample.dta") %>%
  mutate(RowNr = row_number())

nested_indyr_samples <- Sample %>%
  nest(-sic2, -fyear, -IndYearID)
nested_indyr_samples$IndRegNrObs <- map_int(nested_indyr_samples$data, nrow)
head(nested_indyr_samples)



# Simple Modified Jones Coefs ---------------------------------------------
coefs_simp_jones <- function(data) {
  ols_model <- lm(TA ~ InAt + ChCRev + PPE -1, data=data)
  return(tidy(ols_model))
}
nested_indyr_samples$CoefsS <- map(nested_indyr_samples$data, coefs_simp_jones)
# We have broken models for every sic-year combo with less than 4 obs.
# (And of course, even those with slightly higher n are useless )
# Also, the combo sic2 40 years 2016 and 2017 have rank defficient covariates
# InAt is extremely highly correlated with ChCRev and PPE. This is mainly sic
# "4011" 4011 Railroads, Line-Haul Operating
# Establishments primarily engaged in line-haul railroad passenger and freight operations.
# Railways primarily engaged in furnishing passenger transportation confined principally to a single
# municipality, contiguous municipalities, or a municipality and its suburban areas are classified in Major Group 41.
# So it's small railroads, who probably have as assets only trains.



# Industry-year-model Preds -----------------------------------------------
simp_jones_predictions <- function(data) {
  if (nrow(data) < 10) {  # Making sure there are enough obs. fitted does not produce NAs
    fit_y <- rep.int(NA, nrow(data))
  } else {
    ols_model <- felm(TA ~ InAt + ChCRev + PPE -1 |0|0|gvkey, data=data)
    fit_y <- fitted(ols_model)
  }
  new_data <- data %>%
    mutate(FitIYSimp = fit_y)
  return(new_data)
}

ext_jones_predictions <- function(data) {
  if (nrow(data) < 10) {  # Making sure there are enough obs. fitted does not produce NAs
    fit_y <- rep.int(NA, nrow(data))
  } else {
    ols_model <- felm(TA ~ InAt + ChCRev + PPE + RoA + Lev + StdSale -1 |0|0|gvkey, data=data)
    fit_y <- fitted(ols_model)
  }
  new_data <- data %>%
    mutate(FitIYExt = fit_y)
  return(new_data)
}

nested_indyr_samples$SJonesPred <- map(nested_indyr_samples$data,
                                       simp_jones_predictions)
nested_indyr_samples$EJonesPred <- map(nested_indyr_samples$data,
                                       ext_jones_predictions)

simple_jones_pred <- unnest(nested_indyr_samples, SJonesPred) %>%
  select(RowNr, gvkey, fyear, IndYearID, IndRegNrObs, TA, FitIYSimp)
extend_jones_pred <- unnest(nested_indyr_samples, EJonesPred) %>%
  select(RowNr, gvkey, fyear, IndYearID, IndRegNrObs, TA, FitIYExt)



# Saving everything ------------------------------------------------------------
ols_indu_coefs <- unnest(nested_indyr_samples, CoefsS)
ols_indu_coefs %>% readr::write_csv("out/modelfits/OLS-industry-fits.csv")

ols_predictions <- simple_jones_pred %>%
  left_join(extend_jones_pred[, c("RowNr", "FitIYExt")], by="RowNr") %>%
  arrange(RowNr)

ols_predictions %>% readr::write_csv("out/all-OLS-predictions.csv")
