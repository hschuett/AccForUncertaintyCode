# Code for preparing the sample used to fit the accrual models
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports ---------------------------------------------------------------------------
library(dplyr)
library(haven)
library(stringr)


# Definition of helper functions ----------------------------------------------------
demean_var <- function (x) {x - mean(x)}

truncate <- function(variable, perc=0.01, bottom_trunc=T) {
  if (bottom_trunc == TRUE) {
    new_var <- ifelse(variable >= quantile(variable, 1-perc), NA,
                      ifelse(variable <= quantile(variable, perc), NA, variable))
  } else {
    new_var <- ifelse(variable >= quantile(variable, 1-perc), NA,variable)
  }
  return(new_var)
}


roll_rsd <- function(variable) {
  zoo::rollapplyr(data=variable,
                  width=3, FUN=sd, fill=NA, na.rm=T)
}



# Load Data -------------------------------------------------------------------------
raw_sample <- read_dta("data/raw-compustat.dta")
glimpse(raw_sample)



# Construct Final Sample ------------------------------------------------------------
Sample <-
  raw_sample %>%
  # getting rid of financials
  mutate(sich = as.character(sich)) %>%
  mutate(sic = if_else(is.na(sich) == F, sich, sic)) %>%
  mutate(Sic3 = as.integer(substr(sic, 1, 3))) %>%
  filter(Sic3 <= 600 | Sic3 >= 699) %>%
  # setting debt and receivables to zero if missing
  mutate(dltt = if_else(is.na(dltt), 0 , dltt),
         dlc = if_else(is.na(dlc), 0 , dlc),
         rect = if_else(is.na(rect), 0 , rect))

Sample <-
  Sample %>%
  # Computing Regression variables
  mutate(Accrual = ib - oancf) %>%
  arrange(gvkey, datadate) %>%
  group_by(gvkey) %>%
  mutate(TA = Accrual / lag(at),
         ChCRev = (sale - lag(sale) - (rect - lag(rect))) / lag(at),
         PPE = ppegt / lag(at),
         InAt = 1 / lag(at),
         RoA = oiadp / lag(at),
         Lev = lag((dltt + dlc) / at)) %>%
  mutate(lag_at = lag(at)) %>%
  mutate(StdSale = roll_rsd(sale)/lag(at)) %>%
  # keeping those with non-missing
  filter(is.na(lag(at)) == FALSE) %>%
  ungroup() %>%
  filter(complete.cases(TA, ChCRev, PPE, InAt, RoA, Lev, StdSale) & is.infinite(Lev) == F) %>%
  filter(lag_at > 10) %>%
  mutate(AssetGrowth = (at - lag_at)/lag_at)

truncated_Sample <-
  Sample %>%
  # dropping everything but the regression vars
  select(gvkey, fyear, sic, TA, ChCRev, PPE, InAt, RoA, StdSale, Lev, AssetGrowth) %>%
  # de-meaning by firm
  group_by(gvkey) %>%
  mutate_if(is.numeric, funs(demean=demean_var)) %>%
  ungroup() %>%
  #  truncating extreme values
  mutate(TA_demean_trunc = truncate(TA_demean, perc=0.01),
         ChCRev_demean_trunc = truncate(ChCRev_demean, perc=0.01),
         PPE_demean_trunc = truncate(PPE_demean, perc=0.01, bottom_trunc = F),
         InAt_demean_trunc = truncate(InAt_demean, perc=0.01),
         RoA_demean_trunc = truncate(RoA_demean, perc=0.01),
         StdSale_demean_trunc = truncate(StdSale_demean, perc=0.01, bottom_trunc = F))

unscaled_Sample <-
  truncated_Sample %>%
  filter(complete.cases(TA_demean_trunc, ChCRev_demean_trunc, PPE_demean_trunc,
                        InAt_demean_trunc, RoA_demean_trunc, StdSale_demean_trunc, Lev)) %>%
  select(gvkey, fyear, sic, TA, ChCRev, PPE, InAt, RoA, Lev, StdSale, AssetGrowth) %>%
  filter(Lev <= 1)



# Deleting large asset growth firms -------------------------------------
unscaled_Sample <-
  unscaled_Sample %>%
  filter(AssetGrowth <= 2) %>%  # only keep firm-years with less than 200% asset growth
  select(-AssetGrowth)

unscaled_Sample %>% write_dta("out/data/unscaled-sample.dta")



# Scaled Sample with firm and industry indicators ------------------------
Sample2 <-
  unscaled_Sample %>%
  mutate(sic2 = str_sub(sic, 1, 2)) %>%
  mutate(IndYearID = as.integer(as.factor(paste0(sic2, "-", fyear))),
         FirmID = as.integer(as.factor(gvkey)))
Sample2$TA <- scale(Sample2$TA)[,1]
Sample2$PPE <- scale(Sample2$PPE)[,1]
Sample2$ChCRev <- scale(Sample2$ChCRev)[,1]
Sample2$InAt <- scale(Sample2$InAt)[,1]
Sample2$RoA <- scale(Sample2$RoA)[,1]
Sample2$Lev <- scale(Sample2$Lev)[,1]
Sample2$StdSale <- scale(Sample2$StdSale)[,1]

Sample2 %>% write_dta("out/data/final-scaled-sample.dta")

summary(Sample2)
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
# FirmID
# Min.   :    1
# 1st Qu.: 2603
# Median : 6454
# Mean   : 7026
# 3rd Qu.:11143
# Max.   :16753


