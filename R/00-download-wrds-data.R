# Code for downloading raw data from WRDS
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Imports ---------------------------------------------------------------------------
library(RPostgres)
library(dplyr)
library(lubridate)
library(haven)

# Sample data range
begin_date <- "1988-01-01"
end_date   <- "2018-01-01"

# Establishing wrds connection
# WRDS does not allow your passwords to be saved anywhere in clear form, so
# we prompt for it:
user <- rstudioapi::askForPassword(prompt="Enter your WRDS login")
pass <- rstudioapi::askForPassword()
wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  user=user,
                  password=pass,
                  sslmode='require',
                  dbname='wrds')
wrds  # checking if connection exists



# Downloading Compustat Data --------------------------------------------------------
res1 <- dbSendQuery(wrds,
                    paste("SELECT sich, a.gvkey, a.indfmt,
                                  tic, datadate, fyear, a.conm, fyr, at, che, oiadp, act,
                                  lct, dlc, txp, dp, sale, oancf, ib, rect, ppent, ppegt,
                                  b.loc, b.sic, b.spcindcd, b.ipodate, b.stko, dltt
                           FROM COMP.FUNDA as a
                              LEFT JOIN comp.company as b
                                ON a.gvkey = b.gvkey
                           WHERE consol='C'
                             and indfmt='INDL'
                             and datafmt='STD'
                             and popsrc='D'
                             and datadate between '", begin_date, "' and '", end_date,"'
                           ORDER BY a.gvkey, datadate"))
df_funda <- dbFetch(res1, n=-1)
dbClearResult(res1)


glimpse(df_funda)
df_funda %>% write_dta("data/raw-compustat.dta")
