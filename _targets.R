# _targets.R file
library(targets)
tt <- tar_target
library(future)
library(future.callr)
plan(callr)

# global packages
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
glu <- glue::glue
# importing source code for the paper
source("R/functions.R")


# Definition of the pipeline
list(
  # Data generation targets -----
  tt(full_sample,
     gen_sample(begin_date = "1988-01-01",
                end_date = "2018-01-01",
                trunc_lvl = 0.01),
     packages = c("RPostgres", "tidyverse")
  ),
  # Running models --------------
  tt(AQwo,
     run_model(full_sample, modelname = "AQwo")
  ),
  tt(Fi_AQwo,
     run_model(full_sample, modelname = "Fi-AQwo")
  ),
  tt(FrankelSun,
     run_model(full_sample, modelname = "FrankelSun")
  ),
  tt(Fi_FrankelSun,
     run_model(full_sample, modelname = "Fi-FrankelSun")
  ),
  tt(ModJones,
     run_model(full_sample, modelname = "ModJones")
  ),
  tt(Fi_ModJones,
     run_model(full_sample, modelname = "Fi-ModJones")
  ),
  tt(ModJonesROA,
     run_model(full_sample, modelname = "ModJonesROA")
  ),
  tt(Fi_ModJonesROA,
     run_model(full_sample, modelname = "Fi-ModJonesROA")
  ),
  tt(ModJonesAQwo,
     run_model(full_sample, modelname = "ModJonesAQwo")
  ),
  tt(Fi_ModJonesAQwo,
     run_model(full_sample, modelname = "Fi-ModJonesAQwo")
  ),
  tt(Collins,
     run_model(full_sample, modelname = "Collins")
  ),
  tt(Fi_Collins,
     run_model(full_sample, modelname = "Fi-Collins")
  ),
  tt(BallShiva,
     run_model(full_sample, modelname = "BallShiva")
  ),
  tt(Fi_BallShiva,
     run_model(full_sample, modelname = "Fi-BallShiva")
  ),
  tt(loo_AQwo,
     compute_loo_diag(AQwo)
  ),
  tt(loo_Fi_AQwo,
     compute_loo_diag(Fi_AQwo)
  ),
  tt(loo_FrankelSun,
     compute_loo_diag(FrankelSun)
  ),
  tt(loo_Fi_FrankelSun,
     compute_loo_diag(Fi_FrankelSun)
  ),
  tt(loo_ModJones,
     compute_loo_diag(ModJones)
  ),
  tt(loo_Fi_ModJones,
     compute_loo_diag(Fi_ModJones)
  ),
  tt(loo_ModJonesROA,
     compute_loo_diag(ModJonesROA)
  ),
  tt(loo_Fi_ModJonesROA,
     compute_loo_diag(Fi_ModJonesROA)
  ),
  tt(loo_ModJonesAQwo,
     compute_loo_diag(ModJonesAQwo)
  ),
  tt(loo_Fi_ModJonesAQwo,
     compute_loo_diag(Fi_ModJonesAQwo)
  ),
  tt(loo_Collins,
     compute_loo_diag(Collins)
  ),
  tt(loo_Fi_Collins,
     compute_loo_diag(Fi_Collins)
  ),
  tt(loo_BallShiva,
     compute_loo_diag(BallShiva)
  ),
  tt(loo_Fi_BallShiva,
     compute_loo_diag(Fi_BallShiva)
  ),
  tt(loo_weights,
     gen_stacking_weights(list(
       "AQwo" = loo_AQwo,
       "Fi-AQwo" = loo_Fi_AQwo,
       "FrankelSun" = loo_FrankelSun,
       "Fi-FrankelSun" = loo_Fi_FrankelSun,
       "ModJones" = loo_ModJones,
       "Fi-ModJones" = loo_Fi_ModJones,
       "ModJonesROA" = loo_ModJonesROA,
       "Fi-ModJonesROA" = loo_Fi_ModJonesROA,
       "ModJonesAQwo" = loo_ModJonesAQwo,
       "Fi-ModJonesAQwo" = loo_Fi_ModJonesAQwo,
       "Collins" = loo_Collins,
       "Fi-Collins" = loo_Fi_Collins,
       "BallShiva" = loo_BallShiva,
       "Fi-BallShiva" = loo_Fi_BallShiva
     )),
     packages = c("loo")
  ),
  tt(averaged_bayes_oos_preds,
     gen_averaged_oos_preds(weight_table = loo_weights, full_sample)
  ),
  tt(save_ndas,
     write_output(bayes_nda = averaged_bayes_oos_preds,
                  fit_sample = full_sample,
                  out_path = "out/bayes-ndas.csv")
  )
)
