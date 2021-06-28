# helper functions -------
slide_dbl <- slider::slide_dbl

demean_var <- function(x) {x - mean(x)}

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



filter_w_record <- function(dta, note = "", record_store, ..., silence = F) {

  dta_name <- quo_name(enquo(dta))
  dta2 <- dta %>% filter(...)
  if (silence == F) {
    new_entry <- data.frame(Dataset     = dta_name,
                            Firm.years  = nrow(dta2),
                            Firms       = length(unique(dta2$gvkey)),
                            Data.filter = note)

    store_name <- quo_name(enquo(record_store))
    if (!exists(store_name) | startsWith(note, "Start")) {
      message("Create storage data.frame '", store_name, "'")
      first_entry <- data.frame(Dataset     = dta_name,
                                Firm.years  = nrow(dta),
                                Firms       = length(unique(dta$gvkey)),
                                Data.filter = "Raw Data")
      assign(store_name, first_entry, envir = .GlobalEnv)
    }

    assign(store_name, rbind({{ record_store }}, new_entry), envir = .GlobalEnv)
  }

  return(dta2)
}




filter_at_w_record <- function(dta, note = "", record_store, ..., silence = F) {

  dta_name <- quo_name(enquo(dta))
  dta2 <- dta %>% filter_at(...)

  if (silence == F) {
    new_entry <- data.frame(Dataset     = dta_name,
                            Firm.years  = nrow(dta2),
                            Firms       = length(unique(dta2$gvkey)),
                            Data.filter = note)

    store_name <- quo_name(enquo(record_store))
    if (!exists(store_name) | startsWith(note, "Start")) {
      message("Create storage data.frame '", store_name, "'")
      first_entry <- data.frame(Dataset     = dta_name,
                                Firm.years  = nrow(dta),
                                Firms       = length(unique(dta$gvkey)),
                                Data.filter = "Raw Data")
      assign(store_name, first_entry, envir = .GlobalEnv)
    }

    assign(store_name, rbind({{ record_store }}, new_entry), envir = .GlobalEnv)
  }

  return(dta2)
}





# target functions ------
gen_sample <- function(begin_date = "1988-01-01",
                       end_date = "2018-01-01",
                       trunc_lvl = 0.01) {




  # see the keyring documentation for storing username and passwords
  # keyring::default_backend()
  # keyring::key_set("wrds_user")
  # keyring::key_set("wrds_pw")
  # keyring::key_list()

  wrds <- dbConnect(Postgres(),
                    host = "wrds-pgdata.wharton.upenn.edu",
                    port = 9737,
                    user = keyring::key_get("wrds_user"),
                    password = keyring::key_get("wrds_pw"),
                    sslmode = "require",
                    dbname = "wrds"
  )

  res3 <- dbSendQuery(wrds,
                      paste("SELECT sich, a.gvkey, a.indfmt,
                                  tic, datadate, fyear, a.conm, fyr, at, che, oiadp, act,
                                  lct, dlc, txp, dp, sale, oancf, ib, rect, ppent, ppegt,
                                  invt, ap, cogs,
                                  b.loc, b.sic, b.spcindcd, b.ipodate, b.stko, dltt,
                                  mkvalt, ceq
                           FROM COMP.FUNDA as a
                              LEFT JOIN comp.company as b
                                ON a.gvkey = b.gvkey
                           WHERE consol='C'
                             and indfmt='INDL'
                             and datafmt='STD'
                             and popsrc='D'
                             and datadate between '", begin_date, "' and '", end_date,"'
                           ORDER BY a.gvkey, datadate"))
  raw_sample <- dbFetch(res3, n=-1)
  dbClearResult(res3)
  dbDisconnect(wrds)

  model_vars <- c("TA", "ChCRev", "PPE", "InAt", "RoA", "Lev", "StdSale",
                  "SG", "CFOm1", "CFOt0", "CFOp1", "ChCFO", "BBtoM",
                  "NPM", "OpCycle")
  print("sample down")
  prepped_raw_sample <-
    raw_sample %>%
    # getting rid of financials
    mutate(sich = as.character(sich)) %>%
    mutate(sic = if_else(is.na(sich) == F, sich, sic)) %>%
    mutate(Sic3 = as.integer(substr(sic, 1, 3))) %>%
    filter(Sic3 <= 600 | Sic3 >= 699) %>%
    # setting debt and receivables to zero if missing
    mutate(dltt = if_else(is.na(dltt), 0 , dltt),
           dlc  = if_else(is.na(dlc), 0 , dlc),
           rect = if_else(is.na(rect), 0 , rect))
  print("first use f tilde")
  computed_var_sample <-
    prepped_raw_sample %>%
    # Computing Regression variables
    mutate(Accrual = ib - oancf) %>%
    arrange(gvkey, datadate) %>%
    group_by(gvkey) %>%
    mutate(
      lag_fyear = lag(fyear),
      lead_fyear = lead(fyear),
      lag_at = lag(at),
      TA = Accrual / lag_at,
      ChCRev = (sale - lag(sale) - (rect - lag(rect))) / lag_at,
      PPE = ppegt / lag_at,
      InAt = 1 / lag_at,
      RoA = oiadp / lag_at,
      Lev = lag((dltt + dlc) / at),
      CFOm1 = lag(oancf) / lag_at,
      CFOt0 = oancf / lag_at,
      CFOp1 = lead(oancf) / lag_at,
      ChCFO = (oancf - lag(oancf)) / lag_at,
      SG = sale / lag(sale),
      BBtoM = lag(ceq / mkvalt),
      DuCFO = if_else(CFOt0 < 0, 1, 0),
      StdSale = roll_rsd(sale)/lag_at,
      # Following based on Frankel Sun 2018, p. 170
      purchases = cogs + invt - lag(invt),
      NPM = slide_dbl(ib / sale, ~mean(., na.rm = T), .before = 2, .complete = F),
      yearly_opc = 0.5 * (rect + lag(rect)) / sale +
        0.5 * (invt + lag(invt)) / cogs -
        0.5 * (ap + lag(ap)) / purchases,
      OpCycle = slide_dbl(yearly_opc, ~mean(., na.rm = T), .before = 2, .complete = F),
      AssetGrowth = (at - lag_at)/lag_at,

      # Validation test vars
      EP = ib / lag(mkvalt),
      Size = log(lag(at)),
      MSize = log(lag(mkvalt))
    ) %>%
    ungroup()

  filtered_var_sample <-
    computed_var_sample %>%
    filter_w_record("Start: Aligned fyears", t1,
                    fyear - 1 == lag_fyear & fyear + 1 == lead_fyear,
    ) %>%
    filter_w_record("lag_at > 10 (Hribar Collins 2002, p.108)", t1,
                    lag_at > 10) %>%
    filter_w_record("No Acc", t1, is.na(TA) == F) %>%
    filter_w_record("No CFO_t", t1, is.na(CFOt0) == F) %>%
    filter_w_record("No ROA", t1, is.na(RoA) == F) %>%
    filter_w_record("No StdSale", t1,
                    is.na(StdSale) == F,
                    is.infinite(StdSale) == F) %>%
    filter_w_record("No Lev", t1,
                    is.na(Lev) == F,
                    is.infinite(Lev) == F) %>%
    filter_w_record("No PPE", t1, is.na(PPE) == F) %>%
    filter_w_record("No ChCRev", t1, is.na(ChCRev) == F) %>%
    filter_w_record("No InAt", t1, is.na(InAt) == F) %>%
    filter_w_record("No SG", t1,
                    is.na(SG) == F,
                    is.infinite(SG) == F) %>%
    filter_w_record("No CFO t-1", t1, is.na(CFOm1) == F) %>%
    filter_w_record("No ChCFO", t1, is.na(ChCFO) == F) %>%
    filter_w_record("No CFO t + 1", t1, is.na(CFOp1) == F) %>%
    filter_w_record("No BtM", t1, is.na(BBtoM) == F) %>%
    filter_w_record("No NPM", t1,
                    is.na(NPM) == F,
                    is.infinite(NPM) == F) %>%
    filter_w_record("No OpCycle", t1,
                    is.na(OpCycle) == F,
                    is.infinite(OpCycle) == F)


  # truncating at the firm level
  truncated_Sample <-
    filtered_var_sample %>%
    # dropping everything but the regression vars
    select(gvkey, fyear, sic, all_of(model_vars), DuCFO, AssetGrowth) %>%
    # de-meaning by firm
    group_by(gvkey) %>%
    mutate_at(
      vars(all_of(model_vars)),
      list(demean=demean_var)
    ) %>%
    ungroup() %>%
    #  truncating extreme values
    mutate(
      TA_demean_trunc = truncate(TA_demean, perc=trunc_lvl),
      ChCRev_demean_trunc = truncate(ChCRev_demean, perc=trunc_lvl),
      PPE_demean_trunc = truncate(PPE_demean, perc=trunc_lvl, bottom_trunc = F),
      InAt_demean_trunc = truncate(InAt_demean, perc=trunc_lvl),
      RoA_demean_trunc = truncate(RoA_demean, perc=trunc_lvl),
      StdSale_demean_trunc = truncate(StdSale_demean, perc=trunc_lvl, bottom_trunc = F),
      SG_demean_trunc = truncate(SG_demean, perc=trunc_lvl, bottom_trunc = F),
      CFOm1_demean_trunc = truncate(CFOm1_demean, perc=trunc_lvl),
      CFOt0_demean_trunc = truncate(CFOt0_demean, perc=trunc_lvl),
      CFOp1_demean_trunc = truncate(CFOp1_demean, perc=trunc_lvl),
      BBtoM_demean_trunc = truncate(BBtoM_demean, perc=trunc_lvl)
    )


  unscaled_Sample <-
    truncated_Sample %>%
    filter_at_w_record("truncate", t1,
                       vars(ends_with("trunc")),
                       ~is.na(.) == F
    ) %>%
    select(gvkey, fyear, sic, all_of(model_vars), DuCFO, AssetGrowth) %>%
    filter_w_record("Leverage <= 1" , t1, Lev <= 1) %>%
    filter_w_record("No M&A (200% assets growth)", t1, AssetGrowth <= 2) %>%
    select(-AssetGrowth) %>%
    group_by(fyear) %>%
    mutate(  # Following based on Frankel Sun 2018, p. 170
      NPM = percent_rank(NPM),
      OpCycle = percent_rank(OpCycle),
      ChOCFSC = 0.5 * (OpCycle + (1 - NPM))
    ) %>%
    ungroup()


  # Scaled Sample with firm and industry indicators
  fin_unscaled_Sample <-
    unscaled_Sample %>%
    mutate(sic2 = str_sub(sic, 1, 2)) %>%
    mutate(IndYearID = as.integer(as.factor(paste0(sic2, "-", fyear))),
           FirmID = as.integer(as.factor(gvkey)))


  Sample <-
    fin_unscaled_Sample %>%
    mutate_at(all_of(c(model_vars, "ChOCFSC")), ~scale(.)[, 1]) %>%
    mutate(
      SGQ = as.factor(ntile(SG, 5)),
      RoAQ = as.factor(ntile(RoA, 5)),
      BBtoMQ = as.factor(ntile(BBtoM, 5))
    )

  dummies <-
    model.matrix(~ RoAQ + SGQ + BBtoMQ, Sample) %>%
    as_tibble() %>%
    select(-`(Intercept)`)

  Sample <-
    Sample %>%
    select(-SGQ, -RoAQ, -BBtoMQ) %>%
    bind_cols(dummies) %>%
    mutate(
      ChCRevxOpCycle = ChCRev * OpCycle,
      ChOCFSCxChCFO = ChOCFSC * ChCFO,
      DuCFOxCFOt0 = DuCFO * CFOt0,
      Int = 1
    )

  scale_moments <-
    map_dfr(c(model_vars, "ChOCFSC"),
            ~unlist(attributes(scale(unscaled_Sample[[.]])))) %>%
    mutate(Var = c(model_vars, "ChOCFSC"))


  full_sample <- list(
    Sample = Sample,
    scale_moments = scale_moments
  )

  write_csv(t1, "out/sample-filter-steps.csv")

  return(full_sample)
}


run_model <- function(fit_sample,
                      modelname,
                      stan_file = "stan/hier-model.stan") {

  message("Starting model: ", modelname)
  options(mc.cores = parallel::detectCores())

  if (modelname %in% c("Fi-Collins", "Fi-BallShiva")) {
    control_ops <- list(max_treedepth = 15)
  } else {
    control_ops <- NULL
  }

  OUT_PATH = "out/"
  parameters <- c("mu_b", "sigma", "L_Omega", "tau")
  model_par_list <- list(
    "ModJones"        = c("InAt", "ChCRev", "PPE"),
    "Fi-ModJones"     = c("InAt", "ChCRev", "PPE"),
    "ModJonesROA"     = c("InAt", "ChCRev", "PPE", "RoA"),
    "Fi-ModJonesROA"  = c("InAt", "ChCRev", "PPE", "RoA"),
    "ModJonesAQwo"    = c("InAt", "ChCRev", "PPE", "CFOm1", "CFOt0"),
    "Fi-ModJonesAQwo" = c("InAt", "ChCRev", "PPE", "CFOm1", "CFOt0"),
    "AQwo"            = c("CFOm1", "CFOt0"),
    "Fi-AQwo"         = c("CFOm1", "CFOt0"),
    "Collins"         = c("Int", "ChCRev",
                          "RoAQ2", "RoAQ3", "RoAQ4", "RoAQ5",
                          "SGQ2", "SGQ3", "SGQ4", "SGQ5",
                          "BBtoMQ2", "BBtoMQ3", "BBtoMQ4", "BBtoMQ5"),
    "Fi-Collins"      = c("Int", "ChCRev",
                          "RoAQ2", "RoAQ3", "RoAQ4", "RoAQ5",
                          "SGQ2", "SGQ3", "SGQ4", "SGQ5",
                          "BBtoMQ2", "BBtoMQ3", "BBtoMQ4", "BBtoMQ5"),
    "BallShiva"       = c("Int", "ChCRev", "PPE", "CFOt0", "DuCFO" ,"DuCFOxCFOt0"),
    "Fi-BallShiva"    = c("Int", "ChCRev", "PPE", "CFOt0", "DuCFO" ,"DuCFOxCFOt0"),
    "FrankelSun"      = c("Int", "ChCRev", "OpCycle", "ChCRevxOpCycle",
                          "ChCFO", "ChOCFSC", "ChOCFSCxChCFO"),
    "Fi-FrankelSun"   = c("Int", "ChCRev", "OpCycle", "ChCRevxOpCycle",
                          "ChCFO", "ChOCFSC", "ChOCFSCxChCFO")
  )

  # Model data
  model_vars <- model_par_list[[modelname]]
  GroupID <- ifelse(startsWith(modelname, "Fi-"), "FirmID", "IndYearID")
  train_sample <- fit_sample[["Sample"]]
  input_data = list(
    TA = train_sample$TA,
    x = train_sample[, model_vars],
    N = nrow(train_sample),
    J = max(train_sample[[GroupID]]),
    K = length(model_vars),
    GroupID = train_sample[[GroupID]]
  )

  # Logging fitting information
  time_stamp <- gsub("[: ]", "-", substr(Sys.time(), 1, 16))
  logfile <- file(glu("{OUT_PATH}{modelname}-{time_stamp}.txt"), open = "wt")
  sink(logfile ,type = "output")
  sink(logfile, type = "message")
  quiet_stan <- purrr::quietly(stan)
  fit_results <- quiet_stan(
    file = stan_file,
    data = input_data,
    # iter = 50,
    # warmup = 30,
    iter = 3500,
    warmup = 1000,
    chains = 2,
    cores = 2,
    seed = 444,
    control = control_ops,
    verbose = FALSE
  )
  fit <- fit_results[["result"]]
  fit_summary <- list(
    time = Sys.time(),
    model_name = modelname,
    group_id = GroupID,
    num_divergent = get_num_divergent(fit),
    num_max_treedepth = get_num_max_treedepth(fit),
    warnings = fit_results[["warnings"]],
    messages = fit_results[["messages"]],
    runtime = get_elapsed_time(fit)
  )
  options(width=200)
  cat("Timestamp", time_stamp, "\n\n")
  print(fit, par = parameters)
  cat("\n\n")
  print(fit_summary)
  sink(type = "message")
  sink()
  close(logfile)


  # Saving results
  bayes_predicted <- as.data.frame(fit, par = "y_fit_wo")
  saveRDS(bayes_predicted, glu("{OUT_PATH}{modelname}-ta-predicions.rds"))
  print(dim(bayes_predicted))
  rm(bayes_predicted)

  loklik <- as.matrix(fit, par = "log_lik")

  post_draws <- rstan::extract(
    fit,
    permuted = FALSE,
    inc_warmup = FALSE,
    pars = c("mu_b", "sigma", "L_Omega", "tau", "b")
  )

  mcmc_diagnostics <-
    rstan::monitor(post_draws, print = FALSE)|>
    as.data.frame() |>
    rownames_to_column(var = "par")
  write_csv(mcmc_diagnostics, glu("{OUT_PATH}{modelname}-mcmcdiag.csv"))

  for (para in parameters){
    tp1 <- traceplot(fit, pars = c(para))
    ggsave(filename = glu("{OUT_PATH}{modelname}-trace-{para}.pdf"), plot = tp1)
  }

  message("Finished model: ", modelname)
  rm(fit, fit_summary, tp1, train_sample, test_sample, GroupID, model_vars, input_data)

  return(loklik)
}




compute_loo_diag <- function(loklik) {
  log_lik <- as.matrix(loklik)
  model_loo <- loo(log_lik)
  return(model_loo)
}




gen_stacking_weights <- function(loo_list) {
  model_names <- names(loo_list)
  wts <- loo_model_weights(loo_list, method = "stacking")
  moweights <- data.frame(model = model_names, StackWeight = as.vector(wts) )

  return(moweights)
}




gen_averaged_oos_preds <- function(weight_table, fit_sample){
  nrows <- nrow(fit_sample[["Sample"]])
  av_TA_co_preds <- matrix(0, nrow=5000, ncol=nrows)

  for (i in 1:nrow(weight_table)) {
    model_name <- weight_table$model[i]
    print(model_name)
    TA_co_preds <- readRDS(glu("out/{model_name}-ta-predicions.rds"))
    print(dim(TA_co_preds))
    TA_co_preds <- TA_co_preds * weight_table$StackWeight[i]
    av_TA_co_preds  <- av_TA_co_preds + TA_co_preds
    remove(TA_co_preds, model_name)
  }

  return(av_TA_co_preds)
}






write_output <- function(bayes_nda, fit_sample, out_path = "out/bayes-ndas.csv"){

  Sample <- fit_sample[["Sample"]]
  scale_moments <- fit_sample[["scale_moments"]]
  mean_TA <- scale_moments[scale_moments$Var == "TA", ]$`scaled:center`
  sd_TA   <- scale_moments[scale_moments$Var == "TA", ]$`scaled:scale`

  bayes_NDA_draws <- as.matrix( (bayes_nda * sd_TA) + mean_TA )
  bayes_NDA <-
    Sample %>%
    select(gvkey, fyear) %>%
    mutate(NDA    = matrixStats::colMeans2(bayes_NDA_draws),
           NDA_sd = matrixStats::colSds(bayes_NDA_draws)
    )

  write_csv(bayes_NDA, out_path)

  return(out_path)
}
