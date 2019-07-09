# Public Repository for "Accounting for Uncertainty: An Application of Bayesian Methods to Accruals Models"

This repository holds code to compute model averaged predictive distributions for accruals. It provides the necessary materials to replicate and extend the approach as proposed in *"Accounting for Uncertainty: An Application of Bayesian Methods to Accruals Models"*.

## How to use this repository

The code is organized into an RStudio project. While nearly all of the code is written in R, the models are coded in [Stan](https://mc-stan.org/) and all MCMC calculations are performed by Stan.

The folder structure is as follows:

- `R/`: Main code. Numerated in the order scripts should be executed.
- `stan/`: Bayesian models used in *"Accounting for Uncertainty"*
- `data/`: [empty] Location for raw data as downloaded from the first script in `R/`
- `out/`: [empty] Storage location for all intermediate results, data, and final figures
- `doc/`: Additional documentation
- `stata`: Examples for Bayesian computations using STATA

*Note: executing the code will result in roughly 42GB of data. The computations are also RAM hungry (32+GB)*

## Requirements

The R package rstan installs a compatible version of Stan automatically. The only prerequisite is Rtools on windows and a compatible compiler tool chain for MacOS and Linux. See the [rstan documentation](https://mc-stan.org/rstan/) for installation instructions.

In addition, the various other packages are used in different scripts:

    library(broom)
    library(dplyr)
    library(gridExtra)
    library(haven)
    library(lfe)
    library(loo)
    library(lubridate)
    library(matrixStats)
    library(purrr)
    library(readr)
    library(RPostgres)
    library(rstan)
    library(stringr)
    library(tidyr)
    library(viridis)

These are the versions used:

    > sessionInfo()
    R version 3.6.0 (2019-04-26)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 18362)

    Matrix products: default

    locale:
    [1] LC_COLLATE=English_United States.1252
    [2] LC_CTYPE=English_United States.1252
    [3] LC_MONETARY=English_United States.1252
    [4] LC_NUMERIC=C
    [5] LC_TIME=English_United States.1252

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets
    [6] methods   base

    other attached packages:
    [1] viridis_0.5.1         viridisLite_0.3.0
    [3] tidyr_0.8.3           stringr_1.4.0
    [5] rstan_2.18.2          StanHeaders_2.18.1-10
    [7] ggplot2_3.2.0         RPostgres_1.1.1
    [9] readr_1.3.1           purrr_0.3.2
    [11] matrixStats_0.54.0    lubridate_1.7.4
    [13] loo_2.1.0             lfe_2.8-3
    [15] Matrix_1.2-17         haven_2.1.0
    [17] gridExtra_2.3         dplyr_0.8.2
    [19] broom_0.5.2

    loaded via a namespace (and not attached):
    [1] zoo_1.8-6         tidyselect_0.2.5  lattice_0.20-38
    [4] colorspace_1.4-1  generics_0.0.2    stats4_3.6.0
    [7] pkgbuild_1.0.3    blob_1.1.1        rlang_0.4.0
    [10] pillar_1.4.2      glue_1.3.1        withr_2.1.2
    [13] DBI_1.0.0         bit64_0.9-7       munsell_0.5.0
    [16] gtable_0.3.0      inline_0.3.15     forcats_0.4.0
    [19] callr_3.2.0       ps_1.3.0          parallel_3.6.0
    [22] Rcpp_1.0.1        xtable_1.8-4      backports_1.1.4
    [25] scales_1.0.0      bit_1.1-14        hms_0.4.2
    [28] stringi_1.4.3     processx_3.3.1    grid_3.6.0
    [31] cli_1.1.0         tools_3.6.0       sandwich_2.5-1
    [34] magrittr_1.5      lazyeval_0.2.2    tibble_2.1.3
    [37] Formula_1.2-3     crayon_1.3.4      pkgconfig_2.0.2
    [40] prettyunits_1.0.2 assertthat_0.2.1  rstudioapi_0.10
    [43] R6_2.4.0          nlme_3.1-140      compiler_3.6.0
    >

Finally, all scripts are set up to use with a RStudio project in the root folder of this repo. If you do not want to use an RStudio project, you need to manually set the working directory of your R session to the root folder of the project using `setwd()`.

## Installing cmdstan for use with STATA

We provide some code to run Bayesian computations using Stata as well. Getting Stata to play nicely with Stan is not completely trivial, but there are bindings called [StataStan](https://mc-stan.org/users/interfaces/stata-stan.html). We provide a short manual for setting up cmdstan for use with Stata in the `doc/` folder. The rest is explained at the StataStan homepage.
