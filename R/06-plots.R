# Code for replicating shrinkage figure and stacking figure from
# 'Accounting for Uncertainty'
# authors: Mattias Breuer and Harm H. Schuett
# date:    2019-06-30
# WARNING: Running this script sequentially will take a long time (serveral days).
# IMPORTANT! All scripts are set up to use with a Rstudio project in the root folder of
# this repo. If you do not want to use an rstudio project, you need to manually set the
# working directory to the root folder of the project using setwd()

# Setup -----------------------------------------------------------------------------
library(matrixStats)
library(tidyverse)
library(haven)
library(readr)
library(viridis)
grid.arrange <- gridExtra::grid.arrange
arrangeGrob <- gridExtra::arrangeGrob

extract_gglegend <- function(a.gplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Plot Shrinkage Effects -------------------------------------------------------
make_coef_plots <- function(coef_var,
                            dta,
                            alpha_level=0.2,
                            padding=.3,
                            min_group_size=10,
                            xtitle="Group sample size"
) {
  fit_results <-
    dta %>%
    filter(NrObs > min_group_size) %>%
    filter(term == coef_var)
  highest_sd <- max(c(fit_results$std.error, fit_results$PostSd), na.rm = T)
  ymax <- max(c(fit_results$estimate, fit_results$PostMean), na.rm = T) + padding
  ymin <- min(c(fit_results$estimate, fit_results$PostMean), na.rm = T) - padding
  base_plot <- ggplot(data=fit_results, aes(x=NrObs)) +
    labs(x=xtitle,
         y=paste(coef_var, "point estimate")) +
    scale_y_continuous(limits = c(ymin, ymax)) +
    geom_vline(xintercept = min_group_size, color="grey50")
  p1 <- base_plot +
    geom_point(aes(y=estimate), alpha=alpha_level, size=1) +
    labs(subtitle="Standard OLS approach")
  p2 <- base_plot +
    geom_point(aes(y=PostMean), alpha=alpha_level, size=1)+
    labs(subtitle="Hierarchical Bayesian model")
  return(list(p1, p2))
}


# Industry Coefficient data for shrinkage plots
ols_ind_coefs <- read_csv("out/modelfits/OLS-industry-fits.csv")
bayes_ind_coefs <- read_csv("out/modelfits/industryyear-coef.csv.gz") %>%
  select(-X1) %>%
  as.matrix()

indyr_post_moms <-
  tibble(coeffs = colnames(bayes_ind_coefs),
         PostMean = colMeans(bayes_ind_coefs),
         PostSd = colSds(bayes_ind_coefs)) %>%
  mutate(IndYearID = as.integer(str_extract(coeffs, "(?<=\\[)\\d{1,4}")),
         CoeffID = str_extract(coeffs, "(?<=,)\\d{1}")) %>%
  mutate(term = case_when(CoeffID == 1 ~ "InAt",
                          CoeffID == 2 ~ "ChCRev",
                          CoeffID == 3 ~ "PPE"
  )) %>%
  select(-coeffs, -CoeffID)

indyr_coefs_all <-
  ols_ind_coefs %>%
  rename(NrObs = IndRegNrObs) %>%
  left_join(indyr_post_moms, by=c("IndYearID", "term"))

indu_title <- "Industry-year group sample size"
plot_list <- list()
plot_list <- c(plot_list, make_coef_plots("InAt", indyr_coefs_all, xtitle=indu_title))
plot_list <- c(plot_list, make_coef_plots("ChCRev", indyr_coefs_all, xtitle=indu_title))
plot_list <- c(plot_list, make_coef_plots("PPE", indyr_coefs_all, xtitle=indu_title))
fig1 <- grid.arrange(grobs=plot_list, ncol=2)

ggsave("out/figs/fig-coef-shrinkage-industry.pdf", fig1,
       width=16, height=18, units="cm")

rm(ols_ind_coefs, bayes_ind_coefs, indu_title, plot_list, indyr_coefs_all,
   indyr_post_moms, fig1, make_coef_plots)



# Model Stacking Figure --------------------------------------------------------
compsample <- read_dta("out/data/")

weight_table <- read_csv("out/model-averaging-weights.csv")

# Using the Chesapeake, FYE 1999-12-31 case
firm_year_case <- 8195
short_names <- c("IndYr", "IndYrE", "IndYrI",
                 "Firm", "FirmE", "FirmI")

averaged_draws <- rep.int(0, 2000)
post_firmcase <- NULL
for (i in 1:nrow(weight_table)) {
  filename <- weight_table$model[i]
  print(filename)
  postpredicts <- readRDS(paste0("out/modelfits/", filename, "-ta-coef-only-preds.rds"))
  post_firmcase <- rbind(post_firmcase,
                         tibble(Model = filename,
                                PostDraws = postpredicts[, firm_year_case]))
  averaged_draws <- averaged_draws + postpredicts[, firm_year_case] * weight_table$StackWeight[i]
  rm(postpredicts)
}

post_firmcase2 <-
  post_firmcase %>%
  left_join(weight_table, by=c("Model"="model")) %>%
  mutate(Model = case_when(
    Model == "industryyear-ext" ~ "Model 2",
    Model == "industryyear" ~ "Model 1",
    Model == "industryyear-int" ~ "Model 3",
    Model == "firm" ~ "Model 4",
    Model == "firm-ext" ~ "Model 5",
    Model == "firm-int" ~ "Model 6"
  ),
  Panel = "Left")
post_firmcase2 <- rbind(post_firmcase2,
                        tibble(Model = "Stacked Posteriors",
                               PostDraws=averaged_draws,
                               Panel= "Right",
                               StackWeight=1)
                        )
model_colors <- viridisLite::viridis(6)
fig_modave <- ggplot(post_firmcase2,
                     aes(x=PostDraws, group=Model, fill=Model, color=Model)) +
  geom_density(alpha=0.3) +
  facet_wrap(~Panel, labeller=as_labeller(c("Left"="Indivdual model p(TA |M, D)",
                                            "Right"="Stacked p(TA | D)"))) +
  scale_x_continuous(expand=c(.1, .1)) + #limits=c(-0.1, 1.1)
  labs(y= "Posterior density",
       x= "Level of accruals for Chesapeake, FYE 1999-12-31 (std)",
       tag="B") +
  scale_color_manual(name="",
                     values=c(model_colors, "black"),
                     aesthetics=c("colour", "fill"))+
  theme(legend.position="bottom")

fig_leg <- extract_gglegend(fig_modave)

fig_bar_weights <- ggplot(weight_table %>%
                            mutate(Model = case_when(
                              model == "industryyear-ext" ~ "Model 2",
                              model == "industryyear" ~ "Model 1",
                              model == "industryyear-int" ~ "Model 3",
                              model == "firm" ~ "Model 4",
                              model == "firm-ext" ~ "Model 5",
                              model == "firm-int" ~ "Model 6"
                            )),
                          aes(x=Model, y=StackWeight)) +
  geom_bar(stat="identity",
           fill=model_colors,
           color=model_colors,
           alpha=0.5) +
  labs(y = expression("Stacking weight p(M |D)"),
       x = "",
       tag = "A") +
  theme(axis.title.x = element_blank(),
        legend.position="bottom")

fig_stacking_illu <- grid.arrange(fig_bar_weights,
                                  fig_modave + theme(legend.position="None"),
                                  fig_leg,
                                  layout_matrix=matrix(c(1,1,1,
                                                         2,2,2,3), ncol=1))
ggsave("out/figs/fig-stacking-ppd.pdf", fig_stacking_illu,
       width=14, height=15, units="cm")

rm(fig_bar_weights, fig_modave, fig_stacking_illu, post_firmcase,
   post_firmcase2, weight_table, averaged_draws, filename,
   firm_year_case, i, model_colors, short_names)


