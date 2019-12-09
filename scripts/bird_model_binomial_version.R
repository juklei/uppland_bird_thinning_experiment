## model birds in a hierarchical model
## 
## First edit: 20191201
## Last edit: 20191201
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(magrittr)
library(reshape2)
library(parallel)
library(dclone)

## 2. Define or source functions used in this script ---------------------------

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

ldm <- read.csv("data/long_distance_migrants.csv")
bpo <- read.csv("clean/bpo_double.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")
head(bpo)
str(bpo)
head(forest)

## 4. The model ----------------------------------------------------------------

# ## Create matrix with seen at least once per year:
# T1 <- acast(bpo[, c("species", "n_obs", "obs_year")], 
#             species ~ obs_year, 
#             value.var = "n_obs", 
#             fun.aggregate = function(x) sum(x) > 0)
# 
# ## Exclude species which where not seen at least once per year:
# red_names <- names(which(rowSums(T1) == 3))
# bpo <- droplevels(bpo[bpo$species %in% red_names, ])

# ## Only for one observer:
# bpo <- droplevels(bpo[bpo$observer == "jkn", ])
# forest <- droplevels(forest[forest$plot %in% bpo$plot, ])

## Create treatment variable comparing true controls with intervention sites:
forest$t_TCvsNC <- ifelse(forest$treatment == "TC", "TC", "NC")

## Create experiment variable where true control is not part of the experimental 
## evaluation:
forest$e_within <- forest$experiment
forest$e_within[forest$treatment == "TC"] <- "before"

## Make numeric levels for the treatment*experiment matrix:
forest$treat_num <- as.numeric(forest$treatment)
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)
forest$t_TCvsNC_num <- as.numeric(as.factor(forest$t_TCvsNC))
forest$e_within_num <- ifelse(forest$e_within == "before", 1, 2)

## Create data arrays for process model part:
treat <- acast(unique(forest[, c("plot", "treat_num")]), . ~ plot)[1, ]
exp <- acast(forest[, c("plot", "year", "exp_num")], year ~ plot)
t_TCvsNC <- acast(unique(forest[, c("plot", "t_TCvsNC_num")]), . ~ plot)[1, ]
e_within <- acast(forest[, c("plot", "year", "e_within_num")], year ~ plot)

## Create model data set:
data <- list(nobs = nrow(bpo),
             nspecies = max(as.numeric(bpo$species)),
             nyears = dim(exp)[1],
             nsites = dim(exp)[2],
             observed = bpo$n_obs,
             nvisits = bpo$n_visits,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             treat = treat, 
             exp = exp,
             t_TCvsNC = t_TCvsNC,
             e_within = e_within) 

str(data)

## Define initial values:

## For true_occ:
T2 <- array(1, dim = c(data$nspecies, data$nyears, data$nsites))

## Initial values for random slopes:
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, 
                    # u_sd_a_pdet = 0.5,
                    mu_b_pdet_2018 = 0.5, sd_b_pdet_2018 = 5,
                    mu_b_pdet_2019 = 0.5, sd_b_pdet_2019 = 5,
                    # b_observer = 0,
                    mu_obs = 0.5, sd_obs = 5,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(2.1, max(data$treat), max(data$exp)),
                    # u_sd_year = 2,
                    # u_sd_site = 2,
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5),
               list(occ_true = T2,
                    # mu_a_pdet = -0.5, 
                    mu_b_pdet_2018 = 0.7, sd_b_pdet_2018 = 3,
                    mu_b_pdet_2019 = 0.7, sd_b_pdet_2019 = 3,
                    u_sd_a_pdet = 3,
                    # b_observer = 0,
                    mu_obs = -0.5, sd_obs = 1.5,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    # u_sd_year = 1,
                    # u_sd_site = 1,
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 1,
                    mu_b_pocc_2019 = 0.8, sd_b_pocc_2019 = 4),
               list(occ_true = T2,
                    mu_a_pdet = 0, 
                    # u_sd_a_pdet = 1,
                    mu_b_pdet_2018 = 0.1, sd_b_pdet_2018 = 1,
                    mu_b_pdet_2019 = 0.1, sd_b_pdet_2019 = 1,
                    # b_observer = 0,
                    mu_obs = 1, sd_obs = 1,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
                    # u_sd_year = 0.5,
                    # u_sd_site = 0.5,
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 3,
                    mu_b_pocc_2019 = 0.1, sd_b_pocc_2019 = 2)
               )

model <- "scripts/JAGS/bird_JAGS_bpo_bin.R"

start <- Sys.time()

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

jm <- parJagsModel(cl = cl, 
                   name = "bpo_bin",
                   file = model,
                   data = data,
                   n.adapt = 5000, 
                   inits = inits,
                   n.chains = 3) 

parUpdate(cl = cl, object = "bpo_bin", n.iter = 40000)

samples <- 50000
n.thin <- 100

zc1 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_a_pdet", "sd_a_pdet", 
                                         "u_sd_pdet_year",
                                         "mu_b_pdet_2018", "sd_b_pdet_2018",
                                         "mu_b_pdet_2019", "sd_b_pdet_2019",
                                         "mu_obs", "sd_obs", 
                                         "sd_OLRE",
                                         "mu_a_pocc", "sd_a_pocc", 
                                         "mu_b_pocc_2018", "sd_b_pocc_2018",
                                         "mu_b_pocc_2019", "sd_b_pocc_2019",
                                         "u_sd_year", "u_sd_site"),
                      n.iter = samples, thin = n.thin)

zc2 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("a_pdet", "sd_pdet_year", "b_observer",
                                         "a_pocc", "sd_year", "sd_site"),
                      n.iter = samples, thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_hyperparams.txt")
capture.output(summary(zc2), HPDinterval(zc2, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_params.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_binomial_hp.pdf"); plot(zc1); gelman.plot(zc1); dev.off()
pdf("figures/plot_binomial_p.pdf"); plot(zc2); gelman.plot(zc1); dev.off()

capture.output(raftery.diag(zc1), 
               heidel.diag(zc1), 
               gelman.diag(zc1),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., "results/diagnostics_binomial_hyperparams.txt")
capture.output(raftery.diag(zc2), 
               heidel.diag(zc2), 
               gelman.diag(zc2),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., "results/diagnostics_binomial_params.txt")

# ## Produce validation metrics:
# zc3 <- parCodaSamples(cl = cl, model = "bpo_bin",
#                       variable.names = c("mean_obs", "mean_sim", "p_mean",
#                                          "cv_obs", "cv_sim", "p_cv",
#                                          "fit", "fit_sim", "p_fit"),
#                       n.iter = samples,
#                       thin = n.thin)
# 
# zj_val <- data.frame(combine.mcmc(zc3))
# 
# ## Fit of mean:
# plot(zj_val$mean_obs, zj_val$mean_sim,
#      xlab = "mean real", ylab = "mean simulated",
#      cex = .05)
# abline(0, 1)
# text(x = 0.6, y = 0.77, paste0("P=", round(mean(zj_val$p_mean), 3)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_obs, zj_val$cv_sim, 
#      xlab = "cv real", ylab = "cv simulated",
#      cex = .05)
# abline(0,1)
# text(x = 1.5, y = 2.20, paste0("P=", round(mean(zj_val$p_cv), 3)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit, zj_val$fit_sim,
#      xlab = "ssq real", ylab = "ssq simulated",
#      cex = .05)
# abline(0,1)
# text(x = 2700, y = 3200, paste0("P=", round(mean(zj_val$p_fit), 3)), cex = 1.5)

## 6. Extract informtion from posterior ----------------------------------------

## For community metrics:

zc4 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("CI_div_C_r", "CI_ctr_C_r", "BACI_C_r", 
                                         "CI_div_T_r", "CI_ctr_T_r", "BACI_T_r",
                                         "CI_div_URT_r", "CI_ctr_URT_r", "BACI_URT_r",
                                         "CI_div_C_bd", "CI_ctr_C_bd", "BACI_C_bd", 
                                         "CI_div_T_bd", "CI_ctr_T_bd", "BACI_T_bd",
                                         "CI_div_URT_bd", "CI_ctr_URT_bd", "BACI_URT_bd"),
                      n.iter = samples,
                      thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc4), HPDinterval(zc4, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_community.txt")

## For species level metrics:

zc5 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("CI_div_C", "CI_ctr_C", "BACI_C", 
                                         "CI_div_T", "CI_ctr_T", "BACI_T",
                                         "CI_div_URT", "CI_ctr_URT", "BACI_URT",
                                         "CI_div_C_cm", "CI_ctr_C_cm", "BACI_C_cm", 
                                         "CI_div_T_cm", "CI_ctr_T_cm", "BACI_T_cm",
                                         "CI_div_URT_cm", "CI_ctr_URT_cm", "BACI_URT_cm"),
                      n.iter = samples,
                      thin = n.thin)

## Combine MCMC chains:
zc5 <- combine.mcmc(zc5)

## Extract slopes and add ecdf and species names:
BACI_sl <- as.data.frame(summary(zc5)$quantiles[, c("2.5%","50%","97.5%")])
BACI_sl$ecdf <- as.vector(apply(zc5, 2, function(x) 1-ecdf(x)(0)))
BACI_sl$species <- rep(c(levels(bpo$species), "cm"), ((max(data$treat)-1)*3))

## Add treatment*BACI indicator categorisation:
BACI_sl$treatment <- rep(c(rep("C", (data$nspecies+1)), 
                           rep("T", (data$nspecies+1)),
                           rep("URT", (data$nspecies+1))), 3)
BACI_sl$indicator <- c(rep("BACI", (max(data$treat)-1)*(data$nspecies+1)),
                       rep("CI_ctr", (max(data$treat)-1)*(data$nspecies+1)),
                       rep("CI_div", (max(data$treat)-1)*(data$nspecies+1)))

## Export the data set for figures:
write.csv(BACI_sl, "clean/BACI_sl.csv")

end <- Sys.time()
end - start

## -------------------------------END-------------------------------------------
