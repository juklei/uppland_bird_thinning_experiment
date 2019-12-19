## model birds in a hierarchical model
## 
## First edit: 20191218
## Last edit: 20191218
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

## 2. Load and prepare data ----------------------------------------------------

bpo <- read.csv("clean/bpo_double.csv")
forest <- read.csv("clean/forest_experiment_data.csv")

head(forest)

## Reduce the forest data to plots which are unaffected by the experiment:
forest <- forest[forest$experiment != "difference", ]
forest <- forest[forest$experiment == "before", c(1,2,14,15,23:27)]

## Check correlations:
cor(forest[, 3:length(forest)])

## Check order of plots in forest so its alphabetical and the same as in bpo:
order(forest$plot)

## Create model data set:
data <- list(nobs = nrow(bpo),
             nspecies = max(as.numeric(bpo$species)),
             nyears = 3,
             nsites = nlevels(bpo$plot),
             nblocks = nlevels(bpo$block),
             observed = bpo$n_obs,
             nvisits = bpo$n_visits,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             block = as.numeric(bpo$block),
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             umbr = scale(forest$nr_skarm),
             BA = scale(forest$BA),
             BA_spruce = scale(forest$BA_gran),
             BA_dec = scale(forest$BA_lov),
             BA_dw = scale(forest$BA_dv),
             vis = scale(forest$laser_mean)) 

str(data)

## 3. Define and run the model -------------------------------------------------

## For true_occ:
T2 <- array(1, dim = c(data$nspecies, data$nyears, data$nsites))

## Initial values for random slopes:
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, sd_a_pdet = 3,
                    mu_b_pdet_2018 = 0.5, sd_b_pdet_2018 = 5,
                    mu_b_pdet_2019 = 0.5, sd_b_pdet_2019 = 5,
                    mu_obs = 0.5, sd_obs = 5,
                    sd_OLRE = 2,
                    mu_a_pocc = 0.5, sd_a_pocc = 5,
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5,
                    u_sd_block = 2,
                    mu_b_umbr = 0.5, sd_b_umbr = 5,
                    mu_b_BA = 0.5, sd_b_BA = 5,
                    mu_b_BA_spr = 0.5, sd_b_BA_spr = 5,
                    mu_b_BA_dec = 0.5, sd_b_BA_dec = 5,
                    mu_b_BA_dw = 0.5, sd_b_BA_dw = 4, 
                    mu_b_vis = 0.3, sd_b_vis = 1),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, sd_a_pdet = 2,
                    mu_b_pdet_2018 = 0.7, sd_b_pdet_2018 = 3,
                    mu_b_pdet_2019 = 0.7, sd_b_pdet_2019 = 3,
                    mu_obs = -0.5, sd_obs = 1.5,
                    sd_OLRE = 3,
                    mu_a_pocc = 0.1, sd_a_pocc = 1,
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 1,
                    mu_b_pocc_2019 = 0.8, sd_b_pocc_2019 = 4,
                    u_sd_block = 1,
                    mu_b_umbr = 0, sd_b_umbr = 10,
                    mu_b_BA = 0, sd_b_BA = 10,
                    mu_b_BA_spr = 0, sd_b_BA_spr = 10,
                    mu_b_BA_dec = 0, sd_b_BA_dec = 10,
                    mu_b_BA_dw = 0, sd_b_BA_dw = 8, 
                    mu_b_vis = 0, sd_b_vis = 8),
               list(occ_true = T2,
                    mu_a_pdet = 0, sd_a_pdet = 1,
                    mu_b_pdet_2018 = 0.1, sd_b_pdet_2018 = 1,
                    mu_b_pdet_2019 = 0.1, sd_b_pdet_2019 = 1,
                    mu_obs = 1, sd_obs = 1,
                    sd_OLRE = 0.5,
                    mu_a_pocc = 4, sd_a_pocc = 2,
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 3,
                    mu_b_pocc_2019 = 0.1, sd_b_pocc_2019 = 2,
                    u_sd_block = 3,
                    mu_b_umbr = -0.5, sd_b_umbr = 0.5,
                    mu_b_BA = -0.5, sd_b_BA = 0.5,
                    mu_b_BA_spr = -0.5, sd_b_BA_spr = 0.5,
                    mu_b_BA_dec = -0.5, sd_b_BA_dec = 0.5,
                    mu_b_BA_dw = -0.5, sd_b_BA_dw = 0.4, 
                    mu_b_vis = -0.3, sd_b_vis = 0.4)
               )

model <- "scripts/JAGS/bird_JAGS_bpo_bin_non_exp.R"

start <- Sys.time()

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

jm <- parJagsModel(cl = cl, 
                   name = "bpo_bin",
                   file = model,
                   data = data,
                   n.adapt = 500, 
                   inits = inits,
                   n.chains = 3) 

parUpdate(cl = cl, object = "bpo_bin", n.iter = 500)

samples <- 500
n.thin <- 1

zc1 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_a_pdet", "sd_a_pdet",
                                         "mu_b_pdet_2018", "sd_b_pdet_2018",
                                         "mu_b_pdet_2019", "sd_b_pdet_2019",
                                         "mu_obs", "sd_obs", 
                                         "sd_OLRE",
                                         "mu_a_pocc", "sd_a_pocc", 
                                         "mu_b_pocc_2018", "sd_b_pocc_2018",
                                         "mu_b_pocc_2019", "sd_b_pocc_2019",
                                         "u_sd_block",
                                         "mu_b_umbr", "sd_b_umbr",
                                         "mu_b_BA", "sd_b_BA",
                                         "mu_b_BA_spr", "sd_b_BA_spr",
                                         "mu_b_BA_dec", "sd_b_BA_dec", 
                                         "mu_b_BA_dw", "sd_b_BA_dw", 
                                         "mu_b_vis", "sd_b_vis"),
                      n.iter = samples, thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/params_non_exp_hp.txt")

## 4. Validate the model and export validation data and figures ----------------

pdf("figures/non_exp_hp.pdf"); plot(zc1); gelman.plot(zc1); dev.off()

capture.output(raftery.diag(zc1), 
               heidel.diag(zc1), 
               gelman.diag(zc1),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., "results/dign_non_exp_hp.txt")

# ## Produce validation metrics:
# zc2 <- parCodaSamples(cl = cl, model = "bpo_bin",
#                       variable.names = c("mean_obs", "mean_sim", "p_mean",
#                                          "cv_obs", "cv_sim", "p_cv",
#                                          "fit", "fit_sim", "p_fit"),
#                       n.iter = samples,
#                       thin = n.thin)
# 
# zj_val <- data.frame(combine.mcmc(zc2))
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

## 5. Extract informtion from posterior ----------------------------------------

## Predict from posterior:

zc3 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_b_umbr", "b_umbr", 
                                         "mu_b_BA", "b_BA", 
                                         "mu_b_BA_spr", "b_BA_spr", 
                                         "mu_b_BA_dec", "b_BA_dec", 
                                         "mu_b_BA_dw", "b_BA_dw",
                                         "mu_b_vis", "b_vis"),
                      n.iter = samples,
                      thin = n.thin)

## Combine MCMC chains:
zc3 <- combine.mcmc(zc3)

## Extract slopes and add ecdf and species names:
non_exp <- as.data.frame(summary(zc3)$quantiles[, c("2.5%","50%","97.5%")])
non_exp$ecdf <- as.vector(apply(zc3, 2, function(x) 1-ecdf(x)(0)))
non_exp$identity <- c(rep(levels(bpo$species), 
                          length(zc3[1,])/(data$nspecies+1)),
                      rep("cm", length(zc3[1,])/(data$nspecies+1)))
non_exp$variable <- c(sort(rep(names(data)[13:length(data)], data$nspecies)), 
                      sort(names(data)[13:length(data)]))

## Export the data set for figures:
write.csv(non_exp, "clean/non_exp.csv")

end <- Sys.time()
end - start

## -------------------------------END-------------------------------------------
