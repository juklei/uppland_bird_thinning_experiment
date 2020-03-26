## model birds in a hierarchical model
## 
## First edit: 20191201
## Last edit: 20191212
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

## 2. Load and prepare data for analysis ---------------------------------------

bpo <- read.csv("clean/bpo_double.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")
bird_data <- read.csv("data/bird_data.csv")

## Make numeric levels for the treatment and experiment matrices:
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)
forest$treatment_num <- as.numeric(forest$treatment)

## Add numeric of species for use below:
bird_data$numeric <- as.numeric(bird_data$short)

## Create model data set:
data <- list(nobs = nrow(bpo),
             observed = bpo$n_obs,
             nvisits = bpo$n_visits,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             exp = acast(forest, year ~ plot, value.var = "exp_num"),
             treat = acast(forest, plot ~ ., 
                           value.var = "treatment_num", 
                           fun.aggregate = mean),
             eval = c(1, 2, 4), ## Check: levels(forest$treatment)
             ref = 3, ## Check: levels(forest$treatment)
             insect = bird_data$numeric[bird_data$food == "insectivore"],
             omni = bird_data$numeric[bird_data$food == "omnivore"],
             bark = bird_data$numeric[bird_data$foraging == "bark"],
             f_grd = bird_data$numeric[bird_data$foraging == "ground"],
             f_cpy = bird_data$numeric[bird_data$foraging == "canopy"],
             grd_cpy = bird_data$numeric[bird_data$foraging == "ground/canopy"],
             n_grd = bird_data$numeric[bird_data$nesting == "ground"],
             n_cpy = bird_data$numeric[bird_data$nesting == "canopy"],
             hole = bird_data$numeric[bird_data$nesting == "hole"]) 

str(data)

## 3. Define and run the model -------------------------------------------------

## For true_occ:
T2 <- array(1, dim = c(max(data$species), max(data$year), max(data$site)))

## Initial values for random slopes:
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, sd_a_pdet = 3,
                    mu_b_pdet_2018 = 0.5, sd_b_pdet_2018 = 5,
                    mu_b_pdet_2019 = 0.5, sd_b_pdet_2019 = 5,
                    mu_obs = 0.5, sd_obs = 5,
                    sd_OLRE = 2,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(2.1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, sd_a_pdet = 2,
                    mu_b_pdet_2018 = 0.7, sd_b_pdet_2018 = 3,
                    mu_b_pdet_2019 = 0.7, sd_b_pdet_2019 = 3,
                    mu_obs = -0.5, sd_obs = 1.5,
                    sd_OLRE = 3,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 1,
                    mu_b_pocc_2019 = 0.8, sd_b_pocc_2019 = 4),
               list(occ_true = T2,
                    mu_a_pdet = 0, sd_a_pdet = 1,
                    mu_b_pdet_2018 = 0.1, sd_b_pdet_2018 = 1,
                    mu_b_pdet_2019 = 0.1, sd_b_pdet_2019 = 1,
                    mu_obs = 1, sd_obs = 1,
                    sd_OLRE = 4,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
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
                   n.adapt = 10000, 
                   inits = inits,
                   n.chains = 3) 

parUpdate(cl = cl, object = "bpo_bin", n.iter = 10000)

samples <- 10000
n.thin <- 10

zc1 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_a_pdet", "sd_a_pdet",
                                         "mu_b_pdet_2018", "sd_b_pdet_2018",
                                         "mu_b_pdet_2019", "sd_b_pdet_2019",
                                         "mu_obs", "sd_obs", 
                                         "sd_OLRE",
                                         "mu_a_pocc", "sd_a_pocc", 
                                         "mu_b_pocc_2018", "sd_b_pocc_2018",
                                         "mu_b_pocc_2019", "sd_b_pocc_2019",
                                         "u_sd_block"),
                      n.iter = samples, thin = n.thin)

zc2 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("a_pdet", "a_pocc"),
                      n.iter = samples, thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/hyperparams.txt")
# capture.output(summary(zc2), HPDinterval(zc2, prob = 0.95)) %>%
#   write(., "results/params.txt")

## 4. Validate the model and export validation data and figures ----------------

pdf("figures/plot_hparams.pdf"); plot(zc1); gelman.plot(zc1); dev.off()
# pdf("figures/plot_params.pdf"); plot(zc2); gelman.plot(zc2); dev.off()

capture.output(raftery.diag(zc1), 
               heidel.diag(zc1), 
               gelman.diag(zc1),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., "results/diagn_hparams.txt")
# capture.output(raftery.diag(zc2),
#                heidel.diag(zc2),
#                gelman.diag(zc2)) %>% write(., "results/dign_params.txt")

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

## 5. Extract informtion from posterior ----------------------------------------

## Compare true control and treatments before the experiment:
zc4 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("compare_before"),
                      n.iter = samples,
                      thin = n.thin)
## Export differences:
capture.output(summary(zc4), HPDinterval(zc4, prob = 0.95)) %>% 
  write(., "results/difference_before_experiment.txt")

## For community and species level metrics:
zc5 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("CI_div_sl", "CI_ctr_sl", "BACI_sl", 
                                         "CI_div_cm", "CI_ctr_cm", "BACI_cm"),
                      n.iter = samples,
                      thin = n.thin)

## Combine MCMC chains:
zc5 <- combine.mcmc(zc5)

## Extract slopes and add ecdf and species names:
BACI_sl <- as.data.frame(summary(zc5)$quantiles[, c("2.5%","50%","97.5%")])
BACI_sl$ecdf <- as.vector(apply(zc5, 2, function(x) 1-ecdf(x)(0)))
BACI_sl$identity <- c(rep("cm", length(data$eval)),
                      sort(rep(levels(bpo$species), length(data$eval))))

## Add treatment*BACI indicator categorisation:
BACI_sl$treatment <- levels(forest$treatment)[data$eval]
BACI_sl$indicator <- c(rep("BACI", length(data$eval)*(max(data$species)+1)),
                       rep("CI_ctr", length(data$eval)*(max(data$species)+1)),
                       rep("CI_div", length(data$eval)*(max(data$species)+1)))

## Export the data set for figures:
write.csv(BACI_sl, "clean/BACI_sl.csv")

## Calculate BACI indicators for all bird_data:

zc6 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      c("CI_div_bd", "CI_ctr_bd", "BACI_bd",
                        "CI_div_r", "CI_ctr_r", "BACI_r",
                        "CI_div_insect", "CI_ctr_insect", "BACI_insect",
                        "CI_div_omni", "CI_ctr_omni", "BACI_omni",
                        "CI_div_bark", "CI_ctr_bark", "BACI_bark",
                        "CI_div_f_cpy", "CI_ctr_f_cpy", "BACI_f_cpy",
                        "CI_div_f_grd", "CI_ctr_f_grd", "BACI_f_grd",
                        "CI_div_grd_cpy", "CI_ctr_grd_cpy", "BACI_grd_cpy",
                        "CI_div_n_grd", "CI_ctr_n_grd", "BACI_n_grd",
                        "CI_div_n_cpy", "CI_ctr_n_cpy", "BACI_n_cpy",
                        "CI_div_hole", "CI_ctr_hole", "BACI_hole"),
                      n.iter = samples,
                      thin = n.thin)

## Combine MCMC chains:
zc6 <- combine.mcmc(zc6)

## Extract slopes and add ecdf and guild names:
BACI_gl <- as.data.frame(summary(zc6)$quantiles[, c("2.5%","50%","97.5%")])
BACI_gl$ecdf <- as.vector(apply(zc6, 2, function(x) 1-ecdf(x)(0)))
guild_names <- names(data)[(length(data) - 8):length(data)]
BACI_gl$identity <- sort(rep(c(guild_names, "r", "bd"), length(data$eval)))

## Add treatment*BACI indicator categorisation:
BACI_gl$treatment <- levels(forest$treatment)[data$eval]
BACI_gl$indicator <- c(rep("BACI", length(data$eval)*11),
                       rep("CI_ctr", length(data$eval)*11),
                       rep("CI_div", length(data$eval)*11))

## Export the data set for figures:
write.csv(BACI_gl, "clean/BACI_gl.csv")

end <- Sys.time()
end - start

## -------------------------------END-------------------------------------------
