## model birds in a hierarchical model
## 
## First edit: 20191201
## Last edit: 20191210
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

## 2. Chose treatments and define reference treatment. Then create data. -------

## Chose treatment type you want to evaluate:
TT <- "treatment" 
# TT <- "TCvsNC" 
# TT <- "BA_all"  
# TT <- "BA_dec" 
# TT <- "BA_spruce" 
# TT <- "BA_dw" 
# TT <- "V" 

## Chose reference level:
ref <- "TC"
# ref <- "C"

## Calculate data:
source("scripts/data_calc.r")

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
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(2.1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5,
                    u_sd_block = 2),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, sd_a_pdet = 2,
                    mu_b_pdet_2018 = 0.7, sd_b_pdet_2018 = 3,
                    mu_b_pdet_2019 = 0.7, sd_b_pdet_2019 = 3,
                    mu_obs = -0.5, sd_obs = 1.5,
                    sd_OLRE = 3,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 1,
                    mu_b_pocc_2019 = 0.8, sd_b_pocc_2019 = 4,
                    u_sd_block = 1),
               list(occ_true = T2,
                    mu_a_pdet = 0, sd_a_pdet = 1,
                    mu_b_pdet_2018 = 0.1, sd_b_pdet_2018 = 1,
                    mu_b_pdet_2019 = 0.1, sd_b_pdet_2019 = 1,
                    mu_obs = 1, sd_obs = 1,
                    sd_OLRE = 4,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 3,
                    mu_b_pocc_2019 = 0.1, sd_b_pocc_2019 = 2,
                    u_sd_block = 3)
               )

model <- "scripts/JAGS/bird_JAGS_bpo_bin.R"

start <- Sys.time()

## Parallel computing:
cl <- makePSOCKcluster(3) ## On 3 cores

jm <- parJagsModel(cl = cl, 
                   name = "bpo_bin",
                   file = model,
                   data = data,
                   n.adapt = 1000, 
                   inits = inits,
                   n.chains = 3) 

parUpdate(cl = cl, object = "bpo_bin", n.iter = 1000)

samples <- 1000
n.thin <- 2

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
  write(., paste0("results/params_", TT, "_", ref, "_hp.txt"))
capture.output(summary(zc2), HPDinterval(zc2, prob = 0.95)) %>%
  write(., paste0("results/params_", TT, "_", ref, "_p.txt"))

## 4. Validate the model and export validation data and figures ----------------

pdf(paste0("figures/", TT, "_", ref, "_hp.pdf"))
plot(zc1); gelman.plot(zc1) 
dev.off()
pdf(paste0("figures/", TT, "_", ref, "_p.pdf")) 
plot(zc2); gelman.plot(zc2) 
dev.off()

capture.output(raftery.diag(zc1), 
               heidel.diag(zc1), 
               gelman.diag(zc1),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., paste0("results/dign_", TT, "_", ref, "_hp.txt"))
capture.output(raftery.diag(zc2),
               heidel.diag(zc2),
               gelman.diag(zc2)) %>%
  write(., paste0("results/dign_", TT, "_", ref, "_p.txt"))

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

## Store true occurrences for modelling community parameters:

zc4 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = "occ_true_out",
                      n.iter = samples,
                      thin = n.thin)

save(zc4, file = "clean/occ_true_out.rda")

## For community and species level metrics:

zc5 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("CI_div_bd", "CI_ctr_bd", "BACI_bd",
                                         "CI_div_r", "CI_ctr_r", "BACI_r",
                                         "CI_div_sl", "CI_ctr_sl", "BACI_sl", 
                                         "CI_div_cm", "CI_ctr_cm", "BACI_cm"),
                      n.iter = samples,
                      thin = n.thin)

## Combine MCMC chains:
zc5 <- combine.mcmc(zc5)

## Extract slopes and add ecdf and species names:
BACI_sl <- as.data.frame(summary(zc5)$quantiles[, c("2.5%","50%","97.5%")])
BACI_sl$ecdf <- as.vector(apply(zc5, 2, function(x) 1-ecdf(x)(0)))
BACI_sl$identity <- c(rep("beta", length(data$eval)),
                      rep("cm", length(data$eval)),
                      rep("alpha", length(data$eval)),
                      sort(rep(levels(bpo$species), length(data$eval))))

## Add treatment*BACI indicator categorisation:
BACI_sl$treatment <- levels(forest[, TT])[data$eval]
BACI_sl$indicator <- c(rep("BACI", length(data$eval)*(data$nspecies+3)),
                       rep("CI_ctr", length(data$eval)*(data$nspecies+3)),
                       rep("CI_div", length(data$eval)*(data$nspecies+3)))

## Export the data set for figures:
write.csv(BACI_sl, paste0("clean/BACI_sl_", TT, "_", ref, ".csv"))

end <- Sys.time()
end - start

## -------------------------------END-------------------------------------------
