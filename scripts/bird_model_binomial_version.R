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
library(coda)
library(magrittr)
library(reshape2)

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

# ## Create table with nr.obswrved/non-observed per species:
# T1 <- table(bpo[, c("species", "n_obs")])
# 
# ## Exclude species which where seen during at least 20% of all visits:
# red_names <- names(which(T1[, 1] < 150))
# bpo <- droplevels(bpo[bpo$species %in% red_names, ])

# bpo <- droplevels(bpo[bpo$species %in% c("tofss", "trapa"), ])
# bpo <- droplevels(bpo[bpo$observer == "jkn", ])

## Join "T" and "URT" to "T":
# levels(forest$treatment)[4] <- "T"

## Make numeric levels for the treatment*experiment matrix:
forest$treat_num <- as.numeric(forest$treatment)
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)

## Create data arrays for process model part:
treat <- acast(unique(forest[, c("plot", "treat_num")]), . ~ plot)[1, ]
exp <- acast(forest[, c("plot", "year", "exp_num")], year ~ plot)

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
             treat = treat, exp = exp) 

str(data)

## Define initial values:
T2 <- array(1, dim = c(data$nspecies, data$nyears, data$nsites))
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, sd_a_pdet = 0.5,
                    b_observer = 0,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(5, max(data$treat), max(data$exp)),
                    u_sd_year = 5, 
                    u_sd_site = 5),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, sd_a_pdet = 4,
                    b_observer = -0.5,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    u_sd_year = 1, 
                    u_sd_site = 1),
               list(occ_true = T2,
                    mu_a_pdet = 0, sd_a_pdet = 1,
                    b_observer = 0.8,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
                    u_sd_year = 0.5, 
                    u_sd_site = 0.5))

model <- "scripts/JAGS/bird_JAGS_bpo_bin.R"

start <- Sys.time()

jm <- jags.model(model,
                 data = data,
                 n.adapt = 10000, 
                 inits = inits,
                 n.chains = 3) 

burn.in <-  40000

update(jm, n.iter = burn.in) 

samples <- 50000
n.thin <- 100

zc1 <- coda.samples(jm,
                    variable.names = c("mu_a_pdet", "sd_a_pdet", "u_sd_pdet_year",
                                       "mu_a_pocc", "sd_a_pocc", 
                                       "u_sd_year", "u_sd_site"),
                    n.iter = samples, 
                    thin = n.thin)

zc2 <- coda.samples(jm,
                    variable.names = c("a_pdet", "sd_pdet_year", "b_observer",
                                       "a_pocc", "sd_year", "sd_site"),
                    n.iter = samples, 
                    thin = n.thin)

end <- Sys.time()
end - start

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_hyperparams.txt")
capture.output(summary(zc1), HPDinterval(zc2, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_params.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_binomial_hp.pdf"); plot(zc1); dev.off()
pdf("figures/plot_binomial_p.pdf"); plot(zc2); dev.off()

capture.output(raftery.diag(zc1), heidel.diag(zc1)) %>% 
  write(., "results/diagnostics_binomial_hyperparams.txt")
capture.output(raftery.diag(zc2), heidel.diag(zc2)) %>% 
  write(., "results/diagnostics_binomial_params.txt")

# ## Produce validation metrics: 
# zj_val <- jags.samples(jm, 
#                        variable.names = c("mean_nseen", 
#                                           "mean_nseen_sim",
#                                           "p_mean", 
#                                           "cv_nseen", 
#                                           "cv_nseen_sim", 
#                                           "p_cv", 
#                                           "fit", 
#                                           "fit_sim",
#                                           "p_fit"), 
#                        n.iter = samples, 
#                        thin = n.thin)
# 
# ## Fit of mean:
# plot(zj_val$mean_nseen, 
#      zj_val$mean_nseen_sim, 
#      xlab = "mean real", 
#      ylab = "mean simulated", 
#      cex = .05)
# abline(0, 1)
# p <- summary(zj_val$p_mean, mean)
# text(x = 0.6, y = 0.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Fit of variance:
# plot(zj_val$cv_nseen, 
#      zj_val$cv_nseen_sim, 
#      xlab = "cv real", 
#      ylab = "cv simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_cv, mean)
# text(x = 1.4, y = 1.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)
# 
# ## Overall fit:
# plot(zj_val$fit, 
#      zj_val$fit_sim, 
#      xlab = "ssq real", 
#      ylab = "ssq simulated", 
#      cex = .05)
# abline(0,1)
# p <- summary(zj_val$p_fit, mean)
# text(x = 850, y = 700, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Export -------------------------------------------------------------------

zj <- jags.samples(jm,
                   variable.names = c("CI_div_C", "CI_ctr_C", "BACI_C", 
                                      "CI_div_T", "CI_ctr_T", "BACI_T",
                                      "CI_div_URT", "CI_ctr_URT", "BACI_URT",
                                      "CI_div_C_cm", "CI_ctr_C_cm", "BACI_C_cm", 
                                      "CI_div_T_cm", "CI_ctr_T_cm", "BACI_T_cm",
                                      "CI_div_URT_cm", "CI_ctr_URT_cm", "BACI_URT_cm"),
                   n.iter = samples,
                   thin = n.thin)

## -------------------------------END-------------------------------------------
