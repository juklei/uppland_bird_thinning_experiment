## model birds in a hierarchical model
## 
## First edit: 20191031
## Last edit: 20191106
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
bpo <- read.csv("clean/bpo_double_bernoulli.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")
head(bpo)
str(bpo)
head(forest)

## 4. The model ----------------------------------------------------------------

## Create table with nr.obswrved/non-observed per species:
T1 <- table(bpo[, c("species", "observed")])

## Exclude species which where seen > 100 and not seen at > 100 times:
red_names <- names(which(T1[, "1"] > 50 & T1[, "0"] > 50))
bpo <- droplevels(bpo[bpo$species %in% red_names, ])

bpo <- droplevels(bpo[bpo$species %in% c("tofss", "trapa"), ])

## Join "T" and "URT" to "T":
levels(forest$treatment)[4] <- "T"

## Add ldm:
bpo$ldm <- ifelse(bpo$species %in% ldm$species, 2, 1) 

## Scale all continuous values:
bpo$dpm_scaled <- scale(bpo$dp_march)
bpo$mps_scaled <- scale(bpo$min_post_sunrise)
forest$sdbh_scaled <- scale(forest$average_dbh_all_alive)
forest$nr_lov_scaled <- scale(forest$nr_lov)
forest$nr_skarm_scaled <- scale(forest$nr_skarm)
forest$nr_sd_scaled <- scale(forest$nr_staende_dodved)
forest$lm_scaled <- scale(forest$laser_mean)

## Make numeric levels for the treatment*experiment matrix:
forest$treat_num <- as.numeric(forest$treatment)
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)

## Create data arrays:

treat <- acast(unique(forest[, c("plot", "treat_num")]), . ~ plot)[1, ]

exp <- acast(forest[, c("plot", "year", "exp_num")], year ~ plot)

sdbh <- acast(forest[, c("plot", "year", "sdbh_scaled")], year ~ plot)

dec <- acast(forest[, c("plot", "year", "nr_lov_scaled")], year ~ plot)

umbr <- acast(forest[, c("plot", "year", "nr_skarm_scaled")], year ~ plot)

lm <- acast(forest[, c("plot", "year", "lm_scaled")], year ~ plot)

## Create model data set:
data <- list(nobs = nrow(bpo),
             nspecies = max(as.numeric(bpo$species)),
             nyears = dim(sdbh)[1],
             nsites = dim(sdbh)[2],
             observed = bpo$observed,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             dpm = bpo$dpm_scaled,
             mps = bpo$mps_scaled,
             ldm = bpo$ldm,
             treat = treat,
             exp = exp,
             sdbh = sdbh,
             dec = dec,
             umbr = umbr,
             lm = lm) 

str(data)

inits <-  list(list(occ_true = array(1, dim = c(data$nspecies,
                                                data$nyears, 
                                                data$nsites)),
                    mu_a_pdet = 0.5, sd_a_pdet = 0.5,
                    mu_b_pdet_2018 = 0.5, sd_b_pdet_2018 = 5,
                    mu_b_pdet_2019 = 0.5, sd_b_pdet_2019 = 5,
                    b_observer = 0,
                    b_dpm = c(0, 0),
                    b2_dpm = c(0, 0),
                    b_mps = 0,
                    b2_mps = 0,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(5, max(data$treat), max(data$exp)),
                    u_sd_year = 5, #rep(5, data$nspecies),
                    u_sd_site = 5, #rep(5, data$nspecies),
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5,
                    mu_b_sdbh = 0.5, sd_b_sdbh = 5, 
                    mu_b_dec = 0.5, sd_b_dec = 5,
                    mu_b_umbr = 0.5, sd_b_umbr = 5,
                    mu_b_lm = 0.5, sd_b_lm = 5))

model <- "scripts/JAGS/bird_JAGS_bpo.R"

start <- Sys.time()

jm <- jags.model(model,
                 data = data,
                 n.adapt = 1000, 
                 inits = inits,
                 n.chains = 1) 

burn.in <-  1000

update(jm, n.iter = burn.in) 

samples <- 1000
n.thin <- 2

zc1 <- coda.samples(jm,
                    variable.names = c("mu_a_pdet", "sd_a_pdet", 
                                       "mu_b_pdet_2018", "sd_b_pdet_2018", 
                                       "mu_b_pdet_2019", "sd_b_pdet_2019",
                                       "b_observer",
                                       "b_dpm", "b2_dpm", "b_mps", "b2_mps", 
                                       "mu_a_pocc", "sd_a_pocc", 
                                       "u_sd_year", "u_sd_site",
                                       "mu_b_pocc_2018", "sd_b_pocc_2018", 
                                       "mu_b_pocc_2019", "sd_b_pocc_2019",
                                       "mu_b_sdbh", "sd_b_sdbh",
                                       "mu_b_dec", "sd_b_dec", 
                                       "mu_b_umbr", "sd_b_umbr", 
                                       "mu_b_lm", "sd_b_lm"),
                    n.iter = samples, 
                    thin = n.thin)

end <- Sys.time()
end - start

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/parameters_1part.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_1part.pdf")
plot(zc1)
dev.off()

capture.output(raftery.diag(zc1), heidel.diag(zc1)) %>% 
  write(., "results/diagnostics_1part.txt")

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

## 6. Export for graphing ------------------------------------------------------

##...

## -------------------------------END-------------------------------------------
