## model birds in a hierarchical model
## 
## First edit: 20191031
## Last edit: 20191031
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
forest <- read.csv("clean/forest_experiment_data.csv")
head(bpo)
str(bpo)

## 4. The model ----------------------------------------------------------------

## Reduce to two species for trials:
# bpo <- bpo[bpo$species %in% c("bofik", "trapa"), ]

## Scale all continuous values:
bpo$dpm_scaled <- scale(bpo$dp_march)
bpo$mps_scaled <- scale(bpo$min_post_sunrise)
forest$sdbh_scaled <- scale(forest$average_dbh_all_alive)

## Create data arrays:

seen <- acast(bpo[, c("plot", "obs_year", "species", "visit", "observed")],
              formula = species ~ obs_year ~ plot ~ visit,
              value.var = "observed")

dpm <- acast(bpo[, c("plot", "obs_year", "species", "visit", "dpm_scaled")],
             formula = species ~ obs_year ~ plot ~ visit,
             value.var = "dpm_scaled")
dpm[is.na(dpm)] <- 0 ## NA in covariates become dummy mean (0, beacuse scaled) 

mps <- acast(bpo[, c("plot", "obs_year", "species", "visit", "mps_scaled")],
             formula = species ~ obs_year ~ plot ~ visit,
             value.var = "mps_scaled")
mps[is.na(mps)] <- 0 ## NA in covariates become dummy mean (0, beacuse scaled) 

ldm <- ifelse(unlist(dimnames(seen)[1]) %in% ldm$species, 2, 1)
  
# sdbh <- acast(forest[forest$experiment %in% c("before", "after"), 
#                      c("plot", "experiment", "average_dbh_all_alive")],
#               formula = experiment ~ plot,
#               value.var = "average_dbh_all_alive")
# sdbh["after", is.na(sdbh["after", ])] <- sdbh["before", is.na(sdbh["after", ])]
# 
# experiment <-  
#   
  
## Create model data set:
data <- list(nspecies = dim(seen)[1],
             nyears = dim(seen)[2],
             nsites = dim(seen)[3],
             nvisits = dim(seen)[4],
             seen = seen,
             dpm = dpm,
             mps = mps,
             ldm = ldm#,
             # sdbh = t(scale(t(sdbh)))
             ) 

str(data)

inits <-  list(list(occ_true = array(1, dim = c(data$nspecies,
                                                data$nyears, 
                                                data$nsites)),
                    pocc = array(0.5, dim = c(data$nspecies, data$nyears)),
                    a_pdet = rep(0.1, data$nspecies),
                    b_dpm = c(0, 0),
                    b2_dpm = c(0, 0),
                    b_mps = c(0, 0),
                    b2_mps = c(0, 0),
                    mu_a_pdet = 0.5,
                    sd_a_pdet = 0.5,
                    mu_b_dpm = 0.5,
                    sd_b_dpm = 0.5,
                    mu_b2_dpm = 0.5,
                    sd_b2_dpm = 0.5,
                    mu_b_mps = 0.5,
                    sd_b_mps = 0.5,
                    mu_b2_mps = 0.5,
                    sd_b2_mps = 0.5))

model <- "scripts/JAGS/bird_JAGS_bpo_bernoulli.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits,
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 20000
n.thin <- 10

zc <- coda.samples(jm,
                   variable.names = c("pocc",
                                      "a_pdet",
                                      "b_dpm", 
                                      "b2_dpm",
                                      "b_mps",
                                      "b2_mps"),
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_bernoulli.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bernoulli.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_bernoulli.txt")

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

## 6. Produce output and export data for seperate analysis ---------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("richness", "scaled_rb"),
                        n.iter = samples, 
                        thin = n.thin)

zj_pred$plotnames <- levels(bpo$plot)

save(zj_pred, file = "clean/rb_2018.rdata")

## -------------------------------END-------------------------------------------
