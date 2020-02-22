## This script prepares the data from the sticky traps and runs the analysis
## in the BACI experiment
##
## First edit: 20200218
## Last edit: 20200219
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require(data.table)
require(rjags)
require(runjags)
require(coda)
require(dplyr)

## 2. Load and explore data ----------------------------------------------------

dir("data")
i_2017 <- read.csv("data/insects_2017.csv")
i_2018 <- read.csv("data/insects_2018.csv")
i_2019 <- read.csv("data/insects_2019.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")

## This function standardises values to between 0 and 1
standardise <- function(x) {
  (x - min(x[], na.rm = T)) / (max(x[], na.rm = T) - min(x[], na.rm = T))
}


## 3. Prepare insect data and extract acc. cover and date of highest -----------
##    increase and merge with forest data increase

i_comb <- as.data.table(rbind(i_2017, i_2018, i_2019))
i_comb$mi_stand[i_comb$obs_year == 2019] <- standardise(i_comb$mean_incr[i_comb$obs_year == 2019])
i_comb$mi_stand[i_comb$obs_year != 2019] <- standardise(i_comb$mean_incr[i_comb$obs_year != 2019])

## Include only days which are covered across all plots (exclude 120):
temp <- i_comb[i_comb$plot != "plot_120",
               list("min" = min(post_march), "max" = max(post_march)), 
               by = c("plot", "obs_year")]
pm_limits <- c(max(temp$min), min(temp$max))
rm(temp)
i_comb <- i_comb[i_comb$post_march >= pm_limits[1] & 
                 i_comb$post_march <= pm_limits[2], ]

## Calculte sum of increment per plot and year:
i_out <- i_comb[, list("acc_mi" = sum(mean_incr),
                       "acc_mi_st" = sum(mi_stand),
                       "pm_max" = mean(post_march[which(mean_incr == 
                                                        max(mean_incr))])),
                by = c("plot", "obs_year")]

## Merge with forest data required for the BACI analysis:

if_comb <- merge(i_out, forest[, c("plot", "year", "treatment", "experiment")], 
                 by.x = c("plot", "obs_year"), 
                 by.y = c("plot", "year"))

## 4. Prepare the model data, the inits and load the model ---------------------

## Create model data set:
data <- list(nobs = nrow(if_comb),
             nsites = nlevels(if_comb$plot),
             cover = if_comb$acc_mi_st,
             exp = ifelse(if_comb$experiment == "before", 1, 2),
             treat = as.numeric(if_comb$treatment),
             site = as.numeric(if_comb$plot),
             year_2018 = ifelse(if_comb$obs_year == 2018, 1, 0),
             year_2019 = ifelse(if_comb$obs_year == 2019, 1, 0),
             eval = c(1, 2, 4),
             ref = 3) 

## Create inits:

inits <- list(list(a_cov = matrix(1, max(data$treat), max(data$exp)),
                   sigma = 5,
                   b_2018 = 0,
                   b_utb_2019 = 1,
                   sigma_site = 1),
              list(a_cov = matrix(100, max(data$treat), max(data$exp)),
                   sigma = 1,
                   b_2018 = -5,
                   b_utb_2019 = 0.1,
                   sigma_site = 1),
              list(a_cov = matrix(10, max(data$treat), max(data$exp)),
                   sigma = 1,
                   b_2018 = 5,
                   b_utb_2019 = 5,
                   sigma_site = 10)
              )

model <- "scripts/JAGS/insect_JAGS.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 3) 

burn.in <-  50000
update(jm, n.iter = burn.in) 

samples <- 20000
n.thin <- 10

zc <- coda.samples(jm,
                   variable.names = c("a_cov", "sigma", "sigma_site",
                                      "b_2018", "b_utb_2019"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_insect_cover.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_insect_cover.pdf")
plot(zc); gelman.plot(zc) 
dev.off()

capture.output(raftery.diag(zc), 
               heidel.diag(zc), 
               gelman.diag(zc),
               cor(data.frame(combine.mcmc(zc)))) %>% 
  write(., "results/diagnostics_insect_cover.txt")

## Produce validation metrics: 
zj_val <- jags.samples(jm, 
                       variable.names = c("mean_obs", "mean_sim","p_mean", 
                                          "cv_obs", "cv_sim", "p_cv", 
                                          "fit", "fit_sim", "p_fit"), 
                       n.iter = samples, 
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean_obs, 
     zj_val$mean_sim, 
     xlab = "mean real", 
     ylab = "mean simulated", 
     cex = .05)
abline(0, 1)
p <- summary(zj_val$p_mean, mean)
text(x = 0.6, y = 0.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_obs, 
     zj_val$cv_sim, 
     xlab = "cv real", 
     ylab = "cv simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_cv, mean)
text(x = 1, y = .8, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, 
     zj_val$fit_sim, 
     xlab = "ssq real", 
     ylab = "ssq simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_fit, mean)
text(x = 0.2, y = 0.4, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Export BACI indicators ---------------------------------------------------

zj_out <- coda.samples(jm, 
                       variable.names = c("BACI", "CI_div", "CI_ctr"),
                       n.iter = samples, 
                       thin = n.thin)

summary(zj_out)

## -------------------------------END-------------------------------------------
