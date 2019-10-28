## model bird richness in a hierarchical model
## Calculate this script for each year seperately
## 
## First edit: 20190628
## Last edit: 20191023
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

dir.create("results")
dir.create("figures")

## Print all rows for mcmc outputs
options(max.print = 10E5)

## Backscale function
backscale <- function(pred_data, model_input_data) {
  
  pred_data*attr(model_input_data, 'scaled:scale') + 
    attr(model_input_data, 'scaled:center')
  
}

## 3. Load and explore data ----------------------------------------------------

dir("clean")

bpo <- read.csv("clean/bpo_double.csv")
head(bpo)
str(bpo)

## 4. The model ----------------------------------------------------------------

## Do everything for one year and observer each, except 2016, where no ds data 
## is available:
bpo <- droplevels(bpo[bpo$obs_year == 2018, ])

## Create data arrays:

nvisits <- acast(bpo[, c("plot", "species", "n_visits")],
                 formula = plot ~ species, 
                 value.var = "n_visits")

nseen <- acast(bpo[, c("plot", "species", "n_obs")],
               formula = plot ~ species, 
               value.var = "n_obs")

## Create model data set:
data <- list(nsites = nlevels(bpo$plot),
             nspecies = nlevels(bpo$species),
             nvisits = nvisits,
             nseen = nseen,
             R = matrix(c(5,0,0,1), ncol = 2),
             df = 3)

str(data)

## Inits for occ_true
T1 <- data$nseen
T1 <- ifelse(T1 > 0, 1, 0)

inits <-  list(list(occ_true = T1,
                    Omega = matrix(c(1,0,0,1), ncol = 2),
                    eta = matrix(0, nrow = data$nspecies, ncol = 2))
               )

model <- "scripts/JAGS/bird_JAGS_detection_part.R"

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits, 
                 n.chains = 1) 

burn.in <-  10000

update(jm, n.iter = burn.in) 

samples <- 10000
n.thin <- 5

zc <- coda.samples(jm,
                   variable.names = c("mu_eta",
                                      "probs",
                                      "p_occ",
                                      "p_det",
                                      "Sigma"), 
                   n.iter = samples, 
                   thin = n.thin)

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_bird_det_2018.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_bird_det_2018.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_bird_det_2018.txt")

## Produce validation metrics: 
zj_val <- jags.samples(jm, 
                       variable.names = c("mean_nseen", 
                                          "mean_nseen_sim",
                                          "p_mean", 
                                          "cv_nseen", 
                                          "cv_nseen_sim", 
                                          "p_cv", 
                                          "fit", 
                                          "fit_sim",
                                          "p_fit"), 
                       n.iter = samples, 
                       thin = n.thin)

## Fit of mean:
plot(zj_val$mean_nseen, 
     zj_val$mean_nseen_sim, 
     xlab = "mean real", 
     ylab = "mean simulated", 
     cex = .05)
abline(0, 1)
p <- summary(zj_val$p_mean, mean)
text(x = 0.6, y = 0.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_nseen, 
     zj_val$cv_nseen_sim, 
     xlab = "cv real", 
     ylab = "cv simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_cv, mean)
text(x = 1.4, y = 1.7, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, 
     zj_val$fit_sim, 
     xlab = "ssq real", 
     ylab = "ssq simulated", 
     cex = .05)
abline(0,1)
p <- summary(zj_val$p_fit, mean)
text(x = 850, y = 700, paste0("P=", round(as.numeric(p[1]), 4)), cex = 1.5)

## 6. Produce output and export data for seperate analysis ---------------------

## Produce predictions:
zj_pred <- jags.samples(jm, 
                        variable.names = c("richness", "scaled_rb"),
                        n.iter = samples, 
                        thin = n.thin)

zj_pred$plotnames <- levels(bpo$plot)

save(zj_pred, file = "clean/rb_2018.rdata")

## -------------------------------END-------------------------------------------
