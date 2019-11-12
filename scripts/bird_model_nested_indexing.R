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
library(reshape2)
library(data.table)
library(dplyr)

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

## 4. Construct the data set ---------------------------------------------------

# ## Create table with nr.obswrved/non-observed per species:
# T1 <- table(bpo[, c("species", "observed")])
# 
# ## Exclude species which where seen > 100 and not seen at > 100 times:
# red_names <- names(which(T1[, "1"] > 80 & T1[, "0"] > 80))
# bpo <- bpo[bpo$species %in% red_names, ]

## Reduce to two species for trials:
bpo <- droplevels(bpo[bpo$species %in% c("bofik", "trapa"), ])

## Merge bpo, ldm, and forest to one data frame:
bpof <- merge(bpo, 
              forest, 
              by.x = c("plot", "block", "obs_year"),
              by.y = c("plot", "block", "year"))
bpof$ldm <- ifelse(bpof$species %in% ldm$species, 2, 1) 

## Add binary identifier for treatments:
bpof$thinned <- ifelse(bpof$treatment %in% c("T", "URT") & 
                         bpof$experiment == "after", 1, 0)
bpof$control <- ifelse(bpof$treatment == "C" & bpof$experiment == "after", 1, 0)

## Add hierarchical level identifiers (HLI):

bpof <- as.data.table(bpof)

## HLI for plot*year*species:
bpof[ , "pys" := .GRP, by = c("plot", "obs_year", "species")]

## Create plot*year*species unique data set for H2:
dH2 <- unique(bpof[, c("pys", "species", "plot", "obs_year", "thinned", 
                       "control", "average_dbh_all_alive")
                  ])

## HLI for year*species:
dH2[ , "ys" := .GRP, by = c("obs_year", "species")]

## HLI for plot*species:
dH2[ , "ps" := .GRP, by = c("plot", "species")]

## Create a plot*species unique data set for H3:
dH3p <- unique(dH2[, c("ps", "species")])

## Create a year*species unique data set for H3:
dH3y <- unique(dH2[, c("ys", "species")])

## Create model data set:
data <- list(## Observational data:
             nobs = nrow(bpof), 
             pys = bpof$pys,
             sH1 = as.numeric(bpof$species),
             ldm = bpof$ldm,
             observed = bpof$observed,
             dpm = scale(bpof$dp_march)[, 1],
             mps = scale(bpof$min_post_sunrise)[, 1],
             ## Process model:
             npys = max(dH2$pys), 
             ys = dH2$ys,
             ps = dH2$ps,
             sH2 = as.numeric(dH2$species),
             thinned = dH2$thinned,
             control = dH2$control,
             sdbh = scale(dH2$average_dbh_all_alive)[, 1],
             ## Grouping effects:
             nps = max(dH2$ps),
             nys = max(dH2$ys),
             sH3p = as.numeric(dH3p$species),
             sH3y = as.numeric(dH3y$species),
             ## Priors:
             ns = max(as.numeric(bpof$species))) 

str(data)

inits <-  list(list(occ_true = rep(1, data$npys),
                    mu_a_pdet = 0.5,
                    sd_a_pdet = 0.5,
                    b_dpm = c(0, 0),
                    b2_dpm = c(0, 0),
                    b_mps = 0,
                    b2_mps = 0,
                    mu_a_pocc = 0.5,
                    sd_a_pocc = 5,
                    sd_year = rep(5, data$ns),
                    sd_site = rep(5, data$ns),
                    mu_b_thinned = 0.5,
                    sd_b_thinned = 5,
                    mu_b_control = 0.5,
                    sd_b_control = 5,
                    mu_b_sdbh = 0.5,
                    sd_b_sdbh = 5, 
                    mu_b_sdbh_t = 0.5,
                    sd_b_sdbh_t = 5,
                    mu_b_sdbh_c = 0.5,
                    sd_b_sdbh_c = 5))

model <- "scripts/JAGS/bird_JAGS_bpo_ni.R"

start <- Sys.time()

jm <- jags.model(model,
                 data = data,
                 n.adapt = 5000, 
                 inits = inits,
                 n.chains = 1) 

burn.in <-  5000

update(jm, n.iter = burn.in) 

samples <- 5000
n.thin <- 5

zc <- coda.samples(jm,
                   # variable.names = c("a_pdet", "b_dpm", "b2_dpm", "b_mps", 
                   #                    "b2_mps", "a_pocc", "sd_year", "sd_site",
                   #                    "b_thinned", "b_control", "b_sdbh", 
                   #                    "b_sdbh_t", "b_sdbh_c"),
                   variable.names = c("mu_a_pdet", "sd_a_pdet", "b_dpm", 
                                      "b2_dpm", "b_mps", "b2_mps", "mu_a_pocc",
                                      "sd_a_pocc", "sd_year", "sd_site",
                                      "mu_b_thinned", "sd_b_thinned", 
                                      "mu_b_control", "sd_b_control", 
                                      "mu_b_sdbh", "sd_b_sdbh", "mu_b_sdbh_t", 
                                      "sd_b_sdbh_t", "mu_b_sdbh_c", 
                                      "sd_b_sdbh_c"),
                   n.iter = samples, 
                   thin = n.thin)

end <- Sys.time()
end - start

## Export parameter estimates:
capture.output(summary(zc), HPDinterval(zc, prob = 0.95)) %>% 
  write(., "results/parameters_ni.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_zc_ni.pdf")
plot(zc)
dev.off()

capture.output(raftery.diag(zc), heidel.diag(zc)) %>% 
  write(., "results/diagnostics_ni.txt")

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
