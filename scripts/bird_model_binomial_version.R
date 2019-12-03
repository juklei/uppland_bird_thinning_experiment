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

# ## Create table with nr.obswrved/non-observed per species:
# T1 <- table(bpo[, c("species", "n_obs")])
# 
# ## Exclude species which where seen during at least 20% of all visits:
# red_names <- names(which(T1[, 1] < 100))
# bpo <- droplevels(bpo[bpo$species %in% red_names, ])

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
             treat = treat, 
             exp = exp) 

str(data)

## Define initial values:

## For true_occ:
T2 <- array(1, dim = c(data$nspecies, data$nyears, data$nsites))

## Initial values for random slopes:
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, 
                    u_sd_a_pdet = 0.5,
                    b_observer = 0,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(5, max(data$treat), max(data$exp)),
                    u_sd_year = 5,
                    u_sd_site = 5),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, 
                    u_sd_a_pdet = 4,
                    b_observer = -0.5,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    u_sd_year = 1,
                    u_sd_site = 1),
               list(occ_true = T2,
                    mu_a_pdet = 0, 
                    u_sd_a_pdet = 1,
                    b_observer = 0.8,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
                    u_sd_year = 0.5,
                    u_sd_site = 0.5)
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

parUpdate(cl = cl, object = "bpo_bin", n.iter = 15000)

samples <- 10000
n.thin <- 20

zc1 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_a_pdet", 
                                         "sd_a_pdet", "u_sd_pdet_year",
                                         "mu_a_pocc", "sd_a_pocc", 
                                         "u_sd_year", "u_sd_site"),
                      n.iter = samples, thin = n.thin)

zc2 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("a_pdet", "sd_pdet_year", "b_observer",
                                         "a_pocc", "sd_year", "sd_site"),
                      n.iter = samples, thin = n.thin)

end <- Sys.time()
end - start

## Export parameter estimates:
capture.output(summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_hyperparams.txt")
capture.output(summary(zc2), HPDinterval(zc2, prob = 0.95)) %>% 
  write(., "results/parameters_binomial_params.txt")

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/plot_binomial_hp.pdf"); plot(zc1); dev.off()
pdf("figures/plot_binomial_p.pdf"); plot(zc2); dev.off()

capture.output(raftery.diag(zc1), heidel.diag(zc1), gelman.diag(zc1)) %>% 
  write(., "results/diagnostics_binomial_hyperparams.txt")
capture.output(raftery.diag(zc2), heidel.diag(zc2), gelman.diag(zc2)) %>% 
  write(., "results/diagnostics_binomial_params.txt")

## Produce validation metrics:
zc3 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mean_obs", "mean_sim", "p_mean",
                                         "cv_obs", "cv_sim", "p_cv",
                                         "fit", "fit_sim", "p_fit"),
                      n.iter = samples,
                      thin = n.thin)

zj_val <- as.data.frame(rbind(zc3[[1]], zc3[[2]], zc3[[3]]))

## Fit of mean:
plot(zj_val$mean_obs, zj_val$mean_sim,
     xlab = "mean real", ylab = "mean simulated",
     cex = .05)
abline(0, 1)
text(x = 2.5, y = 3, paste0("P=", round(mean(zj_val$p_mean), 3)), cex = 1.5)

## Fit of variance:
plot(zj_val$cv_obs, zj_val$cv_sim, 
     xlab = "cv real", ylab = "cv simulated",
     cex = .05)
abline(0,1)
text(x = 0.7, y = 0.76, paste0("P=", round(mean(zj_val$p_cv), 3)), cex = 1.5)

## Overall fit:
plot(zj_val$fit, zj_val$fit_sim,
     xlab = "ssq real", ylab = "ssq simulated",
     cex = .05)
abline(0,1)
text(x = 850, y = 700, paste0("P=", round(mean(zj_val$p_fit)), 3), cex = 1.5)

## 6. Export -------------------------------------------------------------------

zc4 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("CI_div_C", "CI_ctr_C", "BACI_C", 
                                         "CI_div_T", "CI_ctr_T", "BACI_T",
                                         "CI_div_URT", "CI_ctr_URT", "BACI_URT",
                                         "CI_div_C_cm", "CI_ctr_C_cm", "BACI_C_cm", 
                                         "CI_div_T_cm", "CI_ctr_T_cm", "BACI_T_cm",
                                         "CI_div_URT_cm", "CI_ctr_URT_cm", "BACI_URT_cm",
                                         "CI_div_C_r", "CI_ctr_C_r", "BACI_C_r", 
                                         "CI_div_T_r", "CI_ctr_T_r", "BACI_T_r",
                                         "CI_div_URT_r", "CI_ctr_URT_r", "BACI_URT_r"),
                      n.iter = samples,
                      thin = n.thin)

##############################################################3
## Extract slopes from the linear model:
export_beta_tb <- summary(zc)$quantiles[, c("2.5%","50%","97.5%")][row_select,]

## Extract probability that slope above or below 0 for beta_td:
prob <- data.frame("cat" = colnames(zc[[1]])[row_select], 
                   "prob" = rep(NA, dim(obs)[2]))
for(i in 1:dim(obs)[2]){
  prob$prob[i] <- 1-ecdf(unlist(zc[, paste(prob$cat[i])]))(0)
}

## Add probabilities to export_beta_tb:
export_beta_tb <- cbind(export_beta_tb, prob)

## Add species names:
export_beta_tb$species <- unlist(dimnames(obs)[2])

## Export the data set for figures:
write.csv(export_beta_tb[, -4], "clean/prob&slope_beta_tb.csv", row.names = F)

## -------------------------------END-------------------------------------------
