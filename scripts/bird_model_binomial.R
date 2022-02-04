## model birds in a hierarchical model
## TC is No forestry
## C is Complete retention
## URT is Understory retention thinning
## T is conventional thinning
## 
## First edit: 20191201
## Last edit: 20220204
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(magrittr)
library(reshape)
library(reshape2)
library(parallel)
library(dclone)
library(data.table)

## 2. Load data and define functions -------------------------------------------

occ <- read.csv("data/occ_2016to2019.csv")
occ <- occ[!occ$plot %in% c("plot_30", "plot_118", "plot_120", "plot_121"), ]

bpo <- read.csv("clean/bpo_double.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")

bird_data <- read.csv("data/bird_data.csv")
bird_data <- droplevels(bird_data[bird_data$short %in% bpo$species, ])

## Function to summarise BACI results:
posterior_summary <- function(x){
  c(quantile(x, c(0.025, 0.5, 0.975)), "ecdf" = 1-ecdf(x)(0))
}

## 3. Analyse if some birds always show up after 15 minutes --------------------

occ <- as.data.table(droplevels(occ[occ$obs_year == 2017, ]))

ttest_calc <- function(x){
  if(nrow(x) > 2){
    t.test(x$minutes_to_obs-15, alternative = "greater")$p.value
  } else{as.numeric(NA)}
}

capture.output(occ[occ$species %in% levels(bpo$species),
                   list("mean" = mean(minutes_to_obs),
                        "se" = sd(minutes_to_obs)/sqrt(nrow(.SD)),
                        "ttest_p_value" = ttest_calc(.SD),
                        "N" = nrow(.SD)),
                   by = "species"]) %>%
  write(., "results/time_of_first_occurrence.txt")
  
## 4. Prepare data for experiment analysis -------------------------------------

## Make numeric levels for the treatment, experiment, and block matrices:
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)
forest$treatment_num <- as.numeric(forest$treatment)
forest$block_num <- as.numeric(forest$block)

## Create model data set:
data <- list(nobs = nrow(bpo),
             observed = bpo$n_obs,
             nvisits = bpo$n_visits,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             block = acast(forest, 
                           year ~ plot, 
                           value.var = "block_num")[1, ],
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             visibility = acast(forest, year ~ plot, value.var = "laser_median"),
             exp = acast(forest, year ~ plot, value.var = "exp_num"),
             treat = acast(forest, 
                           year ~ plot, 
                           value.var = "treatment_num")[1, ],
             C_plots = unique(forest$plot_num[forest$treatment == "C"]),
             T_plots = unique(forest$plot_num[forest$treatment == "T"]),
             TC_plots = unique(forest$plot_num[forest$treatment == "TC"]),
             URT_plots = unique(forest$plot_num[forest$treatment == "URT"])
             )

str(data)

## 5. Define and run the model -------------------------------------------------

## For true_occ:
T2 <- array(1, dim = c(max(data$species), max(data$year), max(data$site)))

## Initial values for random slopes:
inits <-  list(list(occ_true = T2,
                    mu_a_pdet = 0.5, sd_a_pdet = 3,
                    mu_b_vis = 1, sd_b_vis = 1,
                    mu_b_pdet_2018 = 0.5, sd_b_pdet_2018 = 5,
                    mu_b_pdet_2019 = 0.5, sd_b_pdet_2019 = 5,
                    mu_obs = 0.5, sd_obs = 5,
                    sd_OLRE = 2,
                    mu_a_pocc = matrix(0.5, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(2.1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.5, sd_b_pocc_2018 = 5,
                    mu_b_pocc_2019 = 0.5, sd_b_pocc_2019 = 5,
                    param_sd_block = 0.5),
               list(occ_true = T2,
                    mu_a_pdet = -0.5, sd_a_pdet = 2,
                    mu_b_vis = -1, sd_b_vis = 1,
                    mu_b_pdet_2018 = 0.7, sd_b_pdet_2018 = 3,
                    mu_b_pdet_2019 = 0.7, sd_b_pdet_2019 = 3,
                    mu_obs = -0.5, sd_obs = 1.5,
                    sd_OLRE = 3,
                    mu_a_pocc = matrix(0, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 1,
                    mu_b_pocc_2019 = 0.8, sd_b_pocc_2019 = 4,
                    param_sd_block = 4.5),
               list(occ_true = T2,
                    mu_a_pdet = 0, sd_a_pdet = 1,
                    mu_b_vis = 0.1, sd_b_vis = 1,
                    mu_b_pdet_2018 = 0.1, sd_b_pdet_2018 = 1,
                    mu_b_pdet_2019 = 0.1, sd_b_pdet_2019 = 1,
                    mu_obs = 1, sd_obs = 1,
                    sd_OLRE = 4,
                    mu_a_pocc = matrix(3, max(data$treat), max(data$exp)),
                    sd_a_pocc = matrix(0.1, max(data$treat), max(data$exp)),
                    mu_b_pocc_2018 = 0.1, sd_b_pocc_2018 = 3,
                    mu_b_pocc_2019 = 0.1, sd_b_pocc_2019 = 2,
                    param_sd_block = 3)
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

parUpdate(cl = cl, object = "bpo_bin", n.iter = 10000)

samples <- 250000
n.thin <- 250

zc1 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mu_a_pdet", "sd_a_pdet",
                                         "mu_b_vis", "sd_b_vis",
                                         "mu_b_pdet_2018", "sd_b_pdet_2018",
                                         "mu_b_pdet_2019", "sd_b_pdet_2019",
                                         "mu_obs", "sd_obs",
                                         "sd_OLRE",
                                         "mu_a_pocc", "sd_a_pocc",
                                         "mu_b_pocc_2018", "sd_b_pocc_2018",
                                         "mu_b_pocc_2019", "sd_b_pocc_2019",
                                         "sd_block"),
                      n.iter = samples, thin = n.thin)

zc2 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = "b_vis",
                      n.iter = samples, thin = n.thin)

## Export parameter estimates:
capture.output(bird_data$long, summary(zc1), HPDinterval(zc1, prob = 0.95)) %>% 
  write(., "results/hparams_vis.txt")
capture.output(bird_data$long, summary(zc2), HPDinterval(zc2, prob = 0.95)) %>%
  write(., "results/det_params_vis.txt")

## 6. Validate the model and export validation data and figures ----------------

pdf("figures/plot_hparams_vis2.pdf"); plot(zc1); gelman.plot(zc1); dev.off()

capture.output(raftery.diag(zc1), 
               heidel.diag(zc1), 
               gelman.diag(zc1),
               cor(data.frame(combine.mcmc(zc1)))) %>% 
  write(., "results/diagn_hparams_vis.txt")

## Produce validation metrics:
zc3 <- parCodaSamples(cl = cl, model = "bpo_bin",
                      variable.names = c("mean_obs", "mean_sim", "p_mean",
                                         "sd_obs_test", "sd_sim", "p_sd",
                                         "fit", "fit_sim", "p_fit"),
                      n.iter = samples,
                      thin = n.thin)

zj_val <- data.frame(combine.mcmc(zc3))

pdf("figures/plot_ppc_vis.pdf")
## Fit of mean:
plot(zj_val$mean_obs, zj_val$mean_sim,
     xlab = "mean real", ylab = "mean simulated",
     cex = .05)
abline(0, 1)
## Fit of variance:
plot(zj_val$sd_obs_test, zj_val$sd_sim,
     xlab = "sd real", ylab = "sd simulated",
     cex = .05)
abline(0,1)
## Overall fit:
plot(zj_val$fit, zj_val$fit_sim,
     xlab = "ssq real", ylab = "ssq simulated",
     cex = .05)
abline(0,1)
dev.off()

## 7. Extract informtion from posterior and create arrays ----------------------

## Extract the thinned MCMC chains for all species and intercepts:
zc_post2 <- parCodaSamples(cl = cl, model = "bpo_bin",
                           variable.names = "a_pocc",
                           n.iter = samples,
                           thin = n.thin)
save(zc_post2, file = "temp/zc_post2_vis.r")
load("temp/zc_post2_vis.r")
zcp2 <- combine.mcmc(zc_post2, collapse.chains = TRUE) 

## Transform values to the probability scale because a logistic GLMM was used:
zcp2_trans <- inv.logit(zcp2)

## To calculate all the results and indicators,create an array with the format 
## iteration*species(k in model)*treatment(m)*experiment(n): 
## The dimensions follow the column names ordering of the output matrix.
## The array is filled along the dimensions (e.g.: row, column, array, ...).
a_zcp2t <- array(as.vector(zcp2_trans), 
                dim = c(nrow(zcp2_trans), max(data$species), 4, 2))

## Adjust treatment and reference here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
levels(forest$treatment)
eval <- c(1, 2, 4); ref <- 3
# eval <- c(2, 4); ref <- 1
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Calculate differences in bird species occurrence between control and  -------
## treatments before the experiment:
if(ref == 3){
  treat_comp <- apply(a_zcp2t[,,,1], c(1,2), function(x) x[ref] - mean(x[eval]))
  treat_comp <- apply(treat_comp, 2, posterior_summary)
  treat_comp <- cbind("species" = bird_data$long, as.data.frame(t(treat_comp)))
  write.csv(treat_comp, "results/treatment_comparison_vis.csv", row.names = FALSE)
}

## Calculate the BACI indicators for all species and for all iterations: -------
BACI_out2 <- array(NA, dim = c(dim(a_zcp2t)[1:2], length(eval), 3))
for(m in 1:length(eval)){
  BACI_out2[,,m,1] <- (a_zcp2t[,,eval[m],2] - a_zcp2t[,,eval[m],1]) -    
                      (a_zcp2t[,,ref,2] - a_zcp2t[,,ref,1])        ## BACI
  BACI_out2[,,m,2] <- abs(a_zcp2t[,,eval[m],2] - a_zcp2t[,,eval[m],1]) - 
                      abs(a_zcp2t[,,ref,2] - a_zcp2t[,,ref,1])     ## CI-ctrl
  BACI_out2[,,m,3] <- abs(a_zcp2t[,,eval[m],2] - a_zcp2t[,,ref,2]) -     
                      abs(a_zcp2t[,,eval[m],1] - a_zcp2t[,,ref,1]) ## CI-div
}

## Keep track of the names:
dimnames(BACI_out2) <- list(NULL, 
                            bird_data$long, 
                            levels(forest$treatment)[eval], 
                            c("BACI", "CI_ctr", "CI_div"))

## Calculate community mean and species level results:
BACI_sl <- BACI_out2
cm <- apply(BACI_out2, c(1, 3, 4), mean)
dim(cm) <- c(nrow(a_zcp2t), 1, length(eval), 3)
dimnames(cm)[2] <- list("Community mean")
BACI_sl <- abind::abind(BACI_sl, cm, along = 2)
BACI_sl <- apply(BACI_sl, c(2, 3, 4), posterior_summary)
BACI_sl <- dcast(melt(BACI_sl), Var2 + Var3 + Var4 ~ Var1)

## Add naming for figures later on:
BACI_sl$ref <- paste0("ref_", levels(forest$treatment)[ref])
colnames(BACI_sl)[1:3] <- c("species", "treatment", "indicator")

## Export BACI_sl and adjust name according to the chosen reference:
write.csv(BACI_sl, 
          paste0("clean/BACI_sl_vis_ref_", 
                 ifelse(ref == 3, "NF", "CR"), 
                 ".csv"),
          row.names = FALSE)

## Calculate guild results: ----------------------------------------------------
bird_data$numeric <- as.numeric(bird_data$short) ## Numeric according to "short"
BACI_gl <- melt(bird_data[, c(2, 4:9)], 
                id.vars = c("long", "numeric"),
                variable.name = "group",
                value.name = "guild")
BACI_gl <- as.data.table(BACI_gl)
gl_calc <- function(x){
  T1 <- apply(BACI_out2[,x$numeric,,], c(1, 3, 4), mean)
  T1 <- apply(T1, c(2, 3), posterior_summary)
  T1 <- dcast(melt(T1), Var2 + Var3 ~ Var1)
  return(T1)
}
BACI_gl <- BACI_gl[, gl_calc(.SD), by = c("group", "guild")]

## Add naming for figures later on:
BACI_gl$ref <- paste0("ref_", levels(forest$treatment)[ref])
colnames(BACI_gl)[3:4] <- c("treatment", "indicator")

## Export BACI_gl and adjust name according to the chosen reference
write.csv(BACI_gl, 
          paste0("clean/BACI_gl_vis_ref_", 
                 ifelse(ref == 3, "NF", "CR"), 
                 ".csv"),
          row.names = FALSE)

end <- Sys.time()
end - start

## -------------------------------END-------------------------------------------

