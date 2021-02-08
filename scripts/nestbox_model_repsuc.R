## Analysis of reproductive success in the great tit
##
## First edit: 20210208
## Last edit: 20210208
##
## Author: Julian Klein

## 1. Clear environment, load libraries, and define funcitons ------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(reshape2)

## Function to summarise BACI results:
posterior_summary <- function(x){
  c(quantile(x, c(0.025, 0.5, 0.975)), "ecdf" = 1-ecdf(x)(0))
}

## 2. Load and explore data ----------------------------------------------------

dir("data")
d_nb <- read.csv("data/nestboxes_reproduction.csv")
d_pos <- read.csv("data/ffmf.csv")
head(d_nb)
head(d_pos)

## Remove data from blocks which were visited only at one occasion after the
## experiment, because many of these nestboxes have an unknown owner:
d_nb <- d_nb[!(d_nb$block %in% c("sodersjon", "stocksatra", "kyrkstigen", 
                                 "stenby", "soderlejde", "halsingdal", 
                                 "brudgatan") & d_nb$year == 2019), ]

## Also plot 33 &36 in 2019, because they were forgotten:
d_nb <- d_nb[!(d_nb$plot %in% c("plot_33", "plot_36") & d_nb$year == 2019), ]

## Reduce to information required for the model:
d_rs <- d_nb[, c("plot", "year", "box", "experiment", "treatment", "fledged")]
d_rs <- na.omit(droplevels(d_rs))

## Reduce d_pos to plot, box and coordinates and add to d_rs:
d_pos <- na.omit(d_pos[, c("plot", "box", "north", "east")])
d_rs <- merge(d_rs, d_pos, by = c("plot", "box"), all.x = TRUE)

## 3. Prepare data for the JAGS module -----------------------------------------

## Create a inverse distance matrix for calculating Moran's I for every year:
## Which year for Moran's I?
y_sel <- 2017
# y_sel <- 2018
# y_sel <- 2019
dm <- d_rs[d_rs$year == y_sel, c("north", "east")]
dm <- 1/as.matrix(dist(dm))
diag(dm) <- 0
seq <- which(d_rs$year == y_sel)

## Compile data for JAGS:
data <- list(nobs = nrow(d_rs),
             repsuc = d_rs$fledged,
             year = d_rs$year,
             treat = as.numeric(d_rs$treatment),
             exp = ifelse(d_rs$experiment == "before", 1, 2),
             dm = dm,
             seq = seq)

## 4. Prepare inits and run the model -------------------------------------------

## Prepare inits:
i_exp_effect <- matrix(1, max(data$treat), max(data$exp))

inits <- list(list(z = rep(1, data$nobs),
                   exp_effect_l = i_exp_effect,
                   exp_effect_p = i_exp_effect,
                   y_effect_l = c(1, 1),
                   y_effect_p = c(1, 1)),
              list(z = rep(1, data$nobs),
                   exp_effect_l = -i_exp_effect,
                   exp_effect_p = -i_exp_effect,
                   y_effect_l = c(-1, -1),
                   y_effect_p = c(-1, -1)),
              list(z = rep(1, data$nobs),
                   exp_effect_l = i_exp_effect*5,
                   exp_effect_p = i_exp_effect*5,
                   y_effect_l = c(5, 5),
                   y_effect_p = c(5, 5)))

## Load the model:
model <- "scripts/JAGS/nestbox_JAGS_repsuccess.R"

## Run and burn in:
jm <- jags.model(model, data, inits, 3, 5000)
update(jm, 5000)

## Sample from the parameter posteriors:
cs_1 <- coda.samples(jm, c("exp_effect_l", "exp_effect_p",
                           "y_effect_l", "y_effect_p"), 
                     50000, 50)

## 5. Validate the model and export validation data and figures ----------------

## Check mixing and convergence:
pdf("figures/nestbox_repsuc.pdf"); plot(cs_1); gelman.plot(cs_1); dev.off()

## Produce validation metrics for posterior predictive checks:
js_2 <- jags.samples(jm, c("mean_obs", "mean_sim",
                           "sd_obs", "sd_sim",
                           "fit", "fit_sim",
                           "I"), 5000, 5)

## Calculate Moran's I if random spatial distribution of residuals: 
E0 <- round(-1/(nrow(dm)*3 - 1), 4)

## Export validations as figure:
pdf("figures/nestbox_repsuc_ppc&MI.pdf")
par(mfrow = c(2, 2))
plot(js_2$mean_obs, js_2$mean_sim, xlab = "mean real", ylab = "mean simulated")
abline(0, 1)
plot(js_2$sd_obs, js_2$sd_sim, xlab = "sd real", ylab = "sd simulated")
abline(0,1)
plot(js_2$fit, js_2$fit_sim, xlab = "ssq real", ylab = "ssq simulated")
abline(0,1)
plot(density(js_2$I), main = paste0("Moran's I, ", y_sel, ", E0 = ", E0))
dev.off()

## 6. Export from posterior for graphing and other results ---------------------

js_3 <- jags.samples(jm, c("lambda_post", "p_post"), 50000, 50)

## Combine mcmc chains. 
post <- abind::abind(apply(js_3$lambda_post, 1:2, c), 
                     apply(js_3$p_post, 1:2, c),
                     along = 4) 

## Calculate BACI indicators: Adjust treatment and reference here !!!!!!!!!!!!!!
levels(d_rs$treatment)
eval <- c(1, 2, 4); ref <- 3
# eval <- c(2, 4); ref <- 1
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Calculate the BACI indicators for both responses and for all iterations:
BACI_out <- array(NA, dim = c(dim(post)[1], length(js_3), length(eval), 3))
for(m in 1:length(eval)){
  BACI_out[,,m,1] <- (post[,eval[m],,2] - post[,eval[m],,1]) -
                     (post[,ref,,2] - post[,ref,,1])        ## BACI
  BACI_out[,,m,2] <- abs(post[,eval[m],,2] - post[,eval[m],,1]) - 
                     abs(post[,ref,,2] - post[,ref,,1])     ## CI-ctrl
  BACI_out[,,m,3] <- abs(post[,eval[m],,2] - post[,ref,,2]) -     
                     abs(post[,eval[m],,1] - post[,ref,,1]) ## CI-div
}

## Keep track of the names:
dimnames(BACI_out) <- list(NULL,
                           names(js_3),
                           levels(d_rs$treatment)[eval],
                           c("BACI", "CI_ctr", "CI_div"))

## Calculate summary statistics:
BACI_rs <- apply(BACI_out, 2:4, posterior_summary)
BACI_rs <- dcast(melt(BACI_rs), Var2 + Var3 + Var4 ~ Var1)

## Add naming for figures later on:
BACI_rs$ref <- paste0("ref_", levels(d_rs$treatment)[ref])
colnames(BACI_rs)[1:3] <- c("response", "treatment", "indicator")

## Export BACI_rs and adjust name according to the chosen reference:
write.csv(BACI_rs, 
          paste0("clean/BACI_rs_ref_", 
                 ifelse(ref == 3, "NF", "CR"), 
                 ".csv"),
          row.names = FALSE)

## -------------------------------END-------------------------------------------
