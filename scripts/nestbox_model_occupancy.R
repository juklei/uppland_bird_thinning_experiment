## Analysis of nestbox occupancy before and after experiment.
## no eggs = empty = no reproductive attempt
##
## First edit: 20210130
## Last edit: 20210201
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
d_nb <- d_nb[d_nb$species %in% c("talgoxe", "sv_flug", "blames", "empty"), 
             c("plot", "year", "box", "experiment", "treatment", "species")]
d_nb <- droplevels(d_nb)

## Reduce d_pos to plot, box and coordinates and add to d_nb:
d_pos <- na.omit(d_pos[, c("plot", "box", "north", "east")])
d_nb <- merge(d_nb, d_pos, by = c("plot", "box"), all.x = TRUE)

## 3. Prepare data for the JAGS module -----------------------------------------

## Create a inverse distance matrix for calculating Moran's I for every year:
## Which year for Moran's I?
y_sel <- 2017
# y_sel <- 2018
# y_sel <- 2019
dm <- d_nb[d_nb$year == y_sel, c("north", "east")]
dm <- 1/as.matrix(dist(dm))
diag(dm) <- 0
seq <- which(d_nb$year == y_sel)

## Response variable matrix with species as rows and observations as columns:
d_nb$id <- 1:nrow(d_nb)
occ <- acast(d_nb, id ~ species)
occ[] <- ifelse(is.na(occ[]), 0, 1)
occ <- occ[, c(1, 3, 4, 2)] ## "empty" as last category

## Compile data for JAGS:
data <- list(nobs = nrow(d_nb),
             occ = occ,
             year = d_nb$year,
             treat = as.numeric(d_nb$treatment),
             exp = ifelse(d_nb$experiment == "before", 1, 2),
             dm = dm,
             seq = seq)

## 4. Prepare inits and run the model -------------------------------------------

## Prepare inits (Remember priors must not overlap, 
## i.e. define non-overlapping inits):
i_y_effect <- array(c(-5, 0, 5), dim = c((ncol(occ)-1), 2))
i_exp_effect <- array(c(-5, 0, 5), 
                      dim = c((ncol(occ)-1), max(data$treat), max(data$exp)))

inits <- list(list(y_effect = i_y_effect, exp_effect = i_exp_effect),
              list(y_effect = i_y_effect/5, exp_effect = i_exp_effect/5),
              list(y_effect = i_y_effect*2, exp_effect = i_exp_effect*2))

## Load the model:
model <- "scripts/JAGS/nestbox_JAGS_species.R"

## Run and burn in:
jm <- jags.model(model, data, inits, 3, 5000)
update(jm, 5000)

## Sample from the parameter posteriors:
cs_1 <- coda.samples(jm, c("y_effect", "exp_effect"), 50000, 50)

## 5. Validate the model and export validation data and figures ----------------

## Check mixing and convergence:
pdf("figures/nestbox_species.pdf"); plot(cs_1); gelman.plot(cs_1); dev.off()

## Produce validation metrics for posterior predictive checks:
js_2 <- jags.samples(jm, c("mean_obs", "mean_sim",
                           "sd_obs", "sd_sim",
                           "fit", "fit_sim",
                           "I"), 5000, 5)

## Calculate Moran's I if random spatial distribution of residuals: 
E0 <- round(-1/(nrow(dm)*3 - 1), 4)

## Export validations as figure:
pdf("figures/nestbox_species_ppc&MI.pdf")
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

js_3 <- jags.samples(jm, "p_post", 50000, 50)
p_post <- apply(js_3$p_post, 1:3, c) ## Combine mcmc chains. 

## Calculate BACI indicators: Adjust treatment and reference here !!!!!!!!!!!!!!
levels(d_nb$treatment)
# eval <- c(1, 2, 4); ref <- 3
eval <- c(2, 4); ref <- 1
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Calculate the BACI indicators for all species and for all iterations: -------
BACI_out <- array(NA, dim = c(dim(p_post)[1:2], length(eval), 3))
for(m in 1:length(eval)){
  BACI_out[,,m,1] <- (p_post[,,eval[m],2] - p_post[,,eval[m],1]) -
                     (p_post[,,ref,2] - p_post[,,ref,1])        ## BACI
  BACI_out[,,m,2] <- abs(p_post[,,eval[m],2] - p_post[,,eval[m],1]) - 
                     abs(p_post[,,ref,2] - p_post[,,ref,1])     ## CI-ctrl
  BACI_out[,,m,3] <- abs(p_post[,,eval[m],2] - p_post[,,ref,2]) -     
                     abs(p_post[,,eval[m],1] - p_post[,,ref,1]) ## CI-div
}

## Keep track of the names:
dimnames(BACI_out) <- list(NULL,
                           colnames(occ),
                           levels(d_nb$treatment)[eval],
                           c("BACI", "CI_ctr", "CI_div"))

## Calculate summary statistics:
BACI_nb <- apply(BACI_out, 2:4, posterior_summary)
BACI_nb <- dcast(melt(BACI_nb), Var2 + Var3 + Var4 ~ Var1)

## Add naming for figures later on:
BACI_nb$ref <- paste0("ref_", levels(d_nb$treatment)[ref])
colnames(BACI_nb)[1:3] <- c("species", "treatment", "indicator")

## Export BACI_nb and adjust name according to the chosen reference:
write.csv(BACI_nb, 
          paste0("clean/BACI_nb_occ_ref_", 
                 ifelse(ref == 3, "NF", "CR"), 
                 ".csv"),
          row.names = FALSE)

## -------------------------------END-------------------------------------------
