## Analysis of nestbox occupancy before and after experiment.
## no eggs = empty = no reproductive attempt
##
## First edit: 20210130
## Last edit: 20210201
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(boot)
library(rjags)
library(runjags)
library(coda)
library(reshape2)

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

## Also plot 33 &3 6 in 2019, because they were forgotten:
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
             year = d_nb$year - 2016,
             treat = as.numeric(d_nb$treatment),
             exp = ifelse(d_nb$experiment == "before", 1, 2),
             dm = dm,
             seq = seq)

## 4. Prepare inits and run the model -------------------------------------------

## Prepare inits (Remember priors must not overlap, 
## i.e. define non-overlapping init):
i_alpha <- array(c(-5, 0, 5), dim = c((ncol(occ)-1), max(data$year)))
i_exp_effect <- array(c(-5, 0, 5), 
                      dim = c((ncol(occ)-1), max(data$treat), max(data$exp)))

inits <- list(list(alpha = i_alpha, exp_effect = i_exp_effect),
              list(alpha = i_alpha/5, exp_effect = i_exp_effect/5),
              list(alpha = i_alpha*2, exp_effect = i_exp_effect*2))

## Load the model:
model <- "scripts/JAGS/nestbox_JAGS_species.R"

## Run and burn in:
jm <- jags.model(model, data, inits, 3, 5000)
update(jm, 15000)

## Sample from the parameter posteriors:
cs_1 <- coda.samples(jm, c("alpha", "exp_effect"), 10000, 10)

## 5. Validate the model and export validation data and figures ----------------

## Check mixing and convergence:
pdf("figures/nestbox_species.pdf"); plot(cs_1); gelman.plot(cs_1); dev.off()

## Produce validation metrics for posterior predictive checks:
cs_2 <- jags.samples(jm, c("mean_obs", "mean_sim",
                           "cv_obs", "cv_sim",
                           "fit", "fit_sim",
                           "I"), 10000, 10)

## Calculate Moran's I if random spatial distribution of residuals: 
E0 <- round(-1/(nrow(dm)*3 - 1), 4)

## Export validations as figure:
pdf("figures/nestbox_species_ppc&MI.pdf")
par(mfrow = c(2, 2))
plot(cs_2$mean_obs, cs_2$mean_sim, xlab = "mean real", ylab = "mean simulated")
abline(0, 1)
plot(cs_2$cv_obs, cs_2$cv_sim, xlab = "cv real", ylab = "cv simulated")
abline(0,1)
plot(cs_2$fit, cs_2$fit_sim, xlab = "ssq real", ylab = "ssq simulated")
abline(0,1)
plot(density(cs_2$I), main = paste0("Moran's I, ", y_sel, ", E0 = ", E0))
dev.off()

## 6. Export from posterior for graphing and other results ---------------------

cs_3 <- coda.samples(jm, "p_post", 10000, 10)
pdf("figures/nestbox_species_post.pdf"); plot(cs_3); dev.off()

## -------------------------------END-------------------------------------------
