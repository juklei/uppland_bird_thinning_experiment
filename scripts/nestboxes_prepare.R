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

head(d_nb)

## Remove data from blocks which were visited only at one occasion after the
## experiment, because many of these nestboxes have an unknown owner:
d_nb <- d_nb[!(d_nb$block %in% c("sodersjon", "stocksatra", "kyrkstigen", 
                                 "stenby", "soderlejde", "halsingdal", 
                                 "brudgatan") & d_nb$year == 2019), ]

## Also plot 33 &3 6 in 2019, because they were forgotten:
d_nb <- d_nb[!(d_nb$plot %in% c("plot_33", "plot_36") & d_nb$year == 2019), ]

## 3. Prepare data for the JAGS module -----------------------------------------

d_nb <- d_nb[d_nb$species %in% c("talgoxe", "sv_flug", "blames", "empty"), 
             c("block", "plot", "year", "experiment", "treatment", "species")]
d_nb <- droplevels(d_nb)

## Response variable matrix with species as rows and observations as columns:
d_nb$id <- 1:nrow(d_nb)
occ <- acast(d_nb, id ~ species)
occ[] <- ifelse(is.na(occ[]), 0, 1)

## Compile data for JAGS:
data <- list(nobs = nrow(d_nb),
             occ = occ,
             block = as.numeric(d_nb$block),
             year = d_nb$year - 2016,
             treat = as.numeric(d_nb$treatment),
             exp = ifelse(d_nb$experiment == "before", 1, 2))

## 4. Prepare inits and run the model -------------------------------------------

## Prepare inits (Remember priors must not overlap, 
## i.e. define non-overlapping init):
i_alpha <- array(c(-5, 0, 5), 
                 dim = c((ncol(occ)-1), max(data$block), max(data$year)))
i_exp_effect <- array(c(-5, 0, 5), 
                      dim = c((ncol(occ)-1), max(data$treat), max(data$exp)))

inits <- list(list(alpha = i_alpha, exp_effect = i_exp_effect),
              list(alpha = i_alpha/5, exp_effect = i_exp_effect/5),
              list(alpha = i_alpha*2, exp_effect = i_exp_effect*2))

## Load the model:
model <- "scripts/JAGS/nestbox_JAGS_species.R"

## Run and burn in:
jm <- jags.model(model, data, inits, 3, 5000)
update(jm, 5000)

## Sample from the parameter posteriors:
cs_1 <- coda.samples(jm, c("alpha", "exp_effect"), 5000, 5)

## 5. Validate the model and export validation data and figures ----------------

pdf("figures/nestbox_species.pdf"); plot(cs_1); gelman.plot(cs_1); dev.off()

## 6. Export from posterior for graphing and other results ---------------------

cs_2 <- coda.samples(jm, "p_post", 5000, 5)

pdf("figures/nestbox_species_post.pdf"); plot(cs_2); dev.off()

## -------------------------------END-------------------------------------------
