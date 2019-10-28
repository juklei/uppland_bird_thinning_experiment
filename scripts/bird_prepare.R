## Create data sets to use for Bayesian models in JAGS defined in a different 
## script.
##
## First edit: 20190628
## Last edit: 20191023
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(reshape)

## 2. Load and explore data ----------------------------------------------------

dir("data")

occ <- read.csv("data/occ_double_2017to2019.csv")
ldm <- read.csv("data/long_distance_migrants.csv")

head(occ); head(ldm)

## 3. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 662. Instead of J indexed on site we want to index it on 
##    the species depending on wether it is long distance migrant or not.

## Exclude predators, birds with large hr and passers from obs:
occ <- occ[!occ$species %in% c("bergk",
                               "duvhk",
                               "ormvk",
                               "ekore",
                               "gravg",
                               "grona", 
                               "korp",
                               "mard",
                               "mindb",
                               "spark"), ]

## Chose relevant columns for this analysis:
occ <- droplevels(occ[, c("observer",
                          "block", 
                          "plot",
                          "visit", 
                          "sampling_period",
                          "obs_year",
                          "species")])

## Make data set so we have all possible combinations of all visits
## and all species seen during the whole survey.
b_occ <- expand.grid.df(unique(occ[, c("block", 
                                       "plot",
                                       "visit",
                                       "sampling_period",
                                       "obs_year")]), 
                        as.data.frame(unique(occ$species)))
colnames(b_occ)[length(b_occ)] <- "species"

## Merge with occ data set to aquire all information of observations
## If a speces was not observed during a visit make obs 0:

## Add obs = 1 to occ:
occ$observed <- 1

## Merge all b_occ with matching b_occ:
b_occ <- merge(b_occ, occ, all.x = TRUE, by = colnames(b_occ))

## NA in observed are actually non-observations, so they will become 0:
b_occ$observed[is.na(b_occ$observed)] <- 0

## All long distance migrants could be seen during at least 3 visits. Some
## were seen before, already during the 2nd one. In that case all 
## plots in block in which that species was seen get one more visit for that 
## species.

## Numerise the visit categories for simplifying calculation
b_occ$visit_num[b_occ$visit == "first"] <- 1
b_occ$visit_num[b_occ$visit == "second"] <- 2
b_occ$visit_num[b_occ$visit == "third"] <- 3
b_occ$visit_num[b_occ$visit == "fourth"] <- 4
b_occ$visit_num[b_occ$visit == "fifth"] <- 5

b_occ <- as.data.table(b_occ)
b_occ[, "n_visits" := ifelse(.SD$species %in% ldm$species, 
                             max(visit_num) - 2, 
                             max(visit_num)),
      by = c("plot", "obs_year")]

## Has the ldm species been seen in block i during the second visit?
T1 <- b_occ[b_occ$visit == "second" & b_occ$species %in% ldm$species, 
            list("observable_second" = ifelse(sum(observed) > 0, 1, 0),
                 "plot" = plot), 
            by = c("obs_year", "block", "species")]

## Reduce to actual observations:
T1 <- unique(T1[T1$observable_second == 1, c("obs_year", "plot", "species")])

## Merge b_occ with T1 to add 1 to rows in T1:
T1$add <- 1
b_occ <- merge(b_occ, T1, all.x = TRUE, by = c("obs_year", "plot", "species"))
b_occ$add[is.na(b_occ$add)] <- 0

## Add together and double n_visits for both sampling periods per visit:
b_occ$n_visits <- (b_occ$n_visits + b_occ$add)*2

## Count the number of times a species was seen per plot and year:
b_occ[, "n_obs" := sum(observed), by = c("obs_year", "plot", "species")]

## Reduce to needed columns:
b_occ <- unique(b_occ[b_occ$species != "no_obs", c("obs_year", 
                                                   "block", 
                                                   "plot", 
                                                   "species", 
                                                   "n_obs",
                                                   "n_visits")])

## Add the observer:
b_occ$observer <- ifelse(b_occ$block %in% c("ravsta", 
                                            "hagaberg", 
                                            "kungshamn", 
                                            "hallsbo" ),
                         "ses",
                         "jkn")

## Export:

dir.create("clean")
write.csv(b_occ, "clean/bpo_double.csv")

## -------------------------------END-------------------------------------------
