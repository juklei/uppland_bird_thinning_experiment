## Create data sets to use for Bayesian models in JAGS defined in a different 
## script.
##
## First edit: 20190628
## Last edit: 20191030
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(reshape)

## 2. Load and explore data ----------------------------------------------------

dir("data")

OCC <- read.csv("data/occ_double_2017to2019.csv")
OCC2 <- read.csv("data/occ_2016to2019.csv")
ldm <- as.vector(na.omit(read.csv("data/bird_data.csv")$ldm))

head(occ); head(ldm)

## 3. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 662. Instead of J indexed on site we want to index it on 
##    the species depending on wether it is long distance migrant or not.

## Its not justified to use plot 30, 118, 120 and 121 as true controls as they
## are located too close to the thinning stands:
OCC <- OCC[!OCC$plot %in% c("plot_30", "plot_118", "plot_120", "plot_121"), ]

## Exclude predators, birds with large hr and passers from obs:
occ <- OCC[!OCC$species %in% c("bergk", "duvhk", "ormvk", "ekore", 
                               "gravg", "grona", "korp", "mard", 
                               "mindb", "spark", "tjadr", "morka", 
                               "spila", "gok", "grong", "jarpe",
                               "kraka", "stort", "ringa", "notsa", "notka"), ]

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
b_occ[, "n_visits" := ifelse(.SD$species %in% ldm, 
                             max(visit_num) - 2, 
                             max(visit_num)),
      by = c("plot", "obs_year")]

## Has the ldm species been seen in block i during the second visit?
T1 <- b_occ[b_occ$visit == "second" & b_occ$species %in% ldm, 
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

## 4. Rearrange bird observation data to fit the model presented in the --------
##    AHM book on p. 690. 

## Exclude predators, birds with large hr and passers from obs:
occ <- droplevels(OCC[!OCC$species %in% c("bergk", "duvhk", "ormvk", "ekore", 
                                          "gravg", "grona", "korp", "mard", 
                                          "mindb", "spark", "tjadr", "morka", 
                                          "spila", "gok", "grong", "jarpe",
                                          "kraka", "stort", "ringa", "notsa", "notka"), ])

## Make data set so we have all possible combinations of all visits
## and all species seen during the whole survey.
b_occ <- expand.grid.df(unique(occ[, c("block", 
                                       "plot",
                                       "visit",
                                       "sampling_period",
                                       "obs_year",
                                       "observer")]), 
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

## Add arrival date:
b_occ <- as.data.table(b_occ)
b_occ$dp_march <- as.numeric(b_occ$dp_march) ## min() had troubles handling integers
b_occ[, "min_dpm" := ifelse(all(is.na(dp_march)), 99, min(dp_march, na.rm = T)), 
      by = c("species", "obs_year")] ## 99 for never observed in a year

## Fill NAs in dp_march and min_post_sunrise:

## Create a function that does that:
F1 <- function(x){
  x$dp_march[is.na(x$dp_march)] <- unique(x$dp_march[!is.na(x$dp_march)])
  x$min_post_sunrise[is.na(x$min_post_sunrise)] <- 
    mean(x$min_post_sunrise[!is.na(x$min_post_sunrise)])
  return(x)
}

## Apply it to all combinatins of plot*obs_year*visit*sampling period:
b_occ <- b_occ[, F1(.SD), by = c("plot", "visit", "sampling_period", "obs_year")]

## Remove ldm non-observations from the first & second visit that have a 
## dp_march smaller than the arrival date (e.g. bird has not arrived yet):
  
b_occ <- b_occ[!(b_occ$species %in% ldm & 
                   b_occ$visit %in% c("first", "second") & 
                   b_occ$dp_march < b_occ$min_dpm), ]

## Remove no_obs from species and select columns to export:
b_occ <- b_occ[b_occ$species != "no_obs", c("obs_year", "block", "observer", 
                                            "plot", "species", "observed", 
                                            "dp_march", "min_post_sunrise")]

## Add a visit indicator for each visit per plot*obs_year*species
b_occ <- b_occ[, "visit" := 1:nrow(.SD), by = c("plot", "obs_year", "species")]

## Export:

dir.create("clean")
write.csv(b_occ, "clean/bpo_double_bernoulli.csv", row.names = FALSE)

## -------------------------------END-------------------------------------------
