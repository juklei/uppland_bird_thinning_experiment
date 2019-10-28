## Take the forest data and create all sorts of explamnatory variables related 
## to the thinning experiment

## First edit: 20191023
## Last edit:  20191028

## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

library(data.table)
library(reshape)

## 2. Load and explore data ----------------------------------------------------

dir("data")
forest <- read.csv("data/forest_data_uppland_plot.csv")

head(forest)

## 3. Reduce data set to needed variables, claculate percentages, 
##    and calculate differences due to treatments.

## Exclude unneeded varables:
f_red <- as.data.table(forest[,c(1:3,19:length(forest))])

## Calculate percentages:
f_red$perc_gran <- f_red$nr_gran/f_red$nr_all_alive
f_red$perc_lov <- f_red$nr_lov/f_red$nr_all_alive
f_red$perc_tall <- f_red$nr_tall/f_red$nr_all_alive
f_red$perc_skarm <- f_red$nr_skarm/f_red$nr_all_alive

## Calculate differences:

## Which plots have after data?
B1 <- f_red$plot[f_red$experiment == "after"]
f_diff <- f_red[f_red$plot %in% B1,]

## Calculate the difference:
diff <- f_diff[f_diff$experiment == "after",] - 
        f_diff[f_diff$experiment == "before",]

## Replace name variables in the new file with the unique values from above:
diff[,c(1:2,16)] <- unique(f_diff[,c(1:2,16)])
diff$experiment <- as.character(diff$experiment)
diff$experiment <- "difference"

## Join before-after with differences:
f_all <- rbind(f_red,diff)

## 4. Export:

dir.create("clean")
write.csv(f_all, "clean/forest_experiment_data.csv", row.names = F)

## -----------------------------------END---------------------------------------