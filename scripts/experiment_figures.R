## Make figures for the experimental data to describe the data and the results
##
## First edit: 20191028
## Last edit:  20191028
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Define or source functions used in this script ---------------------------


## 3. Load and explore data ----------------------------------------------------

forest <- as.data.table(read.csv("clean/forest_experiment_data.csv"))
head(forest)

## 4. Make graphs with the explanatory data ------------------------------------

## Graph with tree numbers:

colnames(forest)

G1data <- melt(forest[,c(2:8,14,16)]) ## tree_nr
# G1data <- melt(forest[,c(3,9:13,16)]) ## dbh
# G1data <- melt(forest[,c(3,15,16)]) ## laser
# G1data <- melt(forest[,c(3,17:20,16)]) ## perc

G1data <- G1data[#forest$experiment != "difference" & 
                 G1data$treatment %in% c("T","URT"),]

G1summary <- G1data[, list("mean"=mean(value, na.rm=T), 
                           "sd"=sd(value, na.rm=T)),
                    by=c("experiment","treatment","variable")]

G1summary$experiment <- factor(G1summary$experiment, 
                               levels=c("before","after","difference"))

G1 <- ggplot(G1summary, aes(x=treatment, y=mean, fill=experiment)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                    position=position_dodge(.9)) +
      facet_wrap(. ~ variable, scales="free_y") + 
      theme_bw(20) +
      theme(legend.position = "top", 
            legend.direction = "horizontal",
            legend.title = element_blank())

png("figures/perc_tree.png", 10000/4, 7000/4, "px", res = 600/4)

G1

dev.off()

## 5. ... ----------------------------------------------------------------------


## -------------------------------END-------------------------------------------
