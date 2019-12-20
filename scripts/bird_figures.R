## Make figures for the experimental data and the bpo model responses
##
## First edit: 20191028
## Last edit:  20191211
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Load and explore data ----------------------------------------------------

forest <- as.data.table(read.csv("clean/forest_experiment_data.csv"))
head(forest)

non_exp <- read.csv("clean/non_exp.csv")
head(non_exp)

BACI_sl <- read.csv("clean/BACI_sl_BA_all_C.csv")
head(BACI_sl)

## 3. Make graphs with the explanatory data ------------------------------------

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

## 4. Make graphs for the non-experimental plots for forest variables ----------

## Add significance grouping:
non_exp$cross <- ifelse(sign(non_exp$X2.5.) == sign(non_exp$X97.5.),
                        "95% CI does not cross zero", 
                        "95% CI crosses zero")

## Chose for which variable you want to make the figure:
var_choice <- "vis"
non_exp_red <- non_exp[non_exp$variable == var_choice, ]

## Graph for slopes with CIs:

## Order:
non_exp_red$identity <- factor(non_exp_red$identity, 
                               non_exp_red$identity[order(non_exp_red$X50.)])

q1 <- ggplot(data = droplevels(non_exp_red[non_exp_red$identity != "cm", ]), 
             aes(x = identity, y = X50., color = cross))
q2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 5, 
                    width = 0, 
                    position = position_dodge(0.5))
q3 <- geom_point(position = position_dodge(1), size = 8, color = "black")
q4 <- geom_hline(yintercept = unlist(non_exp_red[non_exp_red$identity == "cm", 
                                                 c("X2.5.", "X50.", "X97.5.")]), 
                 size = 2, 
                 linetype = c("dashed", "solid", "dashed"),
                 color = "grey")
Q <- q1 +
     geom_hline(yintercept = 0, size = 2) +
     q4 + q2 + q3 + 
     xlab("") + ylab("") + theme_classic(40) + coord_flip() +
     theme(legend.position = "top", #c(0.7, 0.1),
           legend.title = element_blank(),
           legend.key.size = unit(3, 'lines'),
           legend.direction = "horizontal",
           strip.text.y = element_blank())

png(paste0("figures/non_exp_slopes_", var_choice, ".png"), 
           10000/8, 15000/8, 
           "px", 
           res = 600/8)
Q
dev.off()

## 5. Make graphs for bpo predictions ------------------------------------------

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$identity %in% c("cm", "alpha", "beta"), 
                      as.character(BACI_sl$identity),
                      "species")

## Graph for slopes with CIs:

## Order:
O1 <- BACI_sl[BACI_sl$indicator == "BACI" & BACI_sl$treatment == "below", 
              c("X50.", "identity")]
O1 <- O1$identity[order(O1$X50.)]
BACI_sl$identity <- factor(BACI_sl$identity, levels = O1)

## Add significance grouping:
BACI_sl$cross <- ifelse(sign(BACI_sl$X2.5.) == sign(BACI_sl$X97.5.),
                        "95% CI does not cross zero", 
                        "95% CI crosses zero")

g1 <- ggplot(data = droplevels(BACI_sl[BACI_sl$cat %in% c("cm", "species"), ]), 
             aes(x = identity, y = X50., 
                 colour = treatment, 
                 fill = treatment,
                 linetype = cross))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), size = 3, width = 0, position = position_dodge(0.5))
#g3 <- geom_point(position = position_dodge(1), size = 1, color = "black")
g4 <- facet_grid(cat ~ indicator, space = "free", scales = "free")
G <- g1 +
  geom_hline(yintercept = 0, size = 2) + 
  g2 + g4 + #g3 +
  xlab("") + ylab("") +
  theme_bw(40) + coord_flip() +
  theme(legend.position = "top", #c(0.7, 0.1),
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal",
        strip.text.y = element_blank())

png("figures/BACI_slopes_BA_all_C.png", 20000/8, 20000/8, "px", res = 600/8)
G
dev.off()

## Graph for probabilies:

## Order:
O2 <- BACI_sl[BACI_sl$indicator == "BACI" & BACI_sl$treatment == "below", 
              c("ecdf", "identity")]
O2 <- O2$identity[order(O2$ecdf)]
BACI_sl$identity <- factor(BACI_sl$identity, levels = O2)

BACI_sl$P_ecdf <- ifelse(BACI_sl$ecdf < 0.1 | BACI_sl$ecdf > 0.9,
                         "P > .90",
                         "P < .90")

p1 <- ggplot(data = BACI_sl, 
             aes(x = identity, 
                 y = 0, 
                 colour = treatment, 
                 linetype = P_ecdf))
p2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_sl$ecdf), 
                     position = position_dodge(0.5),
                     size = 3, 
                     width = 0)
p2b <- geom_errorbar(aes(ymin = BACI_sl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 3,                  
                     width = 0)
P <- p1 +
  geom_hline(yintercept = 0, size = 2) + 
  scale_y_continuous(breaks = c(-0.95, -0.5, 0, 0.5, 0.95),
                     labels = c(".95", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", ".95"))) +
  p2a + p2b + g4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  theme_bw(40) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.direction = "horizontal",
        strip.text.y = element_blank())

png("figures/BACI_probs_BA_all_C.png", 20000/8, 20000/8, "px", res = 600/8)
P
dev.off()

## -------------------------------END-------------------------------------------
