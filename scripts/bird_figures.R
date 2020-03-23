## Make figures for the experimental data and the bpo model responses
##
## First edit: 20191028
## Last edit:  20200226
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")
require("wesanderson")

## 2. Load and explore data ----------------------------------------------------

bird_names <- read.csv("data/birdnames.csv")
head(bird_names)

non_exp <- read.csv("clean/non_exp.csv")
non_exp <- merge(non_exp, 
                 bird_names, 
                 all.x = TRUE, 
                 by.x = "identity", 
                 by.y = "short")
head(non_exp)

BACI_sl <- read.csv("clean/BACI_sl_treatment_TC.csv")
BACI_sl <- merge(BACI_sl, 
                 bird_names, 
                 all.x = TRUE, 
                 by.x = "identity", 
                 by.y = "short")
levels(BACI_sl$treatment) <- c("Complete retention", 
                               "Conventional thinning", 
                               "Understory retention thinning")
BACI_sl$treatment <- factor(BACI_sl$treatment, 
                            levels = c("Complete retention", 
                                       "Understory retention thinning",
                                       "Conventional thinning"))
levels(BACI_sl$indicator)[2:3] <- c("CI-contribution", "CI-divergence")
head(BACI_sl)

## 3. Make graphs for the non-experimental plots for forest variables ----------

## Add significance grouping:
non_exp$cross <- ifelse(sign(non_exp$X2.5.) == sign(non_exp$X97.5.),
                        "95% CI does not cross zero", 
                        "95% CI crosses zero")

## Chose for which variable you want to make the figure:
var_choice <- "BA_dw"
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

## 4. Make graphs for bpo predictions ------------------------------------------

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$identity %in% c("cm", "alpha", "beta"), 
                      as.character(BACI_sl$identity),
                      "species")

## Graph for slopes with CIs:

## Order:
O1 <- BACI_sl[BACI_sl$indicator == "BACI" & 
                BACI_sl$treatment == "Complete retention", 
              c("X50.", "long")]
O1 <- O1$long[order(O1$X50.)]
BACI_sl$long <- factor(BACI_sl$long, levels = O1)

## Add significance grouping:
BACI_sl$cross <- ifelse(sign(BACI_sl$X2.5.) == sign(BACI_sl$X97.5.),
                        "95% CI does not cross zero", 
                        "95% CI crosses zero")

g1 <- ggplot(data = droplevels(BACI_sl[BACI_sl$cat %in% c("cm", "species"), ]), 
             aes(x = long, y = X50., 
                 colour = treatment, 
                 fill = treatment))#,
                 # linetype = cross))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 3, 
                    width = 0, 
                    position = position_dodge(0.5))
g3 <- geom_point(position = position_dodge(0.5), size = 4, colour = "black")
g4 <- facet_grid(cat ~ indicator, space = "free", scales = "free")
G <- g1 +
  geom_hline(yintercept = 0, size = 2) + 
  g2 + g4 + g3 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        # legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png("figures/BACI_slopes_TC_2.png", 21000/8, 24000/8, "px", res = 600/8)
G
dev.off()

## Graph for probabilies:

## Order:
O2 <- BACI_sl[BACI_sl$indicator == "BACI" & 
                BACI_sl$treatment == "Complete retention", 
              c("ecdf", "long")]
O2 <- O2$long[order(O2$ecdf)]
BACI_sl$long <- factor(BACI_sl$long, levels = O2)

BACI_sl$P_ecdf <- ifelse(BACI_sl$ecdf < 0.11 | BACI_sl$ecdf > 0.89,
                         "P > .89",
                         "P < .89")

p1 <- ggplot(data = BACI_sl, 
             aes(x = long, 
                 y = 0, 
                 colour = treatment))#, 
                 # linetype = P_ecdf))
p2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_sl$ecdf), 
                     position = position_dodge(0.5),
                     size = 3.1, 
                     width = 0)
p2b <- geom_errorbar(aes(ymin = BACI_sl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 3.1,                  
                     width = 0)
P <- p1 +
  geom_hline(yintercept = 0, size = 2) + 
  # scale_linetype_manual(values=c("solid", "dotted")) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1.00"))) +
  p2a + p2b + g4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  # scale_colour_manual(values = wes_palette(n = 3, name = "Rushmore1")) +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        # legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png("figures/BACI_probs_TC_2.png", 21000/8, 25000/8, "px", res = 600/8)
P
dev.off()

## -------------------------------END-------------------------------------------
