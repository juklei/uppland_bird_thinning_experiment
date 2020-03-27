## Make figures for the experimental data and the bpo model responses
##
## First edit: 20191028
## Last edit:  20200325
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Load and prepare data ----------------------------------------------------

bird_data <- read.csv("data/bird_data.csv")
BACI_sl <- read.csv("clean/BACI_sl.csv")
BACI_gl <- read.csv("clean/BACI_gl.csv")

## 3. Make graphs for species level --------------------------------------------

## Prepare data:
BACI_sl <- merge(BACI_sl, 
                 bird_data[, c("short", "long")], 
                 all.x = TRUE, 
                 by.x = "identity", 
                 by.y = "short")
levels(BACI_sl$long) <- c(levels(BACI_sl$long), "Community mean") 
BACI_sl$long[is.na(BACI_sl$long)] <- "Community mean"
levels(BACI_sl$treatment) <- c("Complete retention", 
                               "Conventional thinning", 
                               "Understory retention thinning")
BACI_sl$treatment <- factor(BACI_sl$treatment, 
                            levels = c("Complete retention", 
                                       "Understory retention thinning",
                                       "Conventional thinning"))
levels(BACI_sl$indicator)[2:3] <- c("CI-contribution", "CI-divergence")

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$identity == "cm", "cm", "species")

head(BACI_sl)

## Graph for slopes with CIs:

## Order:
O1 <- BACI_sl[BACI_sl$indicator == "BACI" & 
                BACI_sl$treatment == "Complete retention", 
              c("X50.", "long")]
O1 <- O1$long[order(O1$X50.)]
BACI_sl$long <- factor(BACI_sl$long, levels = O1)

## Make figure:
g1 <- ggplot(data = BACI_sl, 
             aes(x = long, y = X50., colour = treatment, fill = treatment))
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
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png("figures/BACI_sl_slopes.png", 21000/8, 24000/8, "px", res = 600/8)
G
dev.off()

## Graph for probabilies:

## Order:
O2 <- BACI_sl[BACI_sl$indicator == "BACI" & 
                BACI_sl$treatment == "Complete retention", 
              c("ecdf", "long")]
O2 <- O2$long[order(O2$ecdf)]
BACI_sl$long <- factor(BACI_sl$long, levels = O2)

## Make graph:
p1 <- ggplot(data = BACI_sl, aes(x = long, y = 0, colour = treatment))
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
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  p2a + p2b + g4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png("figures/BACI_sl_probs.png", 21000/8, 24000/8, "px", res = 600/8)
P
dev.off()

## 4. Make graphs for guilds & diversity ---------------------------------------

## Prepare the data:
levels(BACI_gl$treatment) <- c("Complete retention", 
                               "Conventional thinning", 
                               "Understory retention thinning")
BACI_gl$treatment <- factor(BACI_gl$treatment, 
                               levels = c("Complete retention", 
                                          "Understory retention thinning",
                                          "Conventional thinning"))
levels(BACI_gl$indicator)[2:3] <- c("CI-contribution", "CI-divergence")
BACI_gl$cat[BACI_gl$identity %in% c("bark", 
                                    "f_cpy", 
                                    "f_grd", 
                                    "grd_cpy")] <- "Foraging"
BACI_gl$cat[BACI_gl$identity %in% c("hole", "n_cpy", "n_grd")] <- "Nesting"
BACI_gl$cat[BACI_gl$identity %in% c("insect", "omni")] <- "Food"
BACI_gl$cat[BACI_gl$identity %in% c("bd", "r")] <- "Diversity"
levels(BACI_gl$identity) <- c("Bark feeder", "Beta (negJaccard)", 
                              "Canopy feeder", "Ground feeder", 
                              "Ground/Canopy feeder", "Hole nester", 
                              "Insectivore", "Canopy nester", "Ground nester", 
                              "Omnivore", "Alpha")
head(BACI_gl)

## Make figure:
h1 <- ggplot(data = droplevels(BACI_gl[BACI_gl$cat != "Diversity", ]), 
             aes(x = identity, y = X50., colour = treatment, fill = treatment))
h2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 5, 
                    width = 0, 
                    position = position_dodge(0.5))
h3 <- geom_point(position = position_dodge(0.5), size = 6.5, colour = "black")
h4 <- facet_grid(cat ~ indicator, space = "free", scales = "free")
H <- h1 +
  geom_hline(yintercept = 0, size = 2) + 
  h2 + h3 + h4 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"))

png("figures/BACI_gl_slopes.png", 21000/8, 13000/8, "px", res = 600/8)
H
dev.off()

## Graph for probabilies:

## Make graph:
q1 <- ggplot(data = BACI_gl, aes(x = identity, y = 0, colour = treatment))
q2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_gl$ecdf), 
                     position = position_dodge(0.5),
                     size = 5, 
                     width = 0)
q2b <- geom_errorbar(aes(ymin = BACI_gl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 5,                  
                     width = 0)
Q <- q1 +
  geom_hline(yintercept = 0, size = 2) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  q2a + q2b + h4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"))

png("figures/BACI_gl_probs.png", 21000/8, 16000/8, "px", res = 600/8)
Q
dev.off()

## -------------------------------END-------------------------------------------
