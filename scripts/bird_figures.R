## Make figures for the experimental data and the bpo model responses
##
## First edit: 20191028
## Last edit:  20200428
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")

## 2. Load and prepare data ----------------------------------------------------

bird_data <- read.csv("data/bird_data.csv")
BACI_sl_ctrl <- read.csv("clean/BACI_sl_ref_control.csv")
BACI_gl_ctrl <- read.csv("clean/BACI_gl_ref_control.csv")
BACI_sl_CR <- read.csv("clean/BACI_sl_ref_CR.csv")
BACI_gl_CR <- read.csv("clean/BACI_gl_ref_CR.csv")

## Chose reference category here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ref <- "control"
ref <- "CR"
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## 3. Make graphs for species level --------------------------------------------

BACI_sl <- rbind(BACI_sl_ctrl, BACI_sl_CR)

## Prepare data:
levels(BACI_sl$treatment) <- c("Complete retention", 
                               "Conventional thinning", 
                               "Understory retention thinning")
BACI_sl$treatment <- factor(BACI_sl$treatment,
                            levels = c("Complete retention", 
                                       "Understory retention thinning",
                                       "Conventional thinning"))
levels(BACI_sl$indicator) <- c("BACI-contrast", 
                               "CI-contribution", 
                               "CI-divergence")

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$species == "Community mean", "cm", "species")

## Reduce data set to chosen reference level:
if(ref == "control"){
  gsl <- BACI_sl[BACI_sl$treatment != "Understory retention thinning" &
                   BACI_sl$ref == "ref_TC", ]
} else{
  gsl <- BACI_sl[BACI_sl$ref == "ref_C", ]
}
gsl <- droplevels(gsl)

## Graph for slopes with CIs:

## Order:
O1 <- gsl[gsl$indicator == "BACI-contrast" & 
          gsl$treatment == "Conventional thinning", 
          c("X50.", "species")]
O1 <- O1$species[order(O1$X50.)]
gsl$species <- factor(gsl$species, levels = O1)

## Make figure:
g1 <- ggplot(gsl, aes(species, X50., colour = treatment, fill = treatment))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 4, 
                    width = 0, 
                    position = position_dodge(0.5))
g3 <- geom_point(position = position_dodge(0.5), size = 5, colour = "black")
g4 <- facet_grid(cat ~ indicator, space = "free", scales = "free")
G <- g1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") + 
  g2 + g4 + g3 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png(paste0("figures/BACI_sl_", ref, "_slopes.png"),
    21000/8, 24000/8, 
    "px", 
    res = 600/8)
G
dev.off()

## Graph for probabilies:

## Order:
O2 <- gsl[gsl$indicator == "BACI-contrast" & 
          gsl$treatment == "Conventional thinning", 
          c("ecdf", "species")]
O2 <- O2$species[order(O2$ecdf)]
gsl$species <- factor(gsl$species, levels = O2)

## Make graph:
p1 <- ggplot(gsl, aes(species, 0, colour = treatment))
p2a <- geom_errorbar(aes(ymin = 0, ymax = gsl$ecdf), 
                     position = position_dodge(0.5),
                     size = 4, 
                     width = 0)
p2b <- geom_errorbar(aes(ymin = gsl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 4,                  
                     width = 0)
P <- p1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  p2a + p2b + g4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png(paste0("figures/BACI_sl_", ref, "_probs.png"),
    21000/8, 24000/8,
    "px", 
    res = 600/8)
P
dev.off()

## 4. Make graphs for guilds & trends ------------------------------------------

BACI_gl <- rbind(BACI_gl_ctrl, BACI_gl_CR)

## Prepare the data:
levels(BACI_gl$treatment) <- c("Complete retention", 
                               "Conventional thinning", 
                               "Understory retention thinning")
BACI_gl$treatment <- factor(BACI_gl$treatment, 
                            levels = c("Complete retention", 
                                       "Understory retention thinning",
                                       "Conventional thinning"))
levels(BACI_gl$indicator) <- c("BACI-contrast", 
                               "CI-contribution", 
                               "CI-divergence")

## Reduce data set to chosen reference level:
if(ref == "control"){
  ggl <- BACI_gl[BACI_gl$treatment != "Understory retention thinning" &
                   BACI_gl$ref == "ref_TC", ]
} else{
  ggl <- BACI_gl[BACI_gl$ref == "ref_C", ]
}
ggl <- droplevels(ggl)

## Make figure:
h1 <- ggplot(ggl, aes(guild, X50., colour = treatment, fill = treatment))
h2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 5, 
                    width = 0, 
                    position = position_dodge(0.5))
h3 <- geom_point(position = position_dodge(0.5), size = 6.5, colour = "black")
h4 <- facet_grid(group ~ indicator, space = "free", scales = "free")
H <- h1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") + 
  h2 + h3 + h4 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"))

png(paste0("figures/BACI_gl_", ref, "_slopes.png"),
    21000/8, 18000/8, 
    "px", 
    res = 600/8)
H
dev.off()

## Graph for probabilies:

## Make graph:
q1 <- ggplot(data = ggl, aes(x = guild, y = 0, colour = treatment))
q2a <- geom_errorbar(aes(ymin = 0, ymax = ggl$ecdf), 
                     position = position_dodge(0.5),
                     size = 5, 
                     width = 0)
q2b <- geom_errorbar(aes(ymin = ggl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 5,                  
                     width = 0)
Q <- q1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  q2a + q2b + h4 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(58) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"))

png(paste0("figures/BACI_gl_", ref, "_probs.png"), 
    21000/8, 18000/8, 
    "px", 
    res = 600/8)
Q
dev.off()

## -------------------------------END-------------------------------------------
