## Make figures for the experimental data and the bpo model responses
##
## First edit: 20191028
## Last edit:  20200511
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

rm(list = ls())

require("ggplot2")
require("data.table")
# devtools::install_github("teunbrand/ggh4x")
require("ggh4x")

## 2. Load and prepare data ----------------------------------------------------

bird_data <- read.csv("data/bird_data.csv")
BACI_sl_NF <- read.csv("clean/BACI_sl_ref_NF.csv")
BACI_gl_NF <- read.csv("clean/BACI_gl_ref_NF.csv")
BACI_sl_CR <- read.csv("clean/BACI_sl_ref_CR.csv")
BACI_gl_CR <- read.csv("clean/BACI_gl_ref_CR.csv")

## Exclude URT from NF controls:
BACI_sl_NF <- droplevels(BACI_sl_NF[BACI_sl_NF$treatment != "URT", ])
BACI_gl_NF <- droplevels(BACI_gl_NF[BACI_gl_NF$treatment != "URT", ])

## Combine all data:
BACI_sl <- rbind(BACI_sl_NF, BACI_sl_CR)
BACI_gl <- rbind(BACI_gl_NF, BACI_gl_CR)

## Change level names:
l1 <- c("Complete retention", 
        "Conventional thinning", 
        "Understory retention thinning")
levels(BACI_sl$treatment) <-  l1
levels(BACI_gl$treatment) <-  l1
l2 <- c("BACI-contrast", "CI-contribution", "CI-divergence")
levels(BACI_sl$indicator) <- l2
levels(BACI_gl$indicator) <- l2
l3 <- c("Control = No forestry", "Control = Complete retention")
levels(BACI_sl$ref) <- l3
levels(BACI_gl$ref) <- l3

## 3. Make combined figure for BACI_sl -----------------------------------------

## Graph with slopes:

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$species == "Community mean", "cm", "species")

## Start species names with "A" on top:
BACI_sl$species <- factor(BACI_sl$species, rev(levels(BACI_sl$species)))

## Make figure:
g1 <- ggplot(BACI_sl, aes(species, X50., colour = treatment, fill = treatment))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 4, 
                    width = 0, 
                    # alpha = 0.8,
                    position = position_dodge(0.6))
g3 <- geom_point(position = position_dodge(0.6), size = 6, colour = "black")
g4 <- facet_nested(vars(cat), vars(ref, indicator), "free", scales = "free")
G <- g1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") + 
  g2 + g4 + g3 +
  xlab("") + ylab("Indicator value for probability of occurrence") + 
  coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#FC4E07", "#E7B800")) +
  theme_light(70) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(5, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white", size = 0.8))

png("figures/BACI_sl_slopes.png", 30000/8, 28000/8, "px", res = 600/8)
G
dev.off()

# ## Graph for probabilies:
# 
# ## Make graph:
# p1 <- ggplot(BACI_sl, aes(species, 0, colour = treatment))
# p2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_sl$ecdf),
#                      position = position_dodge(0.6),
#                      size = 4,
#                      width = 0)
# p2b <- geom_errorbar(aes(ymin = BACI_sl$ecdf - 1, ymax = 0),
#                      position = position_dodge(0.6),
#                      size = 4,
#                      width = 0)
# P <- p1 +
#   geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
#   scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
#                      labels = c("1", ".5", "0", "", ""),
#                      sec.axis = dup_axis(
#                        name = "Probability that the indicator is positive",
#                        labels = c("", "", "0", ".5", "1"))) +
#   p2a + p2b + g4 +
#   xlab("") + ylab("Probability that the indicator is negative") +
#   coord_flip() +
#   scale_colour_manual(values = c("#00AFBB", "#FC4E07", "#E7B800")) +
#   theme_light(70) +
#   theme(legend.position = "top",
#         legend.title = element_blank(),
#         legend.key.size = unit(5, 'lines'),
#         legend.box = "vertical",
#         legend.spacing.y = unit(0, "lines"),
#         strip.text.y = element_blank(),
#         strip.background = element_rect(colour="white", size = 0.8))
# 
# png("figures/BACI_sl_probs.png", 30000/8, 27500/8, "px", res = 600/8)
# P
# dev.off()

## 4. Make graphs for guilds & trends ------------------------------------------

## Graph for probabilies:

q1 <- ggplot(data = BACI_gl, aes(x = guild, y = 0, colour = treatment))
q2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_gl$ecdf), 
                     position = position_dodge(0.5),
                     size = 5, 
                     width = 0)
q2b <- geom_errorbar(aes(ymin = BACI_gl$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 5,                  
                     width = 0)
q3 <- facet_nested(vars(group), vars(ref, indicator), "free", scales = "free")
Q <- q1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  q2a + q2b + q3 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c("#00AFBB", "#FC4E07", "#E7B800")) +
  theme_light(70) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(5, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.background = element_rect(colour = "white", size = 0.8))

png("figures/BACI_gl_probs.png", 30000/8, 20000/8, "px", res = 600/8)
Q
dev.off()

## -------------------------------END-------------------------------------------
