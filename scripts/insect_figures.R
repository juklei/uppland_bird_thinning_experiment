## Make figures for the BACI study results:
##
## First edit: 20200302
## Last edit:  20200302
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("rjags")
require("data.table")
require("wesanderson")
require("ggpubr")
require("cowplot")

## 2. Load and transform data --------------------------------------------------

cover_control <- read.csv("clean/BACI_cover_control.csv")
cover_CR <- read.csv("clean/BACI_cover_CR.csv")
pm_control <- read.csv("clean/BACI_pm_control.csv")
pm_CR <- read.csv("clean/BACI_pm_CR.csv")

## Adjust reference category here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ref <- "control"
ref <- "CR"
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(ref == "control"){
  BACI_i <- rbind(cover_control, pm_control) 
}  else{ 
  BACI_i <- rbind(cover_CR, pm_CR) }

BACI_i$response[1:9] <- "Mean increment (% trap cover) per day"
BACI_i$response[10:18] <- "Day after march with highest increment"

levels(BACI_i$treatment) <- c("Complete retention", 
                              "Conventional thinning", 
                              "Understory retention thinning")
BACI_i$treatment <- factor(BACI_i$treatment, 
                           levels = c("Complete retention", 
                                      "Understory retention thinning",
                                      "Conventional thinning"))
levels(BACI_i$indicator) <- c("BACI-contrast", 
                              "CI-contribution", 
                              "CI-divergence")

## 3. Make graphs --------------------------------------------------------------

## Graph for slopes with CIs:

g1 <- ggplot(data = BACI_i[1:9, ], 
             aes(x = response, y = X50., colour = treatment))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 3.5, 
                    width = 0, 
                    position = position_dodge(0.5))
g3 <- geom_point(aes(fill = treatment), 
                 position = position_dodge(0.5), 
                 size = 4.5, 
                 colour = "black")
g4 <- facet_grid(. ~ indicator, scales = "free")

G <- g1 + geom_hline(yintercept = 0, size = 1, color = "darkgrey") + 
  g2 + g3 + g4 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(30) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        #legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

p1 <- ggplot(data = BACI_i[10:18, ],
             aes(x = response, y = X50., colour = treatment))

P <- p1 + geom_hline(yintercept = 0, size = 1, color = "darkgrey") + 
  g2 + g3 + g4 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(35) +
  theme(legend.position = "none")  

png(paste0("figures/BACI_i_", ref, "_slopes.png"), 
    12000/8, 5000/8, 
    "px", 
    res = 600/8)
plot_grid(G, P, align = "v", nrow = 2, rel_heights = c(4/7, 3/7))
dev.off()

## Graph for probs:

p1 <- ggplot(data = BACI_i, aes(x = response, y = 0, colour = treatment))
p2a <- geom_errorbar(aes(ymin = 0, ymax = BACI_i$ecdf), 
                     position = position_dodge(0.5),
                     size = 4, 
                     width = 0)
p2b <- geom_errorbar(aes(ymin = BACI_i$ecdf - 1, ymax = 0), 
                     position = position_dodge(0.5),
                     size = 4,                  
                     width = 0)
p3 <- facet_grid(. ~ indicator)

P <- p1 + geom_hline(yintercept = 0, size = 1, color = "darkgrey") + 
  p2a + p2b + p3 +
  xlab("") + ylab("") + coord_flip() +
  scale_colour_manual(values = c(ifelse(ref == "control", "#00AFBB", "#E7B800"), 
                                 "#FC4E07")) +
  theme_light(35) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(3, 'lines'),
        #legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank())

png(paste0("figures/BACI_i_", ref, "_probs.png"), 
    12000/8, 3750/8, 
    "px", 
    res = 600/8)
P
dev.off()

## -------------------------------END-------------------------------------------
