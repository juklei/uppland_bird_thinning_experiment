## Make figures for the experimental data and the bpo and nestbox model

# To do: Add number of pbservations per species to the species figure. 
# State why the numbers: Caution with interpreting effects.
# Make a figure like the guilds figure with the nestboxes with response variable
# instead of guilds.

## First edit: 20191028
## Last edit:  20200515
##
## Author: Julian Klein

## 1. Clear environment and load libraries -------------------------------------

# rm(list = ls())

require("ggplot2")
require("data.table")
# devtools::install_github("teunbrand/ggh4x")
require("ggh4x")

## 2. Load and prepare data ----------------------------------------------------

bird_data <- read.csv("data/bird_data.csv")
bpo <- read.csv("clean/bpo_double.csv")
BACI_sl_NF <- read.csv("clean/BACI_sl_red_ref_NF.csv")
BACI_sl_CR <- read.csv("clean/BACI_sl_red_ref_CR.csv")
BACI_gl_NF <- read.csv("clean/BACI_gl_red_ref_NF.csv")
BACI_gl_CR <- read.csv("clean/BACI_gl_red_ref_CR.csv")
BACI_nb_NF <- read.csv("clean/BACI_nb_ref_NF.csv")
BACI_nb_CR <- read.csv("clean/BACI_nb_ref_CR.csv")
BACI_rs_NF <- read.csv("clean/BACI_rs_ref_NF.csv")
BACI_rs_CR <- read.csv("clean/BACI_rs_ref_CR.csv")

## Exclude URT from NF controls:
BACI_sl_NF <- droplevels(BACI_sl_NF[BACI_sl_NF$treatment != "URT", ])
BACI_gl_NF <- droplevels(BACI_gl_NF[BACI_gl_NF$treatment != "URT", ])
BACI_nb_NF <- droplevels(BACI_nb_NF[BACI_nb_NF$treatment != "URT", ])
BACI_rs_NF <- droplevels(BACI_rs_NF[BACI_rs_NF$treatment != "URT", ])

## Combine all data:
BACI_sl <- rbind(BACI_sl_NF, BACI_sl_CR)
BACI_gl <- rbind(BACI_gl_NF, BACI_gl_CR)
BACI_nb <- rbind(BACI_nb_NF, BACI_nb_CR)
BACI_rs <- rbind(BACI_rs_NF, BACI_rs_CR)

## Change level names:
l1 <- c("Complete retention", 
        "Conventional thinning", 
        "Understory retention thinning")
levels(BACI_sl$treatment) <-  l1
levels(BACI_gl$treatment) <-  l1
levels(BACI_nb$treatment) <-  l1
levels(BACI_rs$treatment) <-  l1
l2 <- c("BACI-contrast", "CI-contribution", "CI-divergence")
levels(BACI_sl$indicator) <- l2
levels(BACI_gl$indicator) <- l2
levels(BACI_nb$indicator) <- l2
levels(BACI_rs$indicator) <- l2
l3 <- c("Control = No forestry", "Control = Complete retention")
levels(BACI_sl$ref) <- l3
levels(BACI_gl$ref) <- l3
levels(BACI_nb$ref) <- l3
levels(BACI_rs$ref) <- l3

## Conventional should have same colour in both panels:
l4 <- c("Understory retention thinning",
        "Complete retention", 
        "Conventional thinning")
BACI_sl$treatment <- factor(BACI_sl$treatment, l4) 
BACI_gl$treatment <- factor(BACI_gl$treatment, l4) 
BACI_nb$treatment <- factor(BACI_nb$treatment, l4) 
BACI_rs$treatment <- factor(BACI_rs$treatment, l4) 

## Add raw observation numbers to name:
bpo <- as.data.table(bpo)
nobs <- bpo[, list(sum_nobs = sum(n_obs)), by = "species"]
nobs <- merge(nobs, bird_data, by.x = "species", by.y = "short")
BACI_sl <- merge(BACI_sl, nobs[, c("long", "sum_nobs")], by.x = "species", by.y = "long")
BACI_sl$species_nobs <- as.factor(paste0(BACI_sl$species, " (", BACI_sl$sum_nobs, ")"))

## 3. Make combined figure for BACI_sl -----------------------------------------

## Graph with slopes:

## categorise responses:
BACI_sl$cat <- ifelse(BACI_sl$species == "Community mean", "cm", "species")

## Start species names with "A" on top:
BACI_sl$species <- factor(BACI_sl$species, rev(levels(BACI_sl$species)))

## Make figure:
g1 <- ggplot(BACI_sl, aes(species, X50., colour = treatment, fill = treatment))
g2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 5, 
                    width = 0, 
                    # alpha = 0.8,
                    position = position_dodge(0.6))
g3 <- geom_point(position = position_dodge(0.6), size = 7, colour = "black")
g4 <- facet_nested(vars(cat), vars(ref, indicator), "free")
G <- g1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") + 
  g2 + g4 + g3 +
  xlab("") + ylab("Indicator value for probability of occurrence") + 
  coord_flip() +
  scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07")) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1", "-.5", "0", ".5", "1")) +
  theme_light(85) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(7, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        legend.text = element_text(size = 95),
        axis.title.x = element_text(size = 100),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour = "white", size = 1.5))

png("figures/BACI_sl_red_slopes_new.png", 40000/8, 30000/8, "px", res = 600/8)
G
dev.off()

## Graph for probabilies:

## Start species names with "A" on top:
BACI_sl$species_nobs <- factor(BACI_sl$species_nobs, rev(levels(BACI_sl$species_nobs)))

## Make graph:
p1 <- ggplot(BACI_sl, aes(species_nobs, 0, colour = treatment))
p2 <- geom_errorbar(aes(ymin = ifelse(ecdf < 0.5, ecdf - 1, 0), 
                         ymax = ifelse(ecdf < 0.5, 0, ecdf)),
                     position = position_dodge(0.6),
                     size = 5,
                     width = 0)
P <- p1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  p2 + g4 +
  xlab("") + ylab("Probability that the indicator is negative") +
  coord_flip() +
  scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07")) +
  theme_light(70) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(8, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="white", size = 1.5))

png("figures/BACI_sl_probs2.png", 30000/8, 28500/8, "px", res = 600/8)
P
dev.off()

## 4. Make graphs for guilds & trends ------------------------------------------

## Graph for slopes:

h1 <- ggplot(BACI_gl, aes(guild, X50., colour = treatment, fill = treatment))
h2 <- geom_errorbar(aes(ymin = X2.5., ymax = X97.5.), 
                    size = 5, 
                    width = 0, 
                    position = position_dodge(0.6))
h3 <- geom_point(position = position_dodge(0.6), size = 7, colour = "black")
h4 <- facet_nested(vars(group), vars(ref, indicator), "free", scales = "free_y")
H <- h1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  h2 + h3 + h4 +
  xlab("") + ylab("Indicator value for probability of occurrence") + 
  coord_flip() +
  scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07")) +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1", "-.5", "0", ".5", "1")) +
  theme_light(70) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(5, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.background = element_rect(colour = "white", size = 1.5))

png("figures/BACI_gl_red_slopes_new.png", 31000/8, 19000/8, "px", res = 600/8)
H
dev.off()

## Graph for probabilies:

q1 <- ggplot(data = BACI_gl, aes(x = guild, y = 0, colour = treatment))
q2a <- geom_errorbar(aes(ymin = ifelse(ecdf < 0.5, ecdf - 1, 0), 
                         ymax = ifelse(ecdf < 0.5, 0, ecdf)),
                     position = position_dodge(0.5),
                     size = 5,
                     width = 0)
# q2b <- geom_errorbar(aes(ymin = BACI_gl$ecdf - 1, ymax = 0), 
#                      position = position_dodge(0.5),
#                      size = 5,                  
#                      width = 0)
q3 <- facet_nested(vars(group), vars(ref, indicator), "free", scales = "free_y")
Q <- q1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  q2a + q3 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07")) +
  theme_light(70) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(5, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.background = element_rect(colour = "white", size = 1.5))

png("figures/BACI_gl_red_probs_2.png", 30000/8, 20000/8, "px", res = 600/8)
Q
dev.off()

## 5. Make graph for nestbox data ----------------------------------------------

## Combine nb & rs:
levels(BACI_nb$species)[c(1,3,4)] <- c("Cyanistes careuleus", 
                                       "Ficedula hypoleuca", 
                                       "Parus major")
BACI_nb$response <- "box occupancy"
BACI_rs$species[BACI_rs$response == "lambda_post"] <- "Parus major (no. fledglings)"
BACI_rs$species[BACI_rs$response == "p_post"] <- "Parus major (nestfailure)"
BACI_rs$response <- "reproduction"
BACI_nestbox <- rbind(BACI_nb, BACI_rs[, c(ncol(BACI_rs), 2:ncol(BACI_rs)-1)])

## Graph for probabilies:

## Don't display empty nests:
BACI_nestbox <- droplevels(BACI_nestbox[BACI_nestbox$species != "empty", ])

s1 <- ggplot(BACI_nestbox, aes(x = species, y = 0, colour = treatment))
s2a <- geom_errorbar(aes(ymin = ifelse(ecdf < 0.5, ecdf - 1, 0), 
                         ymax = ifelse(ecdf < 0.5, 0, ecdf)),
                     position = position_dodge(0.5),
                     size = 8,
                     width = 0)
# s2b <- geom_errorbar(aes(ymin = ecdf - 1, ymax = 0), 
#                      position = position_dodge(0.5),
#                      size = 8,                  
#                      width = 0)
s3 <- facet_nested(vars(response), vars(ref, indicator), "free", scales = "free_y")
S <- s1 +
  geom_hline(yintercept = 0, size = 2, color = "darkgrey") +
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("1", ".5", "0", "", ""),
                     sec.axis = dup_axis(
                       name = "Probability that the indicator is positive",
                       labels = c("", "", "0", ".5", "1"))) +
  s2a + s3 +
  xlab("") + ylab("Probability that the indicator is negative") + 
  coord_flip() +
  scale_colour_manual(values = c("#E7B800", "#00AFBB", "#FC4E07")) +
  theme_light(70) +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.key.size = unit(5, 'lines'),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "lines"),
        strip.background = element_rect(colour = "white", size = 1.5))

png("figures/BACI_nb_probs2.png", 30000/8, 14000/8, "px", res = 600/8)
S
dev.off()

## -------------------------------END-------------------------------------------
