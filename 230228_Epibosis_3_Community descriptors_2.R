## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between euendolithic infestation
##              and epibiosis in the brown mussel Perna perna
##              Script 3 - Community descriptors 2 (species richness, diversity and evenness)
##
##
## Purpose of script: 
##
##
## Author: Alexia Dievart
##
## Date Created: 2023-02-28
## Dates Updated:
##
## Copyright (c) Alexia DIEVART 2023
## Email: alexia.dievart@hotmail.fr
##
## ---------------------------
##
## Notes: 
##   - 1 mussel = 2 valves
## ---------------------------

###############################################################################################
# Section: Session setup ----------------------------------------------------------------------
###############################################################################################

# Set working directory
setwd("D:/Little Bull/Etudes/phD - Rhodes/1.3_EPIBIOSIS/STATS")


# Load packages
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(tidyverse,
               ggplot2,
               ggpubr,
               rstatix,
               broom,
               dplyr,
               emmeans,
               DHARMa,
               multcomp)

# Set default ggplot theme
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  legend.position = "none"))

# Create colour vectors
col_inf <- c("coral4", "brown3", "indianred1", "pink3", "darkgrey") # 5 levels: A, B, C, D and E
col_deg <- c("chocolate4", "chocolate", "darkgoldenrod3", "darkgoldenrod1", "bisque2", "cornsilk2") # 6 levels
col_epibio <- c("darkseagreen", "cadetblue3")

colors_position <- c("tomato4", "cadetblue", "darkgrey", "brown3")
colors_position1 <- c("darkgrey", "brown3")
colors_position2 <- c("cadetblue", "darkgrey", "brown3")





###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Load and curate full dataset (long format) ==================================================
epibio <- read.csv("epibiosis.csv", h=T, dec=",", sep=";")
dplyr::glimpse(epibio)
View(epibio)

# Transform character variables into factor variables
epibio_factor <- c("location", "quadrat", "specimen", "mussel.species", "sex", "valve", "inf.lvl", "inf.percent",
                   "higher.group", "epibiont", "position", "basibiont")
epibio[,epibio_factor] <- lapply(epibio[,epibio_factor], factor)

# Recode "inf.percent" into a shell degradation index (from 0 to 6) as in Ma et al. 2023
epibio$inf.percent <- recode_factor(epibio$inf.percent, "0" = "1", "0-25" = "2", "25-50" = "3", "50-75" = "4",
                                    "75-100" = "5", "100" = "6")

# Exclude all Mytilus galloprovincialis mussel shells
epibio_perna <- epibio[! epibio$mussel.species == "Mytilus galloprovincialis",]

epibio_factor <- c("location", "quadrat", "specimen", "mussel.species", "sex", "valve", "inf.lvl", "inf.percent",
                   "higher.group", "epibiont", "position", "basibiont")
epibio_perna[,epibio_factor] <- lapply(epibio_perna[,epibio_factor], factor)

dplyr::glimpse(epibio_perna)
View(epibio_perna)

# Excluse mussel valves without epibionts
epibio_perna1 <- epibio_perna[! epibio_perna$abundance == 0,]
View(epibio_perna1)

# Isolate ABUNDANCE ===========================================================================
abund_long <- epibio_perna1[,1:16]
abund_long <- abund_long[! abund_long$position == "byssus",]
View(abund_long)

# Create shell length categories ####
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
abund_long$length_range <- cut(abund_long$shell.length, breaks = length_range_bins,
                               labels = length_range_labels)



###############################################################################################
# Section: Curate and explore data ------------------------------------------------------------
###############################################################################################

###############################################################################################
# Unique species and common species ===========================================================
###############################################################################################

# Total number of species ####
abund_long %>% summarise(n_distinct(epibiont))

# Total number of sessile species
abund_long1 <- abund_long[abund_long$mobility == "sessile",]
View(abund_long1)

abund_long1 %>% summarise(n_distinct(epibiont))

# Total number of semi-sessile species
abund_long2 <- abund_long[abund_long$mobility == "semi-sessile",]
View(abund_long2)

abund_long2 %>% summarise(n_distinct(epibiont))

# Total number of mobile species
abund_long3 <- abund_long[abund_long$mobility == "mobile",]
View(abund_long3)

abund_long3 %>% summarise(n_distinct(epibiont))

# Number of unique species ####

# Number of unique species per shell length range
abund_long %>%
  group_by(length_range) %>%
  summarise(n_distinct(epibiont))

# Number of unique species per infestation
abund_long %>%
  group_by(inf.lvl) %>%
  summarise(n_distinct(epibiont))

# Number of unique species per degradation
abund_long %>%
  group_by(inf.percent) %>%
  summarise(n_distinct(epibiont))

# Number of unique species per position
abund_long %>%
  group_by(position) %>%
  summarise(n_distinct(epibiont))

# Number of common species ####

# Number and ID of common species across infestation
Reduce(intersect, split(abund_long$epibiont, abund_long$inf.lvl))

# Number and ID of common species across degradation
Reduce(intersect, split(abund_long$epibiont, abund_long$inf.percent))

# Number and ID of common species across position
Reduce(intersect, split(abund_long$epibiont, abund_long$position))





###############################################################################################
# Section 1: Species richness (S) ------------------------------------------------------------
###############################################################################################

# 1.1. Prepare data set =======================================================================
View(abund_long)

# Summarize species richness per length range ####
rich_length <- abund_long %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  summarise(count = n_distinct(epibiont))
View(rich_length)

# Summarize species richness per infestation and degradation ####
rich_inf <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  summarise(count = n_distinct(epibiont))
View(rich_inf)

# Summarize species richness per position ####
rich_pos <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, position) %>%
  summarise(count = n_distinct(epibiont))
View(rich_pos)





# 1.2. Visualizing data set ==================================================================

# BARPLOT - Species richness x Length range ####

# For each quadrat
rich_length1 <- rich_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_length1)

barplot_rich_length1 <- ggplot(rich_length1, aes(x = length_range, y = count_mean,
                                              col = length_range, fill = length_range)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Species richness")
barplot_rich_length1

# Across all quadrats
rich_length2 <- rich_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_length2)

barplot_rich_length2 <- ggplot(rich_length2, aes(x = length_range, y = count_mean, 
                                                          color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(2.5, 3, 4, 5.5, 6, 6.5, 6.5, 8, 6.5),
           label = c("A", "A", "A", "B", "BC", "CD", "CD", "D", "BCD"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean species richness (S)")
barplot_rich_length2

# BARPLOT - Species richness x Infestation ####

# For each quadrat
rich_inf1 <- rich_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_inf1)

barplot_rich_inf1 <- ggplot(rich_inf1, aes(x = inf.lvl, y = count_mean,
                                                 col = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Species richness")
barplot_rich_inf1

# Across all quadrats
rich_inf2 <- rich_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_inf2)

barplot_rich_inf2 <- ggplot(rich_inf2, aes(x = inf.lvl, y = count_mean,
                                           col = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(2.3, 2.5, 3, 3, 3.2),
           label = c("AB", "B", "A", "A", "AB"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean species richness per valve (S)")
barplot_rich_inf2

# BARPLOT - Species richness x Degradation ####

# For each quadrat
rich_deg1 <- rich_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_deg1)

barplot_rich_deg1 <- ggplot(rich_deg1, aes(x = inf.percent, y = count_mean,
                                           col = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Species richness")
barplot_rich_deg1

# Across all quadrats
rich_deg2 <- rich_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_deg2)

barplot_rich_deg2 <- ggplot(rich_deg2, aes(x = inf.percent, y = count_mean,
                                           col = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean species richness per valve (S)")
barplot_rich_deg2

# BARPLOT - Species richness x Position ####

# For each quadrat
rich_pos1 <- rich_pos %>%
  dplyr::group_by(quadrat, position) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_pos1)

barplot_rich_pos1 <- ggplot(rich_pos1, aes(x = position, y = count_mean,
                                           col = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Species richness")
barplot_rich_pos1

# Across all quadrats
rich_pos2 <- rich_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(rich_pos2)

barplot_rich_pos2 <- ggplot(rich_pos2, aes(x = position, y = count_mean,
                                           col = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3),
           y = c(2.8, 3.2, 2.8),
           label = c("A", "B", "A"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean species richness per valve (S)")
barplot_rich_pos2


# 1.3. MODELS: Species richness vs Length range ============================

# Investigate data set ####
View(rich_length)
hist(rich_length$count) # Poisson distribution
# Species richness is considered for each mussel - because each valve has the same length

# Fit model: LM ####
rich_length_lm <- lm(count ~ length_range, data = rich_length)
summary(rich_length_lm)
Anova(rich_length_lm)
AIC(rich_length_lm) # AIC(lm) = 1802.372

# Check model: LM
par(mfrow=c(2,2))
plot(rich_length_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = rich_length_lm, plot = T) # Not a good fit.

# Fit model: GLM (poisson) w/ length range ####
rich_length_glm <- glm(count ~ length_range, data = rich_length,
                        family = "poisson")
summary(rich_length_glm)
Anova(rich_length_glm)
AIC(rich_length_glm) # AIC(glm) = 1715.783

# Check models
DHARMa::simulateResiduals(fittedModel = rich_length_glm, plot = T) # Much better but variance not homogeneous
gratia::appraise(rich_length_glm, method = "simulate")

# Fit model : GLM (negative binomial) w/ length range ####
rich_length_glm2 <- glm.nb(count ~ length_range, data = rich_length)
summary(rich_length_glm2)
Anova(rich_length_glm2)
AIC(rich_length_glm2) # AIC(glm) = 1717.787

# Check models
DHARMa::simulateResiduals(fittedModel = rich_length_glm2, plot = T) # Better than Poisson model.
gratia::appraise(rich_length_glm2, method = "simulate")

# Pairwise comparisons
summary(glht(rich_length_glm, linfct = mcp(length_range = "Tukey")))





# 1.4. MODELS: Species richness vs Infestation x Degradation =================================

# Investigate data set ####
View(rich_inf)
hist(rich_inf$count) # Poisson distribution
# Species richness is considered for each mussel valve - because each valve has different parameters

# Fit model: LM saturated ####
rich_inf_lm <- lm(count ~ inf.lvl * inf.percent, data = rich_inf)
summary(rich_inf_lm)
Anova(rich_inf_lm)
AIC(rich_inf_lm) # AIC(lm) = 2826.975

# Check model: LM saturated
par(mfrow=c(2,2))
plot(rich_inf_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = rich_inf_lm, plot = T) # Not a good fit at all.

# Fit model: LM additive ####
rich_inf_lm1 <- lm(count ~ inf.lvl + inf.percent, data = rich_inf)
summary(rich_inf_lm1)
Anova(rich_inf_lm1)
AIC(rich_inf_lm1) # AIC(lm) = 2826.14

# Check model: LM additive
par(mfrow=c(2,2))
plot(rich_inf_lm1) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = rich_inf_lm1, plot = T) # Not a good fit at all.

# Fit model: GLM (poisson) saturated ####
rich_inf_glm <- glm(count ~ inf.lvl * inf.percent, data = rich_inf,
                     family = "poisson")
summary(rich_inf_glm)
Anova(rich_inf_glm)
AIC(rich_inf_glm) # AIC(glm) = 2697.858

# Check models
DHARMa::simulateResiduals(fittedModel = rich_inf_glm, plot = T) # Better fit than LM

# Fit model: GLM (poisson) additive ####
rich_inf_glm1 <- glm(count ~ inf.lvl + inf.percent, data = rich_inf,
                      family = "poisson")
summary(rich_inf_glm1)
Anova(rich_inf_glm1)
AIC(rich_inf_glm1) # AIC(glm) = 2688.802

# Check models
DHARMa::simulateResiduals(fittedModel = rich_inf_glm1, plot = T) # Not a good fit at all.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = rich_inf_glm1, plot = T)
plotResiduals(simulationOutput, rich_inf$inf.lvl)  # within-group deviations significant for A and B
plotResiduals(simulationOutput, rich_inf$inf.percent) # within-group deviations for 1 to 3

# Fit model : GLM (negative binomial) saturated ####
rich_inf_glm2 <- glm.nb(count ~ inf.lvl * inf.percent, data = rich_inf)
summary(rich_inf_glm2)
Anova(rich_inf_glm2)
AIC(rich_inf_glm2) # AIC(glm) = 2699.868

# Check models
DHARMa::simulateResiduals(fittedModel = rich_inf_glm2, plot = T) # Better fit.
gratia::appraise(rich_inf_glm2, method = "simulate")
simulationOutput <- DHARMa::simulateResiduals(fittedModel = rich_inf_glm2, plot = T)
plotResiduals(simulationOutput, rich_inf$inf.lvl)  # within-group deviations significant for A and B
plotResiduals(simulationOutput, rich_inf$inf.percent) # within-group deviations for 1 and 3

# Fit model : GLM (negative binomial) saturated ####
rich_inf_glm3 <- glm.nb(count ~ inf.lvl + inf.percent, data = rich_inf)
summary(rich_inf_glm3)
Anova(rich_inf_glm3)
AIC(rich_inf_glm3) # AIC(glm) = 2690.811

# Check models
DHARMa::simulateResiduals(fittedModel = rich_inf_glm3, plot = T) # Better fit.
gratia::appraise(rich_inf_glm3, method = "simulate")
simulationOutput <- DHARMa::simulateResiduals(fittedModel = rich_inf_glm3, plot = T)
plotResiduals(simulationOutput, rich_inf$inf.lvl)  # within-group deviations significant for A and B
plotResiduals(simulationOutput, rich_inf$inf.percent) # within-group deviations for 1

# Pairwise comparisons - additive negative binomial
summary(glht(rich_inf_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(rich_inf_glm3, linfct = mcp(inf.percent = "Tukey")))





# 1.3. MODELS: Species richness vs Position ============================

# Investigate data set ####
View(rich_pos)
hist(rich_pos$count) # Poisson distribution
# Species richness is considered for each mussel - because each valve has the same length

# Fit model: LM ####
rich_pos_lm <- lm(count ~ position, data = rich_pos)
summary(rich_pos_lm)
Anova(rich_pos_lm)
AIC(rich_pos_lm) # AIC(lm) = 3563.373

# Check model: LM
par(mfrow=c(2,2))
plot(rich_pos_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = rich_pos_lm, plot = T) # Not a good fit.

# Fit model: GLM (poisson) saturated ####
rich_pos_glm <- glm(count ~ position, data = rich_pos,
                    family = "poisson")
summary(rich_pos_glm)
Anova(rich_pos_glm)
AIC(rich_pos_glm) # AIC(glm) = 3649.136

# Check models
DHARMa::simulateResiduals(fittedModel = rich_pos_glm, plot = T) # Better fit than LM

# Fit model : GLM (negative binomial) saturated ####
rich_pos_glm1 <- glm.nb(count ~ position, data = rich_pos)
summary(rich_pos_glm1)
Anova(rich_pos_glm1)
AIC(rich_pos_glm1) # AIC(glm) = 3651.155

# Check models
DHARMa::simulateResiduals(fittedModel = rich_pos_glm1, plot = T) # Better fit.

# Pairwise comparisons
summary(glht(rich_pos_glm, linfct = mcp(position = "Tukey")))





###############################################################################################
# Section 2: Shannon's Diversity (H') ---------------------------------------------------------
###############################################################################################

# 2.1. Prepare data set =======================================================================
View(abund_long)

# Summarize diversity per shell length  ####
abund_long_length <- abund_long %>%
  dplyr::group_by(location, quadrat, specimen, sex, epibiont, length_range) %>%
  summarise(abundance = sum(abundance))
View(abund_long_length)

abund_wide_length <- spread(abund_long_length, epibiont, abundance)
abund_wide_length[is.na(abund_wide_length)] <- 0
View(abund_wide_length)

shannon_length <- diversity(abund_wide_length[,6:60], "shannon")
View(shannon_length)

abund_div_length <- abund_wide_length[, 1:5] %>%
  add_column(shannon_length, .after = "length_range")
View(abund_div_length)

# Summarize diversity per infestation and degradation ####
abund_long_inf <- abund_long %>%
  dplyr::group_by(location, quadrat, specimen, sex, valve, inf.lvl, inf.percent, epibiont) %>%
  summarise(abundance = sum(abundance))
View(abund_long_inf)

abund_wide_inf <- spread(abund_long_inf, epibiont, abundance)
abund_wide_inf[is.na(abund_wide_inf)] <- 0
View(abund_wide_inf)

shannon_inf <- diversity(abund_wide_inf[,8:62], "shannon")

abund_div_inf<- abund_wide_inf[, 1:7] %>%
  add_column(shannon_inf, .after = "inf.percent")
View(abund_div_inf)

# Summarize diversity per position ####
abund_long_pos <- abund_long %>%
  dplyr::group_by(location, quadrat, specimen, sex, valve, position, epibiont) %>%
  summarise(abundance = sum(abundance))
View(abund_long_pos)

abund_wide_pos <- spread(abund_long_pos, epibiont, abundance)
abund_wide_pos[is.na(abund_wide_pos)] <- 0
View(abund_wide_pos)

shannon_pos <- diversity(abund_wide_pos[,7:61], "shannon")

abund_div_pos <- abund_wide_pos[, 1:6] %>%
  add_column(shannon_pos, .after = "position")
View(abund_div_pos)






# 2.2. Visualizing data set ==================================================================

# BARPLOT - Shannon's diversity x Length range ####

# For each quadrat
shannon_length1 <- abund_div_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  summarise(
    shannon_mean = mean(shannon_length),
    shannon_sd = sd(shannon_length)
  )
View(shannon_length1)

barplot_shannon_length1 <- ggplot(shannon_length1, aes(x = length_range, y = shannon_mean,
                                                 col = length_range, fill = length_range)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Shannon-Wiener's Diversity (H')")
barplot_shannon_length1

# Across all quadrats
shannon_length2 <- abund_div_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    shannon_mean = mean(shannon_length),
    shannon_sd = sd(shannon_length)
  )
View(shannon_length2)

barplot_shannon_length2 <- ggplot(shannon_length2, aes(x = length_range, y = shannon_mean, 
                                                 color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(0.7, 0.8, 1.1, 1.4, 1.6, 1.5, 1.5, 1.7, 1.3),
           label = c("A", "A", "A", "B", "BC", "C", "BC", "C", "ABC"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Shannon-Wiener's Diversity (H')")
barplot_shannon_length2

# BARPLOT - Shannon's diversity x Infestation ####

# For each quadrat
shannon_inf1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_inf1)

barplot_shannon_inf1 <- ggplot(shannon_inf1, aes(x = inf.lvl, y = shannon_mean,
                                                       col = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Shannon-Wiener's Diversity (H')")
barplot_shannon_inf1

# Across all quadrats
shannon_inf2 <- abund_div_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_inf2)

barplot_shannon_inf2 <- ggplot(shannon_inf2, aes(x = inf.lvl, y = shannon_mean, 
                                                       color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(0.7, 0.9, 1.2, 1.3, 1.2),
           label = c("AB", "A", "B", "C", "ABC"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean Shannon-Wiener's Diversity (H')")
barplot_shannon_inf2

# BARPLOT - Shannon's diversity x Degradation ####

# For each quadrat
shannon_deg1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_deg1)

barplot_shannon_deg1 <- ggplot(shannon_deg1, aes(x = inf.percent, y = shannon_mean,
                                                 col = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Shannon-Wiener's Diversity (H')")
barplot_shannon_deg1

# Across all quadrats
shannon_deg2 <- abund_div_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_deg2)

barplot_shannon_deg2 <- ggplot(shannon_deg2, aes(x = inf.percent, y = shannon_mean, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean Shannon-Wiener's Diversity (H')")
barplot_shannon_deg2

# BARPLOT - Shannon's diversity x Position ####

# For each quadrat
shannon_pos1 <- abund_div_pos %>%
  dplyr::group_by(quadrat, position) %>%
  summarise(
    shannon_mean = mean(shannon_pos),
    shannon_sd = sd(shannon_pos)
  )
View(shannon_pos1)

barplot_shannon_pos1 <- ggplot(shannon_pos1, aes(x = position, y = shannon_mean,
                                                 col = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Shannon-Wiener's Diversity (H')")
barplot_shannon_pos1

# Across all quadrats
shannon_pos2 <- abund_div_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    shannon_mean = mean(shannon_pos),
    shannon_sd = sd(shannon_pos)
  )
View(shannon_pos2)

barplot_shannon_pos2 <- ggplot(shannon_pos2, aes(x = position, y = shannon_mean, 
                                                 color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3),
           y = c(0.8, 0.9, 0.8),
           label = c("AB", "A", "B"),
           size = 4) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean Shannon-Wiener's Diversity (H')")
barplot_shannon_pos2





# 2.3. MODELS: Shannon Diversity vs Length range ============================

# Investigate data set ####
View(abund_div_length)
hist(abund_div_length$shannon_length) # Normalish distribution
# Species richness is considered for each mussel - because each valve has the same length

# Fit model: LM ####
shannon_length_lm <- lm(shannon_length ~ length_range, data = abund_div_length)
summary(shannon_length_lm)
Anova(shannon_length_lm)
AIC(shannon_length_lm) # AIC(lm) = 313.831

# Check model: LM
par(mfrow=c(2,2))
plot(shannon_length_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = shannon_length_lm, plot = T) # Not a bad fit at all !

# Pairwise comparisons ####
summary(glht(shannon_length_lm, linfct = mcp(length_range = "Tukey")))

# 2.4. MODELS: Shannon Diversity vs (Infestation x Degradation) =========

# Investigate data set ####
View(abund_div_inf)
hist(abund_div_inf$shannon_inf) # Normalish distribution
# Species richness is considered for each valve

# Fit model: LM saturated ####
shannon_inf_lm <- lm(shannon_inf ~ inf.lvl * inf.percent, data = abund_div_inf)
summary(shannon_inf_lm)
Anova(shannon_inf_lm)
AIC(shannon_inf_lm) # AIC(lm) = 1060.22

# Check model: LM saturated
par(mfrow=c(2,2))
plot(shannon_inf_lm) # Structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = shannon_inf_lm, plot = T)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = shannon_inf_lm, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 3, and 5

# Fit model: LM additive ####
shannon_inf_lm1 <- lm(shannon_inf ~ inf.lvl + inf.percent, data = abund_div_inf)
summary(shannon_inf_lm1)
Anova(shannon_inf_lm1)
AIC(shannon_inf_lm1) # AIC(lm) = 1056.43

# Check model: LM 
par(mfrow=c(2,2))
plot(shannon_inf_lm1) # Structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = shannon_inf_lm1, plot = T)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = shannon_inf_lm1, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 3, and 5

# Fit model: GLM (ziGaussian) ####
library(glmmTMB)

shannon_inf_tmb <- glmmTMB(shannon_inf ~ inf.lvl + inf.percent,
                            ziformula = ~ inf.lvl,
                            family = gaussian(link = identity),
                            data = abund_div_inf)
summary(shannon_inf_tmb)
Anova(shannon_inf_tmb)
AIC(shannon_inf_tmb) # AIC(glmTMB) = 1066.43

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = shannon_inf_tmb, plot = T)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = shannon_inf_tmb, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 3, and 5

# Parwise comparisons w/ LM additive ####
summary(glht(shannon_inf_lm1, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(shannon_inf_lm1, linfct = mcp(inf.percent = "Tukey")))





# 2.5. MODELS: Shannon Diversity vs Position ============================

# Investigate data set ####
View(abund_div_pos)
hist(abund_div_pos$shannon_pos) # Normalish distribution

# Fit model: LM saturated ####
shannon_pos_lm <- lm(shannon_pos ~ position, data = abund_div_pos)
summary(shannon_pos_lm)
Anova(shannon_pos_lm)
AIC(shannon_pos_lm) # AIC(lm) = 1370.979

# Check model: LM saturated
par(mfrow=c(2,2))
plot(shannon_pos_lm) # Structure in the Residuals vs Leverage, distribution not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = shannon_pos_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGaussian) ####
library(glmmTMB)

shannon_pos_tmb <- glmmTMB(shannon_pos ~ position,
                           ziformula = ~ 1,
                           family = gaussian(link = "identity"),
                           data = abund_div_pos)
summary(shannon_pos_tmb)
Anova(shannon_pos_tmb)
AIC(shannon_pos_tmb) # AIC(glmTMB) = 1372.979

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = shannon_pos_tmb, plot = T)

# Fit model: GLM (ziGamma) ####
library(glmmTMB)

shannon_pos_tmb1 <- glmmTMB(shannon_pos ~ position,
                           ziformula = ~ position,
                           family = ziGamma(link = "log"),
                           data = abund_div_pos)
summary(shannon_pos_tmb1)
Anova(shannon_pos_tmb1)
AIC(shannon_pos_tmb1) # AIC(glmTMB) = 1914.741

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = shannon_pos_tmb1, plot = T) # Much better fit.

# Pairwise comparisons ####
summary(glht(shannon_pos_tmb1, linfct = mcp(position = "Tukey")))





###############################################################################################
# Section 3: Simpson Diversity (λ) ---------------------------------------------------------
###############################################################################################

# 3.1. Prepare data set =======================================================================
View(abund_long)

# Summarize diversity per shell length  ####
View(abund_wide_length)

shannon_length <- diversity(abund_wide_length[,6:60], "shannon")
View(shannon_length)

simpson_length <- diversity(abund_wide_length[,6:60], index = "simpson")

View(abund_div_length)
abund_div_length <- abund_div_length %>%
  add_column(simpson_length, .after = "shannon_length")


# Summarize diversity per infestation and degradation ####
View(abund_long_inf)

shannon_inf <- diversity(abund_wide_inf[,8:62], "shannon")

simpson_inf <- diversity(abund_wide_inf[,8:62], "simpson")

View(abund_div_inf)
abund_div_inf <- abund_div_inf %>%
  add_column(simpson_inf, .after = "shannon_inf")


# Summarize diversity per position ####
View(abund_long_pos)

shannon_pos <- diversity(abund_wide_pos[,7:61], "shannon")
simpson_pos <- diversity(abund_wide_pos[,7:61], "simpson")

View(abund_div_pos)
abund_div_pos <- abund_div_pos %>%
  add_column(simpson_pos, .after = "shannon_pos")






# 3.2. Visualizing data set ==================================================================

# BARPLOT - Simpson's diversity x Length range ####

# For each quadrat
simpson_length1 <- abund_div_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  summarise(
    simpson_mean = mean(simpson_length),
    simpson_sd = sd(simpson_length)
  )
View(simpson_length1)

barplot_simpson_length1 <- ggplot(simpson_length1, aes(x = length_range, y = simpson_mean,
                                                       col = length_range, fill = length_range)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Simpson Diversity (λ)")
barplot_simpson_length1

# Across all quadrats
simpson_length2 <- abund_div_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    simpson_mean = mean(simpson_length),
    simpson_sd = sd(simpson_length)
  )
View(simpson_length2)

barplot_simpson_length2 <- ggplot(simpson_length2, aes(x = length_range, y = simpson_mean, 
                                                       color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(0.5, 0.5, 0.7, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0),
           label = c("A", "A", "A", "B", "BC", "C", "BC", "BC", "ABC"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Simpson Diversity (λ)")
barplot_simpson_length2

# BARPLOT - Simpson's diversity x Infestation ####

# For each quadrat
simpson_inf1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_inf1)

barplot_simpson_inf1 <- ggplot(simpson_inf1, aes(x = inf.lvl, y = simpson_mean,
                                                 col = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Simpson Diversity (λ)")
barplot_simpson_inf1

# Across all quadrats
simpson_inf2 <- abund_div_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_inf2)

barplot_simpson_inf2 <- ggplot(simpson_inf2, aes(x = inf.lvl, y = simpson_mean, 
                                                 color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(0.4, 0.5, 0.7, 0.7, 0.7),
           label = c("AB", "A", "B", "C", "ABC"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean Simpson Diversity (λ)")
barplot_simpson_inf2

# BARPLOT - Simpson's diversity x Degradation ####

# For each quadrat
simpson_deg1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_deg1)

barplot_simpson_deg1 <- ggplot(simpson_deg1, aes(x = inf.percent, y = simpson_mean,
                                                 col = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Simpson Diversity (λ)")
barplot_simpson_deg1

# Across all quadrats
simpson_deg2 <- abund_div_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_deg2)

barplot_simpson_deg2 <- ggplot(simpson_deg2, aes(x = inf.percent, y = simpson_mean, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean Simpson Diversity (λ)")
barplot_simpson_deg2

# BARPLOT - Simpson's diversity x Position ####

# For each quadrat
simpson_pos1 <- abund_div_pos %>%
  dplyr::group_by(quadrat, position) %>%
  summarise(
    simpson_mean = mean(simpson_pos),
    simpson_sd = sd(simpson_pos)
  )
View(simpson_pos1)

barplot_simpson_pos1 <- ggplot(simpson_pos1, aes(x = position, y = simpson_mean,
                                                 col = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Simpson Diversity (λ)")
barplot_simpson_pos1

# Across all quadrats
simpson_pos2 <- abund_div_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    simpson_mean = mean(simpson_pos),
    simpson_sd = sd(simpson_pos)
  )
View(simpson_pos2)

barplot_simpson_pos2 <- ggplot(simpson_pos2, aes(x = position, y = simpson_mean, 
                                                 color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
#  annotate("text",
#           x = c(1, 2, 3),
#           y = c(0.8, 0.9, 0.8),
#           label = c("AB", "A", "B"),
#           size = 4) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean Simpson Diversity (λ)")
barplot_simpson_pos2


# 3.3. MODELS: Simpson Diversity vs Length range ============================

# Investigate data set ####
View(abund_div_length)
hist(abund_div_length$simpson_length) # Normalish distribution
# Simpson's diversity is considered for each mussel - because each valve has the same length

# Fit model: LM ####
simpson_length_lm <- lm(simpson_length ~ length_range, data = abund_div_length)
summary(simpson_length_lm)
Anova(simpson_length_lm)
AIC(simpson_length_lm) # AIC(lm) = -8.02964

# Check model: LM
par(mfrow=c(2,2))
plot(simpson_length_lm) # Structure in the Residuals vs Leverage, distribution is normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = simpson_length_lm, plot = T) # Not a bad fit at all !

# Fit model: GLM (ziGaussian) ####
library(glmmTMB)

simpson_length_tmb <- glmmTMB(simpson_length ~ length_range,
                           ziformula = ~ 1,
                           family = gaussian(link = "identity"),
                           data = abund_div_length)
summary(simpson_length_tmb)
Anova(simpson_length_tmb)
AIC(simpson_length_tmb) # AIC(glmTMB) = -6.02964

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = simpson_length_tmb, plot = T) # A bit better, let's go with that.

# Pairwise comparisons ####
summary(glht(simpson_length_tmb, linfct = mcp(length_range = "Tukey")))






# 3.4. MODELS: Simpson Diversity vs (Infestation x Degradation) =========

# Investigate data set ####
View(abund_div_inf)
hist(abund_div_inf$simpson_inf) # Normalish distribution
# Species richness is considered for each valve

# Fit model: LM saturated ####
simpson_inf_lm <- lm(simpson_inf ~ inf.lvl * inf.percent, data = abund_div_inf)
summary(simpson_inf_lm)
Anova(simpson_inf_lm)
AIC(simpson_inf_lm) # AIC(lm) = 108.2243

# Check model: LM saturated
par(mfrow=c(2,2))
plot(simpson_inf_lm) # Structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = simpson_inf_lm, plot = T) # Not a good fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = simpson_inf_lm, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 3, and 5

# Fit model: LM additive ####
simpson_inf_lm1 <- lm(simpson_inf ~ inf.lvl + inf.percent, data = abund_div_inf)
summary(simpson_inf_lm1)
Anova(simpson_inf_lm1)
AIC(simpson_inf_lm1) # AIC(lm) = 101.4261

# Check model: LM saturated
par(mfrow=c(2,2))
plot(simpson_inf_lm1) # Structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = simpson_inf_lm1, plot = T) # Not a good fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = simpson_inf_lm1, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 5

# Fit model: GLM (ziGaussian) ####
library(glmmTMB)

simpson_inf_tmb <- glmmTMB(simpson_inf ~ inf.lvl + inf.percent,
                           ziformula = ~ inf.lvl,
                           family = gaussian(link = "identity"),
                           data = abund_div_inf)
summary(simpson_inf_tmb)
Anova(simpson_inf_tmb)
AIC(simpson_inf_tmb) # AIC(glmTMB) = 111.4261

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = simpson_inf_tmb, plot = T)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = simpson_inf_tmb, plot = T)
plotResiduals(simulationOutput, abund_div_inf$inf.lvl)  # within-group deviations significant for A to D
plotResiduals(simulationOutput, abund_div_inf$inf.percent) # within-group deviations for 1 to 3, and 5

# Fit model: GLM (ziGamma) ####
library(glmmTMB)

simpson_inf_tmb1 <- glmmTMB(simpson_inf ~ inf.lvl + inf.percent,
                            ziformula = ~ inf.lvl,
                            family = ziGamma(link = "log"),
                            data = abund_div_inf)
summary(simpson_inf_tmb1)
Anova(simpson_inf_tmb1)
AIC(simpson_inf_tmb1) # AIC(glmTMB) = 744.9309

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = simpson_inf_tmb1, plot = T) # Much better fit.

# Pairwise comparisons ####
summary(glht(simpson_inf_tmb, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(simpson_inf_tmb, linfct = mcp(inf.percent = "Tukey")))





# 3.5. MODELS: Shannon Diversity vs Position ============================

# Investigate data set ####
View(abund_div_pos)
hist(abund_div_pos$simpson_pos) # Normalish distribution

# Fit model: LM saturated ####
simpson_pos_lm <- lm(simpson_pos ~ position, data = abund_div_pos)
summary(simpson_pos_lm)
Anova(simpson_pos_lm)
AIC(simpson_pos_lm) # AIC(lm) = 75.4147

# Check model: LM saturated
par(mfrow=c(2,2))
plot(simpson_pos_lm) # Structure in the Residuals vs Leverage, distribution not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = simpson_pos_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGaussian) ####
library(glmmTMB)

simpson_pos_tmb <- glmmTMB(simpson_pos ~ position,
                           ziformula = ~ 1,
                           family = gaussian(link = "identity"),
                           data = abund_div_pos)
summary(simpson_pos_tmb)
Anova(simpson_pos_tmb)
AIC(simpson_pos_tmb) # AIC(glmTMB) = 77.4147

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = simpson_pos_tmb, plot = T)

# Fit model: GLM (ziGamma) ####
library(glmmTMB)

simpson_pos_tmb1 <- glmmTMB(simpson_pos ~ position,
                            ziformula = ~ position,
                            family = ziGamma(link = "log"),
                            data = abund_div_pos)
summary(simpson_pos_tmb1)
Anova(simpson_pos_tmb1)
AIC(simpson_pos_tmb1) # AIC(glmTMB) = 1359.7

# Check model: GLM (ziGaussian)
DHARMa::simulateResiduals(fittedModel = simpson_pos_tmb1, plot = T) # Much better fit.

# Pairwise comparisons ####
summary(glht(simpson_pos_tmb1, linfct = mcp(position = "Tukey")))





###############################################################################################
# Section 4: Pielou's evenness (J) ---------------------------------------------------------
###############################################################################################

# 4.1. Prepare data set =======================================================================
View(abund_long)
# J = H'/ln(S)

# Summarize evenness per length range ####
rich_length <- abund_long %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  summarise(count = n_distinct(epibiont))
View(rich_length)

View(abund_div_length)

abund_div_length <- abund_div_length %>%
  add_column(richness = rich_length$count, .after = "simpson_length")

pielou_length <- abund_div_length$shannon_length/log(abund_div_length$richness)

abund_div_length <- abund_div_length %>%
  add_column(pielou_length, .after = "richness")

abund_div_length <- abund_div_length[!is.na(abund_div_length$pielou_length),]

# Summarize evenness per infestation and degradation ####
rich_inf <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  summarise(count = n_distinct(epibiont))
View(rich_inf)

View(abund_div_inf)

abund_div_inf <- abund_div_inf %>%
  add_column(richness = rich_inf$count, .after = "simpson_inf")

pielou_inf <- abund_div_inf$shannon_inf/log(abund_div_inf$richness)

abund_div_inf <- abund_div_inf %>%
  add_column(pielou_inf, .after = "richness")

abund_div_inf <- abund_div_inf[!is.na(abund_div_inf$pielou_inf),]

# Summarize evenness per position ####
rich_pos <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, position) %>%
  summarise(count = n_distinct(epibiont))
View(rich_pos)

View(abund_div_pos)

abund_div_pos <- abund_div_pos %>%
  add_column(richness = rich_pos$count, .after = "simpson_pos")

pielou_pos <- abund_div_pos$shannon_pos/log(abund_div_pos$richness)

abund_div_pos <- abund_div_pos %>%
  add_column(pielou_pos, .after = "richness")

abund_div_pos <- abund_div_pos[!is.na(abund_div_pos$pielou_pos),]




# 4.2. Visualizing data set ==================================================================

# BARPLOT - Pielou's evenness x Length range ####

# For each quadrat
pielou_length1 <- abund_div_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  summarise(
    pielou_mean = mean(pielou_length),
    pielou_sd = sd(pielou_length)
  )
View(pielou_length1)

barplot_pielou_length1 <- ggplot(pielou_length1, aes(x = length_range, y = pielou_mean,
                                                       col = length_range, fill = length_range)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_length1

# Across all quadrats
pielou_length2 <- abund_div_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    pielou_mean = mean(pielou_length),
    pielou_sd = sd(pielou_length)
  )
View(pielou_length2)

barplot_pielou_length2 <- ggplot(pielou_length2, aes(x = length_range, y = pielou_mean, 
                                                       color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Pielou's evenness (J)")
barplot_pielou_length2

# BARPLOT - Pielou's evenness x Infestation ####

# For each quadrat
pielou_inf1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  summarise(
    pielou_mean = mean(pielou_inf),
    pielou_sd = sd(pielou_inf)
  )
View(pielou_inf1)

barplot_pielou_inf1 <- ggplot(pielou_inf1, aes(x = inf.lvl, y = pielou_mean,
                                                 col = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_inf1

# Across all quadrats
pielou_inf2 <- abund_div_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    pielou_mean = mean(pielou_inf),
    pielou_sd = sd(pielou_inf)
  )
View(pielou_inf2)

barplot_pielou_inf2 <- ggplot(pielou_inf2, aes(x = inf.lvl, y = pielou_mean, 
                                                 color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean Pielou's evenness (J)")
barplot_pielou_inf2

# BARPLOT - Pielou's evenness x Degradation ####

# For each quadrat
pielou_deg1 <- abund_div_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  summarise(
    pielou_mean = mean(pielou_inf),
    pielou_sd = sd(pielou_inf)
  )
View(pielou_deg1)

barplot_pielou_deg1 <- ggplot(pielou_deg1, aes(x = inf.percent, y = pielou_mean,
                                                 col = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_deg1

# Across all quadrats
pielou_deg2 <- abund_div_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    pielou_mean = mean(pielou_inf),
    pielou_sd = sd(pielou_inf)
  )
View(pielou_deg2)

barplot_pielou_deg2 <- ggplot(pielou_deg2, aes(x = inf.percent, y = pielou_mean, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean Pielou's evenness (J)")
barplot_pielou_deg2

# BARPLOT - Pielou's evenness x Position ####

# For each quadrat
pielou_pos1 <- abund_div_pos %>%
  dplyr::group_by(quadrat, position) %>%
  summarise(
    pielou_mean = mean(pielou_pos),
    pielou_sd = sd(pielou_pos)
  )
View(pielou_pos1)

barplot_pielou_pos1 <- ggplot(pielou_pos1, aes(x = position, y = pielou_mean,
                                                 col = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_pos1

# Across all quadrats
pielou_pos2 <- abund_div_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    pielou_mean = mean(pielou_pos),
    pielou_sd = sd(pielou_pos)
  )
View(pielou_pos2)

barplot_pielou_pos2 <- ggplot(pielou_pos2, aes(x = position, y = pielou_mean, 
                                                 color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean Pielou's evenness (J)")
barplot_pielou_pos2





# 3.3. MODELS: Pielou's evenness vs Length range ============================

# Investigate data set ####
View(abund_div_length)
hist(abund_div_length$pielou_length) # Normalish distribution
# Simpson's diversity is considered for each mussel - because each valve has the same length

# Fit model: LM ####
pielou_length_lm <- lm(pielou_length ~ length_range, data = abund_div_length)
summary(pielou_length_lm)
Anova(pielou_length_lm)
AIC(pielou_length_lm) # AIC(lm) = -154.4858

# Check model: LM
par(mfrow=c(2,2))
plot(pielou_length_lm) # Structure in the Residuals vs Leverage, distribution is normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = pielou_length_lm, plot = T) # Not too bad of a fit.

# Fit model: GLM (Gamma) ####
pielou_length_glm <- glm(pielou_length ~ length_range, data = abund_div_length,
                    family = "Gamma")
summary(pielou_length_glm)
Anova(pielou_length_glm)
AIC(pielou_length_glm) # AIC(glm) = -81.19874

# Check models
DHARMa::simulateResiduals(fittedModel = pielou_length_glm, plot = T) # Not a good fit at all.

# Pairwise comparisons ####
summary(glht(pielou_length_lm, linfct = mcp(length_range = "Tukey")))





# 3.4. MODELS: Pielou's evenness vs (Infestation x Degradation) ============================

# Investigate data set ####
View(abund_div_inf)
hist(abund_div_inf$pielou_inf) # Normalish distribution
# Simpson's diversity is considered for each mussel - because each valve has the same length

# Fit model: LM saturated ####
pielou_inf_lm <- lm(pielou_inf ~ inf.lvl * inf.percent, data = abund_div_inf)
summary(pielou_inf_lm)
Anova(pielou_inf_lm)
AIC(pielou_inf_lm) # AIC(lm) = -137.163

# Check model: LM saturated
par(mfrow=c(2,2))
plot(pielou_inf_lm) # Structure in the Residuals vs Leverage, distribution is normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = pielou_inf_lm, plot = T) # Not too bad of a fit.

# Fit model: LM additive ####
pielou_inf_lm1 <- lm(pielou_inf ~ inf.lvl + inf.percent, data = abund_div_inf)
summary(pielou_inf_lm1)
Anova(pielou_inf_lm1)
AIC(pielou_inf_lm1) # AIC(lm) = -149.8097

# Check model: LM saturated
par(mfrow=c(2,2))
plot(pielou_inf_lm1) # Structure in the Residuals vs Leverage, distribution is normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = pielou_inf_lm1, plot = T) # Not too bad of a fit.

# Fit model: GLM (Gamma) additive ####
pielou_inf_glm <- glm(pielou_inf ~ inf.lvl + inf.percent, data = abund_div_inf,
                         family = "Gamma")
summary(pielou_inf_glm)
Anova(pielou_inf_glm)
AIC(pielou_inf_glm) # AIC(glm) = -2.006067

# Check models
DHARMa::simulateResiduals(fittedModel = pielou_inf_glm, plot = T) # Not a good fit at all.





# 3.5. MODELS: Pielou's evenness vs Position ============================

# Investigate data set ####
View(abund_div_pos)
hist(abund_div_pos$pielou_pos) # Normalish distribution
# Simpson's diversity is considered for each mussel - because each valve has the same length

# Fit model: LM ####
pielou_pos_lm <- lm(pielou_pos ~ position, data = abund_div_pos)
summary(pielou_pos_lm)
Anova(pielou_pos_lm)
AIC(pielou_pos_lm) # AIC(lm) = -154.4858

# Check model: LM
par(mfrow=c(2,2))
plot(pielou_pos_lm) # Structure in the Residuals vs Leverage, distribution is normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = pielou_pos_lm, plot = T) # Not too bad of a fit.

# Fit model: GLM (Gamma) additive ####
pielou_pos_glm <- glm(pielou_pos ~ position, data = abund_div_pos,
                      family = "Gamma")
summary(pielou_pos_glm)
Anova(pielou_pos_glm)
AIC(pielou_pos_glm) # AIC(glm) = -80.6068

# Check models
DHARMa::simulateResiduals(fittedModel = pielou_pos_glm, plot = T) # Worst.

