## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between euendolithic infestation
##              and epibiosis in the brown mussel Perna perna
##              Script 2 - Community descriptors 1 (abundance, percent cover and biomass)
##
##
## Purpose of script: 
##
##
## Author: Alexia Dievart
##
## Date Created: 2023-02-23
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
View(abund_long)

# Isolate PERCENT COVER =======================================================================
cover_long <- epibio_perna1[,c(1:15, 18)]
View(cover_long)

# Isolate BIOMAS ==============================================================================
biom_long <- epibio_perna1[,c(1:15, 19)]
View(biom_long)





###############################################################################################
# Section 1: ABUNDANCE ------------------------------------------------------------------------
###############################################################################################

# 1.1. Prepare data set =======================================================================

# Create shell length ranges categories ####
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
abund_long$length_range <- cut(abund_long$shell.length, breaks = length_range_bins,
                                 labels = length_range_labels)

# Summarize total abundance for each length range ####
abund_epibio_length <- abund_long %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_epibio_length)

# Summarize total abundance for each infestation and degradation category ####
abund_epibio_inf <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_epibio_inf)

# Summarize total abundance for each position category ####
abund_epibio_pos <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, position) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_epibio_pos)

abund_epibio_pos <- abund_epibio_pos[! abund_epibio_pos$position == "byssus",]





# 1.2. Visualizing data =====================================================================

# BARPLOT - Abundance x Shell length range ####

# For each quadrat
abund_epibio_length1 <- abund_epibio_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_length1)

barplot_abund_length1 <- ggplot(abund_epibio_length1, aes(x = length_range, y = mean_abund, 
                                            color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total abundance per mussel")
barplot_abund_length1

# Across all quadrats
abund_epibio_length2 <- abund_epibio_length %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_length2)

barplot_abund_length2 <- ggplot(abund_epibio_length2, aes(x = length_range, y = mean_abund, 
                                            color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(15, 20, 20, 30, 40, 45, 50, 60, 50),
           label = c("A", "ABC", "ABCD", "BDE", "BEF", "BFG", "BG", "BG", "BG"),
           size = 4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total abundance per mussel")
barplot_abund_length2

# BARPLOT - Abundance x Infestation ####

# For each quadrat
abund_epibio_inf1 <- abund_epibio_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_inf1)

barplot_abund_inf1 <- ggplot(abund_epibio_inf1, aes(x = inf.lvl, y = mean_abund, 
                                                          color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total abundance per valve")
barplot_abund_inf1

# Across all quadrats
abund_epibio_inf2 <- abund_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_inf2)

barplot_abund_inf2 <- ggplot(abund_epibio_inf2, aes(x = inf.lvl, y = mean_abund, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total abundance per valve")
barplot_abund_inf2

# BARPLOT - Abundance x Degradation ####

# For each quadrat
abund_epibio_deg1 <- abund_epibio_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_deg1)

barplot_abund_deg1 <- ggplot(abund_epibio_deg1, aes(x = inf.percent, y = mean_abund, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean total abundance per valve")
barplot_abund_deg1

# Across all quadrats
abund_epibio_deg2 <- abund_epibio_inf %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_deg2)

barplot_abund_deg2 <- ggplot(abund_epibio_deg2, aes(x = inf.percent, y = mean_abund, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean total abundance per valve")
barplot_abund_deg2

# BARPLOT - Abundance x Position ####

# For each quadrat
abund_epibio_pos1 <- abund_epibio_pos %>%
  dplyr::group_by(quadrat, position) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_pos1)

barplot_abund_pos1 <- ggplot(abund_epibio_pos1, aes(x = position, y = mean_abund, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean total abundance per valve")
barplot_abund_pos1

# Across all quadrats
abund_epibio_pos2 <- abund_epibio_pos %>%
  dplyr::group_by(position) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_pos2)

barplot_abund_pos2 <- ggplot(abund_epibio_pos2, aes(x = position, y = mean_abund, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean total abundance per valve")
barplot_abund_pos2





# 1.3. MODELS: Abundance vs Length range ============================

# Investigate data set ####
View(abund_epibio_length)
hist(abund_epibio_length$abundance) # Poisson distribution
# Abundance is considered for each mussel - because each valve has the same length

# Fit model: LM ####
abund_length_lm <- lm(abundance ~ length_range, data = abund_epibio_length)
summary(abund_length_lm)
Anova(abund_length_lm)
AIC(abund_length_lm) # AIC(lm) = 3757.68

# Check model: LM saturated
par(mfrow=c(2,2))
plot(abund_length_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = abund_length_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (poisson) w/ length range ####
abund_length_glm <- glm(abundance ~ length_range, data = abund_epibio_length,
                 family = "poisson")
summary(abund_length_glm)
Anova(abund_length_glm)
AIC(abund_length_glm) # AIC(glm) = 6371.033

# Check models
DHARMa::simulateResiduals(fittedModel = abund_length_glm, plot = T) # Not a good fit at all.

# Fit model: GLM (quasi-poisson) w/ length range ####
abund_length_glm1 <- glm(abundance ~ length_range, data = abund_epibio_length,
                        family = "quasipoisson")
summary(abund_length_glm1)
Anova(abund_length_glm1)
AIC(abund_length_glm1) # Does not work with quasi-models
# I tried to calculate a qAIC but none of the methods I used worked. 

# Check models
DHARMa::simulateResiduals(fittedModel = abund_length_glm1, plot = T) # Quasipoisson not implemented
gratia::appraise(abund_length_glm1, method = "simulate") # Not a good fit.

# Fit model : GLM (negative binomial) w/ length range ####
abund_length_glm2 <- glm.nb(abundance ~ length_range, data = abund_epibio_length)
summary(abund_length_glm2)
Anova(abund_length_glm2)
AIC(abund_length_glm2) # AIC(glm) = 3209.095

# Check models
DHARMa::simulateResiduals(fittedModel = abund_length_glm2, plot = T) # Much better fit than the other models.
gratia::appraise(abund_length_glm2, method = "simulate")

# Pairwise comparison
summary(glht(abund_length_glm2, linfct = mcp(length_range = "Tukey")))





# 1.4. MODELS: Abundance vs Infestation x Degradation ============================

# Investigate data set ####
View(abund_epibio_inf)
hist(abund_epibio_inf$abundance) # Poisson distribution
# Abundance is considered for each mussel valve - because each valve has different parameters

# Fit model: LM saturated ####
abund_inf_lm <- lm(abundance ~ inf.lvl * inf.percent, data = abund_epibio_inf)
summary(abund_inf_lm)
Anova(abund_inf_lm)
AIC(abund_inf_lm) # AIC(lm) = 8515.926

# Check model: LM saturated
par(mfrow=c(2,2))
plot(abund_inf_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = abund_inf_lm, plot = T) # Not a good fit at all.

# Fit model: LM additive ####
abund_inf_lm1 <- lm(abundance ~ inf.lvl + inf.percent, data = abund_epibio_inf)
summary(abund_inf_lm1)
Anova(abund_inf_lm1)
AIC(abund_inf_lm1) # AIC(lm) = 5941.921

# Check model: LM additive
par(mfrow=c(2,2))
plot(abund_inf_lm1) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = abund_inf_lm1, plot = T) # Not a good fit at all.

# Fit model: GLM (poisson) saturated ####
abund_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_epibio_inf,
                        family = "poisson")
summary(abund_inf_glm)
Anova(abund_inf_glm)
AIC(abund_inf_glm) # AIC(glm) = 8935.121

# Check models
DHARMa::simulateResiduals(fittedModel = abund_inf_glm, plot = T) # Not a good fit at all.

# Fit model: GLM (poisson) additive ####
abund_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_epibio_inf,
                     family = "poisson")
summary(abund_inf_glm1)
Anova(abund_inf_glm1)
AIC(abund_inf_glm1) # AIC(glm) = 9051.979

# Check models
DHARMa::simulateResiduals(fittedModel = abund_inf_glm1, plot = T) # Not a good fit at all.

# Fit model : GLM (negative binomial) saturated ####
abund_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_epibio_inf)
summary(abund_inf_glm2)
Anova(abund_inf_glm2)
AIC(abund_inf_glm2) # AIC(glm) = 4963.625

# Check models
DHARMa::simulateResiduals(fittedModel = abund_inf_glm2, plot = T) # Better fit.
gratia::appraise(abund_inf_glm2, method = "simulate")

# Fit model : GLM (negative binomial) saturated ####
abund_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_epibio_inf)
summary(abund_inf_glm3)
Anova(abund_inf_glm3)
AIC(abund_inf_glm3) # AIC(glm) = 4958.362

# Check models
DHARMa::simulateResiduals(fittedModel = abund_inf_glm3, plot = T) # Better fit.
gratia::appraise(abund_inf_glm3, method = "simulate")

simulationOutput <- DHARMa::simulateResiduals(fittedModel = abund_inf_glm3, plot = T)
plotResiduals(simulationOutput, abund_epibio_inf$inf.lvl)  # within-group deviations significant for C
plotResiduals(simulationOutput, abund_epibio_inf$inf.percent) # within-group deviations for 5

# Pairwise comparison
summary(glht(abund_inf_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(abund_inf_glm3, linfct = mcp(inf.percent = "Tukey")))





# 1.5. MODELS: Abundance vs Position ============================

# Investigate data set ####
View(abund_epibio_pos)
hist(abund_epibio_pos$abundance) # Poisson distribution
# Abundance is considered for each mussel valve - because each valve has different parameters

# Fit model: LM ####
abund_pos_lm <- lm(abundance ~ position, data = abund_epibio_pos)
summary(abund_pos_lm)
Anova(abund_pos_lm)
AIC(abund_pos_lm) # AIC(lm) = 8470.069

# Check model: LM
par(mfrow=c(2,2))
plot(abund_pos_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = abund_pos_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (poisson) w/ position ####
abund_pos_glm <- glm(abundance ~ position, data = abund_epibio_pos,
                        family = "poisson")
summary(abund_pos_glm)
Anova(abund_pos_glm)
AIC(abund_pos_glm) # AIC(glm) = 10372.65

# Check models
DHARMa::simulateResiduals(fittedModel = abund_pos_glm, plot = T) # Worst !

# Fit model : GLM (negative binomial) w/ length range ####
abund_pos_glm2 <- glm.nb(abundance ~ position, data = abund_epibio_pos)
summary(abund_pos_glm2)
Anova(abund_pos_glm2)
AIC(abund_pos_glm2) # AIC(glm) = 6842.088

# Check models
DHARMa::simulateResiduals(fittedModel = abund_pos_glm2, plot = T) # Much better fit than the other models.
gratia::appraise(abund_pos_glm2, method = "simulate")

# Pairwise comparison
summary(glht(abund_pos_glm2, linfct = mcp(position = "Tukey")))





###############################################################################################
# Section 2: PERCENT COVER --------------------------------------------------------------------
###############################################################################################

# 2.1. Prepare data set =======================================================================

# Create shell length ranges categories ####
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
cover_long$length_range <- cut(cover_long$shell.length, breaks = length_range_bins,
                               labels = length_range_labels)

# Exclude data ####
cover_long1 <- cover_long[! cover_long$position == "byssus",] # Exclude epibionts in byssus
cover_long1 <- cover_long1[! cover_long1$position == "II epibiosis",] # Exclude IIary epibiosis
cover_long1 <- cover_long1[! is.na(cover_long1$cover.percent),] # Exclude NA
cover_long1 <- cover_long1[! cover_long1$mobility == "mobile",] # Exclude mobile organisms
View(cover_long1)
glimpse(cover_long1)

# Summarize total percent cover for each length range ####
cover_epibio_length <- cover_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_epibio_length)

cover_epibio_length1 <- cover_epibio_length %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  dplyr::summarise(
    cover = mean(cover)
  )
View(cover_epibio_length1)
# Per mussel and not per valve.

# Summarize total percent cover for each infestation and degradation category ####
cover_epibio_inf <- cover_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_epibio_inf)

cover_epibio_inf1 <- cover_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_epibio_inf1)

# Summarize total percent cover for each position category ####
cover_epibio_pos <- cover_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_epibio_pos)

# 2.2. Visualizing data =====================================================================

# BARPLOT - Percent cover x Shell length range ####

# For each quadrat
cover_epibio_length2 <- cover_epibio_length1 %>%
  dplyr::group_by(quadrat, length_range) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_length2)

barplot_cover_length2 <- ggplot(cover_epibio_length2, aes(x = length_range, y = mean_cover, 
                                                          color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel")
barplot_cover_length2

# Across all quadrats
cover_epibio_length3 <- cover_epibio_length1 %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_length3)

barplot_cover_length3 <- ggplot(cover_epibio_length3, aes(x = length_range, y = mean_cover, 
                                                          color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total percent cover per mussel")
barplot_cover_length3

# BARPLOT - Percent cover x Infestation ####

# For each quadrat
cover_epibio_inf2 <- cover_epibio_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_inf2)

barplot_cover_inf2 <- ggplot(cover_epibio_inf2, aes(x = inf.lvl, y = mean_cover, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean percent cover per valve")
barplot_cover_inf1

# Across all quadrats
cover_epibio_inf3 <- cover_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_inf3)

barplot_cover_inf3 <- ggplot(cover_epibio_inf3, aes(x = inf.lvl, y = mean_cover, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean percent cover per mussel valve")
barplot_cover_inf3

# BARPLOT - Percent cover x Degradation ####

# For each quadrat
cover_epibio_deg1 <- cover_epibio_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_deg1)

barplot_cover_deg1 <- ggplot(cover_epibio_deg1, aes(x = inf.percent, y = mean_cover, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean percent cover per valve")
barplot_cover_deg1

# Across all quadrats
cover_epibio_deg2 <- cover_epibio_inf %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_deg2)

barplot_cover_deg2 <- ggplot(cover_epibio_deg2, aes(x = inf.percent, y = mean_cover, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean percent cover per valve")
barplot_cover_deg2

# BARPLOT - Percent cover x Position ####

# For each quadrat
cover_epibio_pos1 <- cover_epibio_pos %>%
  dplyr::group_by(quadrat, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_pos1)

barplot_cover_pos1 <- ggplot(cover_epibio_pos1, aes(x = position, y = mean_cover, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean percent cover per valve")
barplot_cover_pos1

# Across all quadrats
cover_epibio_pos2 <- cover_epibio_pos %>%
  dplyr::group_by(position) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_pos2)

barplot_cover_pos2 <- ggplot(cover_epibio_pos2, aes(x = position, y = mean_cover, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_cover - sd_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean percent cover per valve")
barplot_cover_pos2





# 2.3. MODELS: Percent cover vs Length range ============================

# Investigate data set ####
View(cover_epibio_length1)
hist(cover_epibio_length1$cover) # Gamma distribution
# Percent cover is considered for each mussel

# Fit model: LM ####
cover_length_lm <- lm(cover ~ length_range, data = cover_epibio_length1)
summary(cover_length_lm)
Anova(cover_length_lm)
AIC(cover_length_lm) # AIC(lm) = 3650.644

# Check model: LM saturated
par(mfrow=c(2,2))
plot(cover_length_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = cover_length_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGamma) w/ log-link ####
library(glmmTMB)

cover_length_tmb <- glmmTMB(cover ~ length_range,
                    ziformula = ~ length_range,
                    family = ziGamma(link = "log"),
                    data = cover_epibio_length1)
summary(cover_length_tmb)
Anova(cover_length_tmb)
AIC(cover_length_tmb) # AIC(glmTMB) = 2829.021

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = cover_length_tmb, plot = T) # Better fit, but variance non homogeneous
# Tried same models with logit-link, that was a worst match.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = cover_length_tmb, plot = T)
plotResiduals(simulationOutput, cover_epibio_length1$length_range)  # within-group deviations for 40, 50 and 60

# Pairwise comparison
summary(glht(cover_length_tmb, linfct = mcp(length_range = "Tukey")))





# 2.3. MODELS: Percent cover vs Infestation x Degradation x Position ============================

# Investigate data set ####
View(cover_epibio_inf1)
hist(cover_epibio_inf1$cover) # Gamma distribution
# Percent cover is considered for each valve

# Fit model: LM saturated ####
cover_inf_lm <- lm(cover ~ inf.lvl * inf.percent * position, data = cover_epibio_inf1)
summary(cover_inf_lm)
Anova(cover_inf_lm)
AIC(cover_inf_lm) # AIC(lm) = 8846.096

# Check model: LM saturated
par(mfrow=c(2,2))
plot(cover_inf_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = cover_inf_lm, plot = T) # Not a good fit at all.

# Fit model: LM additive w/ some interactions ####
cover_inf_lm1 <- lm(cover ~ inf.lvl + inf.percent + position + inf.lvl:position + inf.percent:position,
                   data = cover_epibio_inf1)
summary(cover_inf_lm1)
Anova(cover_inf_lm1)
AIC(cover_inf_lm1) # AIC(lm) = 8816.593

# Check model: LM additive
par(mfrow=c(2,2))
plot(cover_inf_lm1) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = cover_inf_lm1, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGamma) saturated w/ log-link ####
cover_inf_tmb <- glmmTMB(cover ~ inf.lvl + inf.percent + position + inf.lvl:position + inf.percent:position,
                         ziformula = ~ inf.lvl + inf.percent + position,
                         family = ziGamma(link = "log"),
                         data = cover_epibio_inf1)
summary(cover_inf_tmb)
Anova(cover_inf_tmb)
AIC(cover_inf_tmb) # AIC(glmTMB) = 6040.311

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = cover_inf_tmb, plot = T) # Better fit, but variance non homogeneous
# Tried same models with logit-link, that was a worst match.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = cover_inf_tmb, plot = T)
plotResiduals(simulationOutput, cover_epibio_inf1$inf.lvl)  # within-group deviations for A and B
plotResiduals(simulationOutput, cover_epibio_inf1$inf.percent)  # within-group deviations for 1 to 3

# Fit model: GLM (ziGamma) saturated w/ log-link - inf.percent:position ####
cover_inf_tmb1 <- glmmTMB(cover ~ inf.lvl + inf.percent + position + inf.percent:position,
                            ziformula = ~ inf.lvl + inf.percent,
                            family = ziGamma(link = "log"),
                            data = cover_epibio_inf1)
summary(cover_inf_tmb1)
Anova(cover_inf_tmb1)
AIC(cover_inf_tmb1) # AIC(glmTMB) = 6039.732

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = cover_inf_tmb1, plot = T) # Better fit, but variance non homogeneous
# Tried same models with logit-link, that was a worst match.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = cover_inf_tmb1, plot = T)
plotResiduals(simulationOutput, cover_epibio_inf1$inf.lvl)  # within-group deviations for A and B
plotResiduals(simulationOutput, cover_epibio_inf1$inf.percent)  # within-group deviations not significant

# Fit model: GLM (ziGamma) saturated w/ log-link - all interactions ####
cover_inf_tmb2 <- glmmTMB(cover ~ inf.lvl + inf.percent + position,
                          ziformula = ~ inf.lvl + inf.percent + position,
                          family = ziGamma(link = "log"),
                          data = cover_epibio_inf1)
summary(cover_inf_tmb2)
Anova(cover_inf_tmb2)
AIC(cover_inf_tmb2) # AIC(glmTMB) = 6057.466

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = cover_inf_tmb2, plot = T) # Better fit, but variance non homogeneous
# Tried same models with logit-link, that was a worst match.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = cover_inf_tmb2, plot = T)
plotResiduals(simulationOutput, cover_epibio_inf1$inf.lvl)  # within-group deviations for A and B
plotResiduals(simulationOutput, cover_epibio_inf1$inf.percent)  # within-group deviations for 1 to 3

# Fit model: GLM (ziGamma) saturated w/ log-link - all interactions - position ####
cover_inf_tmb3 <- glmmTMB(cover ~ inf.lvl + inf.percent,
                          ziformula = ~ inf.lvl + inf.percent,
                          family = ziGamma(link = "log"),
                          data = cover_epibio_inf1)
summary(cover_inf_tmb3)
Anova(cover_inf_tmb3)
AIC(cover_inf_tmb3) # AIC(glmTMB) = 6099.168

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = cover_inf_tmb3, plot = T) # Better fit, but variance non homogeneous
# Tried same models with logit-link, that was a worst match.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = cover_inf_tmb3, plot = T)
plotResiduals(simulationOutput, cover_epibio_inf1$inf.lvl)  # within-group deviations for A, B and D
plotResiduals(simulationOutput, cover_epibio_inf1$inf.percent)  # within-group deviations for 1 to 3

# Compare models ####
anova(cover_inf_tmb, cover_inf_tmb1, cover_inf_tmb2, cover_inf_tmb3)
anova(cover_inf_tmb, cover_inf_tmb1)
# Choose models with interactions - cover_inf_tmb1

# Pairwise comparisons ####
emm_cover_inf_tmb <- emmeans(cover_inf_tmb1, pairwise ~ inf.percent * position)
emm_cover_inf_tmb





###############################################################################################
# Section 3: BIOMASS --------------------------------------------------------------------------
###############################################################################################

# 3.1. Prepare data set =======================================================================
View(biom_long)

# Create shell length ranges categories ####
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
biom_long$length_range <- cut(biom_long$shell.length, breaks = length_range_bins,
                               labels = length_range_labels)

# Exclude data ####
biom_long1 <- biom_long[! biom_long$position == "byssus",] # Exclude epibionts in byssus
biom_long1 <- biom_long1[! is.na(biom_long1$biomass),] # Exclude NA
View(biom_long1)
glimpse(biom_long1)

# Summarize total biomass for each length range ####
biom_epibio_length <- biom_long1 %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_epibio_length)

# Summarize total biomass for each infestation and degradation category ####
biom_epibio_inf <- biom_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_epibio_inf)

# Summarize total biomass for each position category ####
biom_epibio_pos <- biom_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_epibio_pos)




# 3.2. Visualizing data =====================================================================

# BARPLOT - Biomass x Shell length range ####

# For each quadrat
biom_epibio_length1 <- biom_epibio_length %>%
  dplyr::group_by(quadrat, length_range) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_length1)

barplot_biom_length1 <- ggplot(biom_epibio_length1, aes(x = length_range, y = mean_biom, 
                                                          color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Total biomass per mussel valve")
barplot_biom_length1

# Across all quadrats
biom_epibio_length2 <- biom_epibio_length %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_length2)

barplot_biom_length2 <- ggplot(biom_epibio_length2, aes(x = length_range, y = mean_biom, 
                                                          color = length_range, fill = length_range)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total abundance per mussel valve")
barplot_biom_length2

# BARPLOT - Biomass x Infestation ####

# For each quadrat
biom_epibio_inf1 <- biom_epibio_inf %>%
  dplyr::group_by(quadrat, inf.lvl) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_inf1)

barplot_biom_inf1 <- ggplot(biom_epibio_inf1, aes(x = inf.lvl, y = mean_biom, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - sd_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total biomass per valve")
barplot_biom_inf1

# Across all quadrats
biom_epibio_inf2 <- biom_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_inf2)

barplot_biom_inf2 <- ggplot(biom_epibio_inf2, aes(x = inf.lvl, y = mean_biom, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total biomass per mussel valve")
barplot_biom_inf2

# BARPLOT - Biomass x Degradation ####

# For each quadrat
biom_epibio_deg1 <- biom_epibio_inf %>%
  dplyr::group_by(quadrat, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_deg1)

barplot_biom_deg1 <- ggplot(biom_epibio_deg1, aes(x = inf.percent, y = mean_biom, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - sd_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean total biomass per valve")
barplot_biom_deg1

# Across all quadrats
biom_epibio_deg2 <- biom_epibio_inf %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_deg2)

barplot_biom_deg2 <- ggplot(biom_epibio_deg2, aes(x = inf.percent, y = mean_biom, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean total biomass per valve")
barplot_biom_deg2

# BARPLOT - Biomass x Position ####

# For each quadrat
biom_epibio_pos1 <- biom_epibio_pos %>%
  dplyr::group_by(quadrat, position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_pos1)

barplot_biom_pos1 <- ggplot(biom_epibio_pos1, aes(x = position, y = mean_biom, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - sd_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat, ncol = 7) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean total biomass per valve")
barplot_biom_pos1

# Across all quadrats
biom_epibio_pos2 <- biom_epibio_pos %>%
  dplyr::group_by(position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_pos2)

barplot_biom_pos2 <- ggplot(biom_epibio_pos2, aes(x = position, y = mean_biom, 
                                                    color = position, fill = position)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean percent cover per valve")
barplot_biom_pos2




# 3.3. MODELS: Biomass vs Length range ============================

# Investigate data set ####
View(biom_epibio_length)
hist(biom_epibio_length$biomass) # Gamma distribution
hist(log(biom_epibio_length$biomass))
# Biomass is considered for each mussel - because each valve has the same length

# Fit model: LM ####
biom_length_lm <- lm(biomass ~ length_range, data = biom_epibio_length)
summary(biom_length_lm)
Anova(biom_length_lm)
AIC(biom_length_lm) # AIC(lm) = 6669.401

# Check model: LM saturated
par(mfrow=c(2,2))
plot(biom_length_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = biom_length_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGamma) w/ log-link ####
library(glmmTMB)

biom_length_tmb <- glmmTMB(biomass ~ length_range,
                            ziformula = ~ length_range,
                            family = ziGamma(link = "log"),
                            data = biom_epibio_length)
summary(biom_length_tmb)
Anova(biom_length_tmb)
AIC(biom_length_tmb) # AIC(glmTMB) = 4156.659

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = biom_length_tmb, plot = T) # Better fit, but variance non homogeneous
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_length_tmb, plot = T)
plotResiduals(simulationOutput, biom_length_tmb$length_range)  # within-group deviations for 120

# Pairwise comparisons
summary(glht(biom_length_tmb, linfct = mcp(length_range = "Tukey")))






# 3.3. MODELS: Biomass vs Infestation x Degradation============================

# Investigate data set ####
View(biom_epibio_inf)
hist(biom_epibio_inf$biomass) # Gamma distribution
# Biomass is considered for each valve

# Fit model: LM saturated ####
biom_inf_lm <- lm(biomass ~ inf.lvl * inf.percent, data = biom_epibio_inf)
summary(biom_inf_lm)
Anova(biom_inf_lm)
AIC(biom_inf_lm) # AIC(lm) = 10874.55

# Check model: LM saturated
par(mfrow=c(2,2))
plot(biom_inf_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = biom_inf_lm, plot = T) # Not a good fit at all.

# Fit model: LM additive ####
biom_inf_lm1 <- lm(biomass ~ inf.lvl + inf.percent, data = biom_epibio_inf)
summary(biom_inf_lm1)
Anova(biom_inf_lm1)
AIC(biom_inf_lm1) # AIC(lm) = 10858.7

# Check model: LM saturated
par(mfrow=c(2,2))
plot(biom_inf_lm1) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = biom_inf_lm1, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGamma) w/ log-link ####
biom_inf_tmb <- glmmTMB(biomass ~ inf.lvl + inf.percent,
                          ziformula = ~ inf.lvl + inf.percent,
                          family = ziGamma(link = "log"),
                          data = biom_epibio_inf)
summary(biom_inf_tmb)
Anova(biom_inf_tmb)
AIC(biom_inf_tmb) # AIC(glmTMB) = 6326.979

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = biom_inf_tmb, plot = T) # Better fit, but variance non homogeneous
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_inf_tmb, plot = T)
plotResiduals(simulationOutput, biom_epibio_inf$inf.lvl)  # within-group deviations for D
plotResiduals(simulationOutput, biom_epibio_inf$inf.percent)  # within-group deviations for 4

# Pairwise comparisons
summary(glht(biom_inf_tmb, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(biom_inf_tmb, linfct = mcp(inf.percent = "Tukey")))





# 3.3. MODELS: Biomass vs Position ============================

# Investigate data set ####
View(biom_epibio_pos)
hist(biom_epibio_pos$biomass) # Gamma distribution
# Biomass is considered for each valve

# Fit model: LM ####
biom_pos_lm <- lm(biomass ~ position, data = biom_epibio_pos)
summary(biom_pos_lm)
Anova(biom_pos_lm)
AIC(biom_pos_lm) # AIC(lm) = 16473.79

# Check model: LM saturated
par(mfrow=c(2,2))
plot(biom_pos_lm) # Structure in the Residuals vs Leverage, distribution is not normal
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = biom_pos_lm, plot = T) # Not a good fit at all.

# Fit model: GLM (ziGamma) w/ log-link ####
biom_pos_tmb <- glmmTMB(biomass ~ position,
                        ziformula = ~ position,
                        family = ziGamma(link = "log"),
                        data = biom_epibio_pos)
summary(biom_pos_tmb)
Anova(biom_pos_tmb)
AIC(biom_pos_tmb) # AIC(glmTMB) = 8547.751

# Check model: GLM (ziGamma)
DHARMa::simulateResiduals(fittedModel = biom_pos_tmb, plot = T) # Better fit, but variance non homogeneous
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_pos_tmb, plot = T)
plotResiduals(simulationOutput, biom_pos_tmb$position)  # within-group deviations for D

# Pairwise comparisons
summary(glht(biom_pos_tmb, linfct = mcp(position = "Tukey")))
