## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between euendolithic infestation
##              and epibiosis in the brown mussel Perna perna
##              Script 1 - Shell parameters (including graphs) 
##
##
## Purpose of script: 
##      1.1. Infestation x Degradation - Pearson X2 test
##      1.2. Shell length vs Infestation x Degradation - LM
##      1.3. Epibiosis vs Shell length or Infestation x Degradation - GLM binomial
##
## Author: Alexia Dievart
##
## Date Created: 2023-02-22
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



###############################################################################################
# Section 1.1: Infestation x Degradation ------------------------------------------------------
###############################################################################################

# Prepare dataset =============================================================================

# Summarize infestation and degradation levels for each valve of each specimen in each quadrat
infxdeg <- epibio_perna %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    n = n_distinct(specimen)
  )
View(infxdeg)

# Create a contingency table
table(infxdeg$inf.lvl, infxdeg$inf.percent)



# Visualizing data ============================================================================

# Summarize infestation and degradation for each valve of each specimen
infxdeg1 <- infxdeg %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    N = sum(n)
  )
View(infxdeg1)

# Histogram of the number of valves for each shell degradation and infestation levels #########
hist_infxdeg1 <- ggplot(infxdeg1, aes (x = inf.lvl, y = N, fill = inf.percent)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = col_deg) +
  labs(y = "Proportion of mussel valves", x = "Infestation levels",
       fill = "SD Index") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.text=element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
hist_infxdeg1



# Pearson Chisq Test ==========================================================================
chisq.test(infxdeg$inf.lvl, infxdeg$inf.percent, simulate.p.value = TRUE)
# Infestation and degradation are independent variables
chisq.test(infxdeg$inf.lvl, infxdeg$inf.percent) # To obtain the degrees of freedom.





###############################################################################################
# Section 1.2: Shell length and (euendolithic infestation x shell degradation) ----------------
###############################################################################################

# Summarizing data ============================================================================

# Summarize shell length, infestation and degradation, for each valve of each specimen in each quadrat
length_infxdeg <- epibio_perna %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    shell_length = mean(shell.length)
  )
View(length_infxdeg)

# Summarize shell length, infestation and degradation for each quadrat
length_infxdeg1 <- length_infxdeg %>%
  dplyr::group_by(quadrat, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length),
    n = n(),
    se_length  = sd_length / sqrt(n)
  )
View(length_infxdeg1)

# Summarize shell length, infestation and degradation across all quadrats
length_infxdeg2 <- length_infxdeg %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length)
  )
View(length_infxdeg2)

# Summarize shell length and infestation across all quadrats
length_infxdeg3 <- length_infxdeg %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length)
  )
View(length_infxdeg3)

# Summarize shell length and degradation across all quadrats
length_infxdeg4 <- length_infxdeg %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length)
  )
View(length_infxdeg4)



# Visualizing data =========================================================================

# Barplot of mean shell length vs infestation x degradation, for each quadrat ####
barplot_length_infxdeg1 <- ggplot(length_infxdeg1, aes(x = inf.percent, y = mean_length, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - se_length, ymax = mean_length + se_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ quadrat) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Degradation", y = "Mean shell length")
barplot_length_infxdeg1

# Barplot of mean shell length vs infestation x degradation, across all quadrats ####
barplot_length_infxdeg2 <- ggplot(length_infxdeg2, aes(x = inf.percent, y = mean_length, 
                                           color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap(~ inf.lvl, ncol = 5) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = col_deg) +
  scale_fill_manual(values = col_deg) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Degradation", y = "Mean shell length")
barplot_length_infxdeg2

# Barplot of mean shell length vs infestation, across all quadrats ####
barplot_length_inf3 <- ggplot(length_infxdeg3, aes(x = inf.lvl, y = mean_length, 
                                                       color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = c(75, 80, 95, 105, 100),
           label = c("A", "A", "A", "B", "B"),
           size = 4) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Infestation", y = "Mean shell length")
barplot_length_inf3

# Barplot of mean shell length vs degradation, across all quadrats ####
barplot_length_inf4 <- ggplot(length_infxdeg4, aes(x = inf.percent, y = mean_length, 
                                                   color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  scale_color_manual(values = col_deg) +
  scale_fill_manual(values = col_deg) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Degradation", fill = "Degradation",
       x = "Degradation", y = "Mean shell length")
barplot_length_inf4





# MODELS: Shell length vs Infestation x Degradation =========================================

# Linear Model ##############################################################################

# Assumptions
hist(length_infxdeg$shell_length) # Data is normally distributed (but only positive)
plot(shell_length ~ inf.lvl * inf.percent, data = length_infxdeg)
# Relationship between shell length and infestation looks exponential, and variances are not constant.
# Same for the relationship between shell length and degradation.

# Fit model: LM saturated ###################################################################
length_lm <- lm(shell_length ~ inf.lvl * inf.percent, data = length_infxdeg)
summary(length_lm)
Anova(length_lm) # Interaction is nearly significant
AIC(length_lm) # AIC(lm) = 7232.53
# This is equivalent to a classic ANOVA with an unbalanced design.

# Check model: LM saturated
par(mfrow=c(2,2))
plot(length_lm) # Structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = length_lm, plot = T) # Quantile deviations detected
simulationOutput <- DHARMa::simulateResiduals(fittedModel = length_lm, plot = F) 
plotResiduals(simulationOutput, length_infxdeg$inf.lvl)  # within-group deviation detected for A
plotResiduals(simulationOutput, length_infxdeg$inf.percent) # within-group deviation detected for 1

# Fit model: LM additive ###################################################################
length_lm1 <- lm(shell_length ~ inf.lvl + inf.percent, data = length_infxdeg)
summary(length_lm1)
Anova(length_lm1)
AIC(length_lm1) # AIC(lm1) = 7228.259 < AIC(lm)

# Check model: LM additive
par(mfrow=c(2,2))
plot(length_lm1) # No structure in the Residuals vs Leverage
par(mfrow=c(1,1))

DHARMa::simulateResiduals(fittedModel = length_lm1, plot = T) # Quantile deviations detected
simulationOutput <- DHARMa::simulateResiduals(fittedModel = length_lm1, plot = F) 
plotResiduals(simulationOutput, length_infxdeg$inf.lvl)  # within-group deviation detected for A
plotResiduals(simulationOutput, length_infxdeg$inf.percent) # within-group deviation detected for 1

# Compare saturated and additive LMs #######################################################
anova(length_lm, length_lm1) # No real statistical difference
# I chose the model with the best AIC (lm1)

# Pairwise comparisons with the LM additive model (better AIC, length_lm1) #################
length_posthoc_inf <- glht(length_lm, linfct = mcp(inf.lvl = 'Tukey'))
summary(length_posthoc_inf)

length_posthoc_deg <- glht(length_lm1, linfct = mcp(inf.percent = 'Tukey'))
summary(length_posthoc_deg)





###############################################################################################
# Section 1.3: Epibiosis vs (shell length range x degradation x infestation x position) -------
###############################################################################################

# Summarizing data ============================================================================
View(epibio_perna)

# Create shell length ranges categories ####
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
epibio_perna$length_range <- cut(epibio_perna$shell.length, breaks = length_range_bins,
                                      labels = length_range_labels)

# Summarize prevalence of epibiosis for each shell length range, degradation and infestation for each quadrat ####
prev_epibio <- epibio_perna %>%
  dplyr::group_by(quadrat, specimen, valve, length_range, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    epibiosis = sum(abundance)
  )
View(prev_epibio)

prev_epibio$epibiosis[prev_epibio$epibiosis != 0] <- 1
prev_epibio$epibiosis <- as.factor(prev_epibio$epibiosis)

# Summarize prevalence of epibiosis and shell length for each degradation and infestation ####
prev_epibio_wide1 <- epibio_perna %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell.length),
    epibiosis = sum(abundance)
  )
View(prev_epibio_wide1)

prev_epibio_wide1$epibiosis[prev_epibio_wide1$epibiosis != 0] <- 1
prev_epibio_wide1$epibiosis <- as.factor(prev_epibio_wide1$epibiosis)





# Visualizing data =========================================================================

# Density plot of the prevalence of epibionts 
prev_dens <- ggplot(data = prev_epibio_wide1, aes(x = mean_length, y = after_stat(count), group = epibiosis)) +
  geom_density(aes(fill = epibiosis), position = position_fill()) +
  scale_fill_manual(values = colors_epibiosis, labels = c("No", "Yes")) +
  scale_x_continuous(expand = expansion(0),
                     breaks = seq(from = 30, to = 112, by = 10)) +
  scale_y_continuous(expand = expansion(0)) +
  labs(y = "Proportion of mussel valves", x = "Shell length (in mm)",
       fill = "Epibiosis") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
prev_dens




# MODELS: Epibiosis vs Length range / Infestation x Degradation ============================

# Fit model: GLM (binomial) w/ length range ####
prev_glm3 <- glm(epibiosis ~ length_range, data = prev_epibio,
                 family = "binomial")
summary(prev_glm3)
anova(prev_glm3, test = "Chisq")
AIC(prev_glm3) # AIC(glm2) = 516.4496

# Check models
DHARMa::simulateResiduals(fittedModel = prev_glm3, plot = T)

# Post-hoc comparisons
prev_glm3_emm <- emmeans(prev_glm3, pairwise ~ length_range)
prev_glm3_emm

# Fit model: GLM (binomial) w/ inf.lvl + inf.percent ####
prev_glm4 <- glm(epibiosis ~ inf.lvl + inf.percent, data = prev_epibio,
                 family = "binomial")
summary(prev_glm4)
anova(prev_glm4, test = "Chisq")
AIC(prev_glm4) # AIC(glm2) = 563.1809

# Check models
DHARMa::simulateResiduals(fittedModel = prev_glm4, plot = T)

# Post-hoc comparisons
epibio_posthoc_inf <- glht(prev_glm4, linfct = mcp(inf.lvl = 'Tukey'))
summary(epibio_posthoc_inf)

epibio_posthoc_deg <- glht(prev_glm4, linfct = mcp(inf.percent = 'Tukey'))
summary(epibio_posthoc_deg)





