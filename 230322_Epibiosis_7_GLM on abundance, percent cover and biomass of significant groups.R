## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##              Script 7 - Subsequent analysises between epibiotic groups and species
##
##
## Purpose of script: 
##   
##
## Author: Alexia Dievart
##
## Date Created: 2023-01-24
## Dates Updated:
##
## Copyright (c) Alexia DIEVART 2023
## Email: alexia.dievart@hotmail.fr
##
## ---------------------------
##
## Notes: 
##   - 
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
               DHARMa,
               lme4,
               mgcv,
               gratia,
               ggplot2,
               ggtext,
               tidyr,
               dplyr,
               MuMIn,
               glue,
               glmmTMB,
               ggeffects,
               ggpubr,
               ggrepel,
               rstatix,
               multcomp,
               car,
               MASS,
               emmeans,
               vegan,
               mvabund,
               MuMIn,
               gllvm,
               corrplot,
               gclus, 
               lme4)

# Set default ggplot theme
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  legend.position = "none"))

# Create colour vectors
colors_inflvl <- c("coral4", "brown3", "indianred1", "pink3", "darkgrey") # Inf. lvl: A, B, C, D and E
colors_epibiosis <- c("darkseagreen", "cadetblue3")
colors_degrad <- c("chocolate4", "chocolate", "darkgoldenrod3", "darkgoldenrod1", "bisque2", "cornsilk2")
colors_position <- c("tomato4", "cadetblue", "darkgrey", "brown3")
colors_position1 <- c("darkgrey", "brown3")
colors_position2 <- c("cadetblue", "darkgrey", "brown3")





###############################################################################################
# Section 1: Load and clean data ----------------------------------------------------------------
###############################################################################################

# Load raw data 
epibio <- read.csv("epibiosis.csv", h=T, dec=",", sep=";")
head(epibio)

# Clean raw dataset
factor_epibio <- c("location", "quadrat", "specimen", "mussel.species", "sex", "valve", "inf.lvl", "inf.percent",
                   "higher.group", "epibiont", "position", "basibiont")
epibio[,factor_epibio] <- lapply(epibio[,factor_epibio], factor)
epibio$inf.percent <- recode_factor(epibio$inf.percent, "0" = "1", "0-25" = "2", "25-50" = "3", "50-75" = "4",
                                    "75-100" = "5", "100" = "6")
epibio <- epibio[! epibio$mussel.species == "Mytilus galloprovincialis",]
epibio <- epibio[! epibio$abundance == 0,]
View(epibio)

# Create a new "shell length range" variable in the raw dataset
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
epibio$length_range <- cut(epibio$shell.length, breaks = length_range_bins,
                           labels = length_range_labels)






###############################################################################################
# Section 2: Abundance of higher groups -------------------------------------------------------
###############################################################################################

# Summarize abundance of significant epibiotic higher groups
abund_long <- epibio[,c(1:16,22)]
View(abund_long)

abund_group_length <- abund_long[abund_long$higher.group %in% c("Bryozoa", "Cirripedia", "Ochrophyta", "Sedentaria",
                                                                "Bivalvia", "Rhodophyta"),]
View(abund_group_length)





# 2.1. Cirripedia ==============================================================================

# 2.1.1. Abundance x Shell length ####
abund_cirri <- abund_group_length[abund_group_length$higher.group %in% "Cirripedia",]

abund_cirri1 <- abund_cirri %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_cirri1)

hist(abund_cirri1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
cirriab_length_aov <- aov(abundance ~ length_range, data = abund_cirri1)
summary(cirriab_length_aov)

# Check ANOVA assumptions
plot(cirriab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_cirri1) # But, variances are actually homogeneous

plot(cirriab_length_aov, 2) # Normal distribution of the data NOT assumed
cirriab_aov_residuals <- residuals(object = cirriab_length_aov)
shapiro.test(x = cirriab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
cirriab_length_glm <- glm(abundance ~ length_range, data = abund_cirri1,
                       family = "poisson")
summary(cirriab_length_glm)
Anova(cirriab_length_glm)
AIC(cirriab_length_glm) # AIC(glm) = 3443.539

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_length_glm, plot = T) # Not a good fit.
gratia::appraise(cirriab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
cirriab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_cirri1)
summary(cirriab_length_glm1)
Anova(cirriab_length_glm1)
AIC(cirriab_length_glm1) # AIC(glm) = 2092.168

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_length_glm1, plot = T) # Much better fit.
gratia::appraise(cirriab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Cirripedia ####
summary(glht(cirriab_length_glm1, linfct = mcp(length_range = "Tukey")))






# 2.1.2. Abundance x (Infestation + Degradation) ####
abund_cirri <- abund_group_length[abund_group_length$higher.group %in% "Cirripedia",]

abund_cirri2 <- abund_cirri %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_cirri2)

hist(abund_cirri2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
cirriab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_cirri2)
summary(cirriab_inf_aov)

# Check ANOVA assumptions
plot(cirriab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_cirri2) # But, variances are actually homogeneous

plot(cirriab_inf_aov, 2) # Normal distribution of the data NOT assumed
cirriab_inf_residuals <- residuals(object = cirriab_inf_aov)
shapiro.test(x = cirriab_inf_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
cirriab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_cirri2)
summary(cirriab_inf_aov1)

# Check ANOVA assumptions
plot(cirriab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_cirri2) # But, variances are actually homogeneous

plot(cirriab_inf_aov1, 2) # Normal distribution of the data NOT assumed
cirriab_inf_residuals1 <- residuals(object = cirriab_inf_aov1)
shapiro.test(x = cirriab_inf_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
cirriab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_cirri2,
                          family = "poisson")
summary(cirriab_inf_glm)
Anova(cirriab_inf_glm)
AIC(cirriab_inf_glm) # AIC(glm) = 3826.816

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(cirriab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
cirriab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_cirri2,
                       family = "poisson")
summary(cirriab_inf_glm1)
Anova(cirriab_inf_glm1)
AIC(cirriab_inf_glm1) # AIC(glm) = 3858.379

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(cirriab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
cirriab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_cirri2)
summary(cirriab_inf_glm2)
Anova(cirriab_inf_glm2)
AIC(cirriab_inf_glm2) # AIC(glm) = 2718.501

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_inf_glm2, plot = T) # Much better fit.
gratia::appraise(cirriab_inf_glm2, method = "simulate")

# Model #6: Additive GLM (negative binomial) ####
cirriab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_cirri2)
summary(cirriab_inf_glm3)
Anova(cirriab_inf_glm3)
AIC(cirriab_inf_glm3) # AIC(glm) = 2711.114

# Check models
DHARMa::simulateResiduals(fittedModel = cirriab_inf_glm3, plot = T) # Much better fit.
gratia::appraise(cirriab_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in Cirripedia ####
summary(glht(cirriab_inf_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(cirriab_inf_glm3, linfct = mcp(inf.percent = "Tukey")))





# 2.2. Ochrophyta ==============================================================================

# 2.2.1. Abundance x Shell length ####
abund_ochro <- abund_group_length[abund_group_length$higher.group %in% "Ochrophyta",]

abund_ochro1 <- abund_ochro %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_ochro1)

hist(abund_ochro1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
ochroab_length_aov <- aov(abundance ~ length_range, data = abund_ochro1)
summary(ochroab_length_aov)

# Check ANOVA assumptions
plot(ochroab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_ochro1) # But, variances are actually homogeneous

plot(ochroab_length_aov, 2) # Normal distribution of the data NOT assumed
ochroab_aov_residuals <- residuals(object = ochroab_length_aov)
shapiro.test(x = ochroab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
ochroab_length_glm <- glm(abundance ~ length_range, data = abund_ochro1,
                          family = "poisson")
summary(ochroab_length_glm)
Anova(ochroab_length_glm)
AIC(ochroab_length_glm) # AIC(glm) = 480.4857

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_length_glm, plot = T) # Not a bad fit.
gratia::appraise(ochroab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
ochroab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_ochro1)
summary(ochroab_length_glm1)
Anova(ochroab_length_glm1)
AIC(ochroab_length_glm1) # AIC(glm) = 479.94

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_length_glm1, plot = T) # Much better fit.
gratia::appraise(ochroab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Ochrophyta ####
summary(glht(ochroab_length_glm1, linfct = mcp(length_range = "Tukey")))










# 2.2.2. Abundance x (Infestation + Degradation) ####
abund_ochro <- abund_group_length[abund_group_length$higher.group %in% "Ochrophyta",]

abund_ochro2 <- abund_ochro %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_ochro2)

hist(abund_ochro2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
ochroab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_ochro2)
summary(ochroab_inf_aov)

# Check ANOVA assumptions
plot(ochroab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_ochro2) # But, variances are NOT homogeneous

plot(ochroab_inf_aov, 2) # Normal distribution of the data NOT assumed
ochroab_inf_aov_residuals <- residuals(object = ochroab_inf_aov)
shapiro.test(x = ochroab_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
ochroab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_ochro2)
summary(ochroab_inf_aov1)

# Check ANOVA assumptions
plot(ochroab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_ochro2) # But, variances are NOT homogeneous

plot(ochroab_inf_aov1, 2) # Normal distribution of the data NOT assumed
ochroab_inf_aov_residuals1 <- residuals(object = ochroab_inf_aov1)
shapiro.test(x = ochroab_inf_aov_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
ochroab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_ochro2,
                          family = "poisson")
summary(ochroab_inf_glm)
Anova(ochroab_inf_glm)
AIC(ochroab_inf_glm) # AIC(glm) = 557.9235

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(ochroab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
ochroab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_ochro2,
                       family = "poisson")
summary(ochroab_inf_glm1)
Anova(ochroab_inf_glm1)
AIC(ochroab_inf_glm1) # AIC(glm) = 552.5534

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(ochroab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
ochroab_length_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_ochro2)
summary(ochroab_length_glm2)
Anova(ochroab_length_glm2)
AIC(ochroab_length_glm2) # AIC(glm) = 559.9251

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_length_glm2, plot = T) # Much better fit.
gratia::appraise(ochroab_length_glm2, method = "simulate")

# Model #6: Saturated GLM (negative binomial) ####
ochroab_length_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_ochro2)
summary(ochroab_length_glm3)
Anova(ochroab_length_glm3)
AIC(ochroab_length_glm3) # AIC(glm) = 554.5544

# Check models
DHARMa::simulateResiduals(fittedModel = ochroab_length_glm3, plot = T) # Much better fit.
gratia::appraise(ochroab_length_glm2, method = "simulate")

# Pairwise comparisons on abundance in Ochrophyta ####
summary(glht(ochroab_length_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(ochroab_length_glm3, linfct = mcp(inf.percent = "Tukey")))





# 2.3. Sedentaria ==============================================================================
abund_seden <- abund_group_length[abund_group_length$higher.group %in% "Sedentaria",]

abund_seden1 <- abund_seden %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_seden1)

hist(abund_seden1$abundance) # Poisson distribution

# 2.3.1. Abundance x Shell length ####

# Model #1: Regular ANOVA ####
sedenab_length_aov <- aov(abundance ~ length_range, data = abund_seden1)
summary(sedenab_length_aov)

# Check ANOVA assumptions
plot(sedenab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_seden1) # But, variances are actually homogeneous

plot(sedenab_length_aov, 2) # Normal distribution of the data NOT assumed
sedenab_aov_residuals <- residuals(object = sedenab_length_aov)
shapiro.test(x = sedenab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
sedenab_length_glm <- glm(abundance ~ length_range, data = abund_seden1,
                          family = "poisson")
summary(sedenab_length_glm)
Anova(sedenab_length_glm)
AIC(sedenab_length_glm) # AIC(glm) = 3715.029

# Check models
DHARMa::simulateResiduals(fittedModel = sedenab_length_glm, plot = T) # Not a good fit.
gratia::appraise(sedenab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
sedenab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_seden1)
summary(sedenab_length_glm1)
Anova(sedenab_length_glm1)
AIC(sedenab_length_glm1) # AIC(glm) = 1384.658

# Check models
DHARMa::simulateResiduals(fittedModel = sedenab_length_glm1, plot = T) # Much better fit.
gratia::appraise(sedenab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Sedentaria ####
summary(glht(sedenab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 2.4. Bivalvia ==============================================================================
abund_bival <- abund_group_length[abund_group_length$higher.group %in% "Bivalvia",]

abund_bival1 <- abund_bival %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_bival1)

hist(abund_bival1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
bivalab_length_aov <- aov(abundance ~ length_range, data = abund_bival1)
summary(bivalab_length_aov)

# Check ANOVA assumptions
plot(bivalab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_bival1) # But, variances are actually homogeneous

plot(bivalab_length_aov, 2) # Normal distribution of the data NOT assumed
bivalab_aov_residuals <- residuals(object = bivalab_length_aov)
shapiro.test(x = bivalab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
bivalab_length_glm <- glm(abundance ~ length_range, data = abund_bival1,
                          family = "poisson")
summary(bivalab_length_glm)
Anova(bivalab_length_glm)
AIC(bivalab_length_glm) # AIC(glm) = 134.5676

# Check models
DHARMa::simulateResiduals(fittedModel = bivalab_length_glm, plot = T) # Not a good fit.
gratia::appraise(bivalab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
bivalab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_bival1)
summary(bivalab_length_glm1)
Anova(bivalab_length_glm1)
AIC(bivalab_length_glm1) # AIC(glm) = 136.5679

# Check models
DHARMa::simulateResiduals(fittedModel = bivalab_length_glm1, plot = T) # Much better fit.
gratia::appraise(bivalab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Bivalvia ####
summary(glht(bivalab_length_glm, linfct = mcp(length_range = "Tukey")))





# 2.5. Rhodophyta ==============================================================================

# 2.5.1. Abundance x Shell length ####
abund_rhodo <- abund_group_length[abund_group_length$higher.group %in% "Rhodophyta",]

abund_rhodo1 <- abund_rhodo %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_rhodo1)

hist(abund_rhodo1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
rhodoab_length_aov <- aov(abundance ~ length_range, data = abund_rhodo1)
summary(rhodoab_length_aov)

# Check ANOVA assumptions
plot(rhodoab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_rhodo1) # But, variances are actually homogeneous

plot(rhodoab_length_aov, 2) # Normal distribution of the data NOT assumed
rhodoab_aov_residuals <- residuals(object = rhodoab_length_aov)
shapiro.test(x = rhodoab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
rhodoab_length_glm <- glm(abundance ~ length_range, data = abund_rhodo1,
                          family = "poisson")
summary(rhodoab_length_glm)
Anova(rhodoab_length_glm)
AIC(rhodoab_length_glm) # AIC(glm) = 988.7863

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_length_glm, plot = T) # Not a good fit.
gratia::appraise(rhodoab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
rhodoab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_rhodo1)
summary(rhodoab_length_glm1)
Anova(rhodoab_length_glm1)
AIC(rhodoab_length_glm1) # AIC(glm) = 795.6389

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_length_glm1, plot = T) # Much better fit.
gratia::appraise(rhodoab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Rhodophyta ####
summary(glht(rhodoab_length_glm1, linfct = mcp(length_range = "Tukey")))






# 2.5.2. Abundance x (Infestation + Degradation) ####
abund_rhodo <- abund_group_length[abund_group_length$higher.group %in% "Rhodophyta",]

abund_rhodo2 <- abund_rhodo %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_rhodo2)

hist(abund_rhodo2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
rhodoab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_rhodo2)
summary(rhodoab_inf_aov)

# Check ANOVA assumptions
plot(rhodoab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_rhodo2) # But, variances are actually homogeneous

plot(rhodoab_inf_aov, 2) # Normal distribution of the data NOT assumed
rhodoab_inf_aov_residuals <- residuals(object = rhodoab_inf_aov)
shapiro.test(x = rhodoab_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
rhodoab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_rhodo2)
summary(rhodoab_inf_aov1)

# Check ANOVA assumptions
plot(rhodoab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl + inf.percent, data = abund_rhodo2) # But, variances are actually homogeneous

plot(rhodoab_inf_aov1, 2) # Normal distribution of the data NOT assumed
rhodoab_inf_aov_residuals1 <- residuals(object = rhodoab_inf_aov1)
shapiro.test(x = rhodoab_inf_aov_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
rhodoab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_rhodo2,
                          family = "poisson")
summary(rhodoab_inf_glm)
Anova(rhodoab_inf_glm)
AIC(rhodoab_inf_glm) # AIC(glm) = 1054.072

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(rhodoab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
rhodoab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_rhodo2,
                       family = "poisson")
summary(rhodoab_inf_glm1)
Anova(rhodoab_inf_glm1)
AIC(rhodoab_inf_glm1) # AIC(glm) = 1081.929

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(rhodoab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
rhodoab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_rhodo2)
summary(rhodoab_inf_glm2)
Anova(rhodoab_inf_glm2)
AIC(rhodoab_inf_glm2) # AIC(glm) = 946.8692

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_inf_glm2, plot = T) # Much better fit.
gratia::appraise(rhodoab_inf_glm2, method = "simulate")

# Model #6: Additive GLM (negative binomial) ####
rhodoab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_rhodo2)
summary(rhodoab_inf_glm3)
Anova(rhodoab_inf_glm3)
AIC(rhodoab_inf_glm3) # AIC(glm) = 950.0343

# Check models
DHARMa::simulateResiduals(fittedModel = rhodoab_inf_glm3, plot = T) # Much better fit.
gratia::appraise(rhodoab_inf_glm3, method = "simulate")

anova(rhodoab_inf_glm2, rhodoab_inf_glm3)

# Pairwise comparisons on abundance in Rhodophyta ####
rhodoab_inf_glm2.emm <- emmeans(rhodoab_inf_glm2, ~ inf.lvl | inf.percent)
pairs(rhodoab_inf_glm2.emm)




# 2.6. Bryozoa ==============================================================================
abund_bryo <- abund_group_length[abund_group_length$higher.group %in% "Bryozoa",]

abund_bryo1 <- abund_bryo %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_bryo1)

hist(abund_bryo1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
bryoab_length_aov <- aov(abundance ~ length_range, data = abund_bryo1)
summary(bryoab_length_aov)

# Check ANOVA assumptions
plot(bryoab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_bryo1) # But, variances are actually homogeneous

plot(bryoab_length_aov, 2) # Normal distribution of the data NOT assumed
bryoab_aov_residuals <- residuals(object = bryoab_length_aov)
shapiro.test(x = bryoab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
bryoab_length_glm <- glm(abundance ~ length_range, data = abund_bryo1,
                          family = "poisson")
summary(bryoab_length_glm)
Anova(bryoab_length_glm)
AIC(bryoab_length_glm) # AIC(glm) = 270.4471

# Check models
DHARMa::simulateResiduals(fittedModel = bryoab_length_glm, plot = T) # Not too bad of a fit.
gratia::appraise(bryoab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
bryoab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_bryo1)
summary(bryoab_length_glm1)
Anova(bryoab_length_glm1)
AIC(bryoab_length_glm1) # AIC(glm) = 272.4483

# Check models
DHARMa::simulateResiduals(fittedModel = bryoab_length_glm1, plot = T) # Much better fit.
gratia::appraise(bryoab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in Bryozoa (poisson) ####
summary(glht(bryoab_length_glm1, linfct = mcp(length_range = "Tukey")))






###############################################################################################
# Section 3: Abundance of epibiotic species ---------------------------------------------------
###############################################################################################

# Summarize abundance of significant epibiotic higher groups
abund_long <- epibio[,c(1:16,22)]
View(abund_long)

abund_sp_length <- abund_long[abund_long$epibiont %in% c("Tetraclita serrata", "Juvenile barnacle",
                                                           "Ralfsia verrucosa", "Spirorbis spp.",
                                                           "Corallinales"),]
View(abund_sp_length)

abund_sp_inf <- abund_long[abund_long$epibiont %in% c("Tetraclita serrata", "Juvenile barnacle",
                                                         "Ralfsia verrucosa", "Hildenbrandia lecannellieri",
                                                         "Corallinales"),]
View(abund_sp_inf)



# 3.1. Tetraclita serrata =================================================================

# 3.1.1. Abundance x Shell length ####
abund_tetser <- abund_sp_length[abund_sp_length$epibiont %in% "Tetraclita serrata",]

abund_tetser1 <- abund_tetser %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_tetser1)

hist(abund_tetser1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
tetserab_length_aov <- aov(abundance ~ length_range, data = abund_tetser1)
summary(tetserab_length_aov)

# Check ANOVA assumptions
plot(tetserab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_tetser1) # But, variances are actually homogeneous

plot(tetserab_length_aov, 2) # Normal distribution of the data NOT assumed
tetserab_aov_residuals <- residuals(object = tetserab_length_aov)
shapiro.test(x = tetserab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
tetserab_length_glm <- glm(abundance ~ length_range, data = abund_tetser1,
                         family = "poisson")
summary(tetserab_length_glm)
Anova(tetserab_length_glm)
AIC(tetserab_length_glm) # AIC(glm) = 785.797

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_length_glm, plot = T) # Not a good fit.
gratia::appraise(tetserab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
tetserab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_tetser1)
summary(tetserab_length_glm1)
Anova(tetserab_length_glm1)
AIC(tetserab_length_glm1) # AIC(glm) = 630.4342

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_length_glm1, plot = T) # Much better fit.
gratia::appraise(tetserab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in TETSER (negative binomial) ####
summary(glht(tetserab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 3.1.1. Abundance x (Infestation + Degradation) ####
abund_tetser <- abund_sp_length[abund_sp_length$epibiont %in% "Tetraclita serrata",]

abund_tetser2 <- abund_tetser %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent, epibiont) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_tetser2)

hist(abund_tetser2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
tetserab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_tetser2)
summary(tetserab_inf_aov)

# Check ANOVA assumptions
plot(tetserab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_tetser2) # But, variances are actually homogeneous

plot(tetserab_inf_aov, 2) # Normal distribution of the data NOT assumed
tetserab_inf_residuals <- residuals(object = tetserab_inf_aov)
shapiro.test(x = tetserab_inf_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
tetserab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_tetser2)
summary(tetserab_inf_aov1)

# Check ANOVA assumptions
plot(tetserab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_tetser2) # But, variances are actually homogeneous

plot(tetserab_inf_aov1, 2) # Normal distribution of the data NOT assumed
tetserab_inf_residuals1 <- residuals(object = tetserab_inf_aov1)
shapiro.test(x = tetserab_inf_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
tetserab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_tetser2,
                           family = "poisson")
summary(tetserab_inf_glm)
Anova(tetserab_inf_glm)
AIC(tetserab_inf_glm) # AIC(glm) = 851.9744

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(tetserab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
tetserab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_tetser2,
                        family = "poisson")
summary(tetserab_inf_glm1)
Anova(tetserab_inf_glm1)
AIC(tetserab_inf_glm1) # AIC(glm) = 854.9371

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(tetserab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
tetserab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_tetser2)
summary(tetserab_inf_glm2)
Anova(tetserab_inf_glm2)
AIC(tetserab_inf_glm2) # AIC(glm) = 731.6987

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_inf_glm2, plot = T) # Much better fit.
gratia::appraise(tetserab_inf_glm2, method = "simulate")

# Model #6: Additive GLM (negative binomial) ####
tetserab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_tetser2)
summary(tetserab_inf_glm3)
Anova(tetserab_inf_glm3)
AIC(tetserab_inf_glm3) # AIC(glm) = 724.6908

# Check models
DHARMa::simulateResiduals(fittedModel = tetserab_inf_glm3, plot = T) # Much better fit.
gratia::appraise(tetserab_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in TETSER (negative binomial) ####
summary(glht(tetserab_inf_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(tetserab_inf_glm3, linfct = mcp(inf.percent = "Tukey")))



# 3.2. Juvenile barnacle =================================================================

# 3.2.1. Abundance x Shell length ####
abund_juvbar <- abund_sp_length[abund_sp_length$epibiont %in% "Juvenile barnacle",]

abund_juvbar1 <- abund_juvbar %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_juvbar1)

hist(abund_juvbar1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
juvbarab_length_aov <- aov(abundance ~ length_range, data = abund_juvbar1)
summary(juvbarab_length_aov)

# Check ANOVA assumptions
plot(juvbarab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_juvbar1) # Variances are not homogeneous

plot(juvbarab_length_aov, 2) # Normal distribution of the data NOT assumed
juvbarab_aov_residuals <- residuals(object = juvbarab_length_aov)
shapiro.test(x = juvbarab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
juvbarab_length_glm <- glm(abundance ~ length_range, data = abund_juvbar1,
                           family = "poisson")
summary(juvbarab_length_glm)
Anova(juvbarab_length_glm)
AIC(juvbarab_length_glm) # AIC(glm) = 2578.478

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_length_glm, plot = T) # Not a good fit.
gratia::appraise(juvbarab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
juvbarab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_juvbar1)
summary(juvbarab_length_glm1)
Anova(juvbarab_length_glm1)
AIC(juvbarab_length_glm1) # AIC(glm) = 1675.637

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_length_glm1, plot = T) # Much better fit.
gratia::appraise(juvbarab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in JUVBAR (negative binomial) ####
summary(glht(juvbarab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 3.2.1. Abundance x (Infestation + Degradation) ####
abund_juvbar <- abund_sp_length[abund_sp_length$epibiont %in% "Juvenile barnacle",]

abund_juvbar2 <- abund_juvbar %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent, epibiont) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_juvbar2)

hist(abund_juvbar2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
juvbarab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_juvbar2)
summary(juvbarab_inf_aov)

# Check ANOVA assumptions
plot(juvbarab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_juvbar2) # Variances are actually homogeneous

plot(juvbarab_inf_aov, 2) # Normal distribution of the data NOT assumed
juvbarab_inf_residuals <- residuals(object = juvbarab_inf_aov)
shapiro.test(x = juvbarab_inf_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
juvbarab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_juvbar2)
summary(juvbarab_inf_aov1)

# Check ANOVA assumptions
plot(juvbarab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_juvbar2) # Variances are actually homogeneous

plot(juvbarab_inf_aov1, 2) # Normal distribution of the data NOT assumed
juvbarab_inf_residuals1 <- residuals(object = juvbarab_inf_aov1)
shapiro.test(x = juvbarab_inf_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
juvbarab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_juvbar2,
                           family = "poisson")
summary(juvbarab_inf_glm)
Anova(juvbarab_inf_glm)
AIC(juvbarab_inf_glm) # AIC(glm) = 2660.919

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(juvbarab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
juvbarab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_juvbar2,
                        family = "poisson")
summary(juvbarab_inf_glm1)
Anova(juvbarab_inf_glm1)
AIC(juvbarab_inf_glm1) # AIC(glm) = 2680.377

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(juvbarab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
juvbarab_length_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_juvbar2)
summary(juvbarab_length_glm2)
Anova(juvbarab_length_glm2)
AIC(juvbarab_length_glm2) # AIC(glm) = 2051.3

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_length_glm2, plot = T) # Much better fit.
gratia::appraise(juvbarab_length_glm2, method = "simulate")

# Model #6: Additive GLM (negative binomial) ####
juvbarab_length_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_juvbar2)
summary(juvbarab_length_glm3)
Anova(juvbarab_length_glm3)
AIC(juvbarab_length_glm3) # AIC(glm) = 2042.327

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarab_length_glm3, plot = T) # Much better fit.
gratia::appraise(juvbarab_length_glm3, method = "simulate")

# Pairwise comparisons on abundance in JUVBAR (negative binomial) ####
summary(glht(juvbarab_length_glm3, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(juvbarab_length_glm3, linfct = mcp(inf.percent = "Tukey")))





# 3.3. Spirorbis spp. =================================================================

# 3.3.1. Abundance x Shell length ####
abund_spiror <- abund_sp_length[abund_sp_length$epibiont %in% "Spirorbis spp.",]

abund_spiror1 <- abund_spiror %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_spiror1)

hist(abund_spiror1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
spirorab_length_aov <- aov(abundance ~ length_range, data = abund_spiror1)
summary(spirorab_length_aov)

# Check ANOVA assumptions
plot(spirorab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_spiror1) # Variances are actually homogeneous

plot(spirorab_length_aov, 2) # Normal distribution of the data NOT assumed
spirorab_aov_residuals <- residuals(object = spirorab_length_aov)
shapiro.test(x = spirorab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
spirorab_length_glm <- glm(abundance ~ length_range, data = abund_spiror1,
                           family = "poisson")
summary(spirorab_length_glm)
Anova(spirorab_length_glm)
AIC(spirorab_length_glm) # AIC(glm) = 3004.359

# Check models
DHARMa::simulateResiduals(fittedModel = spirorab_length_glm, plot = T) # Not a good fit.
gratia::appraise(spirorab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
spirorab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_spiror1)
summary(spirorab_length_glm1)
Anova(spirorab_length_glm1)
AIC(spirorab_length_glm1) # AIC(glm) = 1128.483

# Check models
DHARMa::simulateResiduals(fittedModel = spirorab_length_glm1, plot = T) # Much better fit.
gratia::appraise(spirorab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in SPIROR (negative binomial) ####
summary(glht(spirorab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 3.4. Corallinales =================================================================

# 3.4.1. Abundance x Shell length ####
abund_coral <- abund_sp_length[abund_sp_length$epibiont %in% "Corallinales",]

abund_coral1 <- abund_coral %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_coral1)

hist(abund_coral1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
coralab_length_aov <- aov(abundance ~ length_range, data = abund_coral1)
summary(coralab_length_aov)

# Check ANOVA assumptions
plot(coralab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_coral1) # Variances are actually homogeneous

plot(coralab_length_aov, 2) # Normal distribution of the data NOT assumed
coralab_aov_residuals <- residuals(object = coralab_length_aov)
shapiro.test(x = coralab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
coralab_length_glm <- glm(abundance ~ length_range, data = abund_coral1,
                           family = "poisson")
summary(coralab_length_glm)
Anova(coralab_length_glm)
AIC(coralab_length_glm) # AIC(glm) = 448.8019

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_length_glm, plot = T) # Not a bad fit.
gratia::appraise(coralab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
coralab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_coral1)
summary(coralab_length_glm1)
Anova(coralab_length_glm1)
AIC(coralab_length_glm1) # AIC(glm) = 446.8992

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_length_glm1, plot = T) # Much better fit.
gratia::appraise(coralab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in SPIROR (negative binomial) ####
summary(glht(coralab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 3.4.1. Abundance x (Infestation + Degradation) ####
abund_coral <- abund_sp_length[abund_sp_length$epibiont %in% "Corallinales",]

abund_coral2 <- abund_coral %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_coral2)

hist(abund_coral2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
coralab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_coral2)
summary(coralab_inf_aov)

# Check ANOVA assumptions
plot(coralab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_coral2) # Variances are actually NOT homogeneous

plot(coralab_inf_aov, 2) # Normal distribution of the data NOT assumed
coralab_inf_aov_residuals <- residuals(object = coralab_inf_aov)
shapiro.test(x = coralab_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
coralab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_coral2)
summary(coralab_inf_aov1)

# Check ANOVA assumptions
plot(coralab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_coral2) # Variances are actually NOT homogeneous

plot(coralab_inf_aov1, 2) # Normal distribution of the data NOT assumed
coralab_inf_aov_residuals1 <- residuals(object = coralab_inf_aov1)
shapiro.test(x = coralab_inf_aov_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
coralab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_coral2,
                          family = "poisson")
summary(coralab_inf_glm)
Anova(coralab_inf_glm)
AIC(coralab_inf_glm) # AIC(glm) = 498.9419

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(coralab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
coralab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_coral2,
                       family = "poisson")
summary(coralab_inf_glm1)
Anova(coralab_inf_glm1)
AIC(coralab_inf_glm1) # AIC(glm) = 490.6389

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(coralab_inf_glm1, method = "simulate")

# Model #3: Saturated GLM (negative binomial) ####
coralab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_coral2)
summary(coralab_inf_glm2)
Anova(coralab_inf_glm2)
AIC(coralab_inf_glm2) # AIC(glm) = 500.943

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_inf_glm2, plot = T) # Much better fit.
gratia::appraise(coralab_length_glm1, method = "simulate")

# Model #3: Additive (negative binomial) ####
coralab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_coral2)
summary(coralab_inf_glm3)
Anova(coralab_inf_glm3)
AIC(coralab_inf_glm3) # AIC(glm) = 492.6401

# Check models
DHARMa::simulateResiduals(fittedModel = coralab_inf_glm3, plot = T) # Much better fit.
gratia::appraise(coralab_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in Corallinales (Poisson) ####
summary(glht(coralab_inf_glm1, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(coralab_inf_glm1, linfct = mcp(inf.percent = "Tukey")))





# 3.5. Ralfsia verrucosa =================================================================

# 3.5.1. Abundance x Shell length ####
abund_ralver <- abund_sp_length[abund_sp_length$epibiont %in% "Ralfsia verrucosa",]

abund_ralver1 <- abund_ralver %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_ralver1)

hist(abund_ralver1$abundance) # Poisson distribution

# Model #1: Regular ANOVA ####
ralverab_length_aov <- aov(abundance ~ length_range, data = abund_ralver1)
summary(ralverab_length_aov)

# Check ANOVA assumptions
plot(ralverab_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ length_range, data = abund_ralver1) # Variances are actually homogeneous

plot(ralverab_length_aov, 2) # Normal distribution of the data NOT assumed
ralverab_aov_residuals <- residuals(object = ralverab_length_aov)
shapiro.test(x = ralverab_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (poisson) w/ length range ####
ralverab_length_glm <- glm(abundance ~ length_range, data = abund_ralver1,
                          family = "poisson")
summary(ralverab_length_glm)
Anova(ralverab_length_glm)
AIC(ralverab_length_glm) # AIC(glm) = 464.8049

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_length_glm, plot = T) # Not a bad fit.
gratia::appraise(ralverab_length_glm, method = "simulate")

# Model #3: GLM (negative binomial) w/ length range ####
ralverab_length_glm1 <- glm.nb(abundance ~ length_range, data = abund_ralver1)
summary(ralverab_length_glm1)
Anova(ralverab_length_glm1)
AIC(ralverab_length_glm1) # AIC(glm) = 466.4278

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_length_glm1, plot = T) # Much better fit.
gratia::appraise(ralverab_length_glm1, method = "simulate")

# Pairwise comparisons on abundance in RALVER (negative binomial) ####
summary(glht(ralverab_length_glm1, linfct = mcp(length_range = "Tukey")))





# 3.5.2. Abundance x (Infestation + Degradation) ####
abund_ralver <- abund_sp_length[abund_sp_length$epibiont %in% "Ralfsia verrucosa",]

abund_ralver2 <- abund_ralver %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_ralver2)

hist(abund_ralver2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
ralverab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_ralver2)
summary(ralverab_inf_aov)

# Check ANOVA assumptions
plot(ralverab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_ralver2) # Variances are actually NOT homogeneous

plot(ralverab_inf_aov, 2) # Normal distribution of the data NOT assumed
ralverab_inf_aov_residuals <- residuals(object = ralverab_inf_aov)
shapiro.test(x = ralverab_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
ralverab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_ralver2)
summary(ralverab_inf_aov1)

# Check ANOVA assumptions
plot(ralverab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_ralver2) # Variances are actually NOT homogeneous

plot(ralverab_inf_aov1, 2) # Normal distribution of the data NOT assumed
ralverab_inf_aov_residuals1 <- residuals(object = ralverab_inf_aov1)
shapiro.test(x = ralverab_inf_aov_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
ralverab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_ralver2,
                           family = "poisson")
summary(ralverab_inf_glm)
Anova(ralverab_inf_glm)
AIC(ralverab_inf_glm) # AIC(glm) = 548.0595

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(ralverab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
ralverab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_ralver2,
                        family = "poisson")
summary(ralverab_inf_glm1)
Anova(ralverab_inf_glm1)
AIC(ralverab_inf_glm1) # AIC(glm) = 536.2538

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(ralverab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
ralverab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_ralver2)
summary(ralverab_inf_glm2)
Anova(ralverab_inf_glm2)
AIC(ralverab_inf_glm2) # AIC(glm) = 550.0615

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_inf_glm2, plot = T) # Not much better.
gratia::appraise(ralverab_inf_glm2, method = "simulate")

# Model #6: Additive GLM (negative binomial) ####
ralverab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_ralver2)
summary(ralverab_inf_glm3)
Anova(ralverab_inf_glm3)
AIC(ralverab_inf_glm3) # AIC(glm) = 538.2558

# Check models
DHARMa::simulateResiduals(fittedModel = ralverab_inf_glm3, plot = T) # Not much better.
gratia::appraise(ralverab_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in RALVER (additive Poisson) ####
summary(glht(ralverab_inf_glm1, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(ralverab_inf_glm1, linfct = mcp(inf.percent = "Tukey")))





# 3.6. Hildenbrandia lecannellieri ============================================================

# 3.6.2. Abundance x (Infestation + Degradation) ####
abund_hillec <- abund_sp_inf[abund_sp_inf$epibiont %in% "Hildenbrandia lecannellieri",]
View(abund_hillec)

abund_hillec2 <- abund_hillec %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_hillec2)

hist(abund_hillec2$abundance) # Poisson distribution

# Model #1: Saturated ANOVA ####
hillecab_inf_aov <- aov(abundance ~ inf.lvl * inf.percent, data = abund_hillec2)
summary(hillecab_inf_aov)

# Check ANOVA assumptions
plot(coralab_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_hillec2) # Variances are actually homogeneous

plot(hillecab_inf_aov, 2) # Normal distribution of the data NOT assumed
hillecab_inf_aov_residuals <- residuals(object = hillecab_inf_aov)
shapiro.test(x = hillecab_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Additive ANOVA ####
hillecab_inf_aov1 <- aov(abundance ~ inf.lvl + inf.percent, data = abund_hillec2)
summary(hillecab_inf_aov1)

# Check ANOVA assumptions
plot(hillecab_inf_aov1, 1) # Homogeneity of variances NOT assumed
leveneTest(abundance ~ inf.lvl * inf.percent, data = abund_hillec2) # Variances are actually homogeneous

plot(hillecab_inf_aov1, 2) # Normal distribution of the data NOT assumed
hillecab_inf_aov_residuals1 <- residuals(object = hillecab_inf_aov1)
shapiro.test(x = hillecab_inf_aov_residuals1)  # Residuals not normally distributed

# Model #3: Saturated GLM (poisson) ####
hillecab_inf_glm <- glm(abundance ~ inf.lvl * inf.percent, data = abund_hillec2,
                       family = "poisson")
summary(hillecab_inf_glm)
Anova(hillecab_inf_glm)
AIC(hillecab_inf_glm) # AIC(glm) = 223.8799

# Check models
DHARMa::simulateResiduals(fittedModel = hillecab_inf_glm, plot = T) # Not a good fit.
gratia::appraise(hillecab_inf_glm, method = "simulate")

# Model #4: Additive GLM (poisson) ####
hillecab_inf_glm1 <- glm(abundance ~ inf.lvl + inf.percent, data = abund_hillec2,
                        family = "poisson")
summary(hillecab_inf_glm1)
Anova(hillecab_inf_glm1)
AIC(hillecab_inf_glm1) # AIC(glm) = 217.9964

# Check models
DHARMa::simulateResiduals(fittedModel = hillecab_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(hillecab_inf_glm1, method = "simulate")

# Model #5: Saturated GLM (negative binomial) ####
hillecab_inf_glm2 <- glm.nb(abundance ~ inf.lvl * inf.percent, data = abund_hillec2)
summary(hillecab_inf_glm2)
Anova(coralab_inf_glm2)
AIC(hillecab_inf_glm2) # AIC(glm) = 225.8808

# Check models
DHARMa::simulateResiduals(fittedModel = hillecab_inf_glm2, plot = T) # Much better fit.
gratia::appraise(hillecab_inf_glm2, method = "simulate")

# Model #6: Additive (negative binomial) ####
hillecab_inf_glm3 <- glm.nb(abundance ~ inf.lvl + inf.percent, data = abund_hillec2)
summary(hillecab_inf_glm3)
Anova(hillecab_inf_glm3)
AIC(hillecab_inf_glm3) # AIC(glm) = 219.9973

# Check models
DHARMa::simulateResiduals(fittedModel = hillecab_inf_glm3, plot = T) # Much better fit.
gratia::appraise(hillecab_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in HILLEC (Poisson) ####
summary(glht(hillecab_inf_glm1, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(hillecab_inf_glm1, linfct = mcp(inf.percent = "Tukey")))





###############################################################################################
# Section 4: Percent cover of higher groups -------------------------------------------------------
###############################################################################################

# Summarize percent cover of significant epibiotic higher groups
cover_long <- epibio[,c(1:15, 18,22)]
cover_long <- cover_long[! cover_long$cover.percent == 0,]
cover_long <- cover_long[! cover_long$position == "II epibiosis",]
cover_long <- cover_long[! cover_long$position == "byssus",]
View(cover_long)

cover_group_length <- cover_long[cover_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Sedentaria",
                                                                "Rhodophyta"),]
View(cover_group_length)





# 4.1. Cirripedia ==============================================================================

# 4.1.1. Percent cover x Shell length ####
cover_cirri <- cover_group_length[cover_group_length$higher.group %in% "Cirripedia",]

cover_cirri1 <- cover_cirri %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_cirri1)

hist(cover_cirri1$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
cirricov_length_aov <- aov(cover ~ length_range, data = cover_cirri1)
summary(cirricov_length_aov)

# Check ANOVA assumptions
plot(cirricov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_cirri1) # But, variances are actually NOT homogeneous

plot(cirricov_length_aov, 2) # Normal distribution of the data NOT assumed
cirricov_aov_residuals <- residuals(object = cirricov_length_aov)
shapiro.test(x = cirricov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
cirricov_length_glm <- glm(cover ~ length_range, data = cover_cirri1,
                          family = "Gamma")
summary(cirricov_length_glm)
Anova(cirricov_length_glm)
AIC(cirricov_length_glm) # AIC(glm) = 1393.08

# Check models
DHARMa::simulateResiduals(fittedModel = cirricov_length_glm, plot = T) # Not a bad fit.
gratia::appraise(cirricov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Cirripedia ####
summary(glht(cirricov_length_glm, linfct = mcp(length_range = "Tukey")))





# 4.1.2. Percent cover x (Infestation x Degradation x Position) ####
cover_cirri <- cover_group_length[cover_group_length$higher.group %in% "Cirripedia",]

cover_cirri2 <- cover_cirri %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_cirri2)

hist(cover_cirri2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
cirricov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_cirri2)
summary(cirricov_inf_aov)
AIC(cirricov_inf_aov) # AIC(aov) = 4015.124

# Check ANOVA assumptions
plot(cirricov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_cirri2) # But, variances are actually NOT homogeneous

plot(cirricov_inf_aov, 2) # Normal distribution of the data NOT assumed
cirricov_inf_aov_residuals <- residuals(object = cirricov_inf_aov)
shapiro.test(x = cirricov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
cirricov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_cirri2,
                           family = "Gamma")
summary(cirricov_inf_glm)
Anova(cirricov_inf_glm) # Significant: inf.percent and inf.percent:position
AIC(cirricov_inf_glm) # AIC(glm) = 1884.403

# Check models
DHARMa::simulateResiduals(fittedModel = cirricov_inf_glm, plot = T) # Not a good fit.
gratia::appraise(cirricov_length_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ interaction ####
cirricov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.percent:position, data = cover_cirri2,
                        family = "Gamma")
summary(cirricov_inf_glm1)
Anova(cirricov_inf_glm1) # Significant: inf.percent and inf.percent:position
AIC(cirricov_inf_glm1) # AIC(glm) = 1870.042

# Check models
DHARMa::simulateResiduals(fittedModel = cirricov_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(cirricov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) w/out interaction ####
cirricov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_cirri2,
                         family = "Gamma")
summary(cirricov_inf_glm2)
Anova(cirricov_inf_glm2) # Significant: inf.percent
AIC(cirricov_inf_glm2) # AIC(glm) = 1873.484

# Check models
DHARMa::simulateResiduals(fittedModel = cirricov_inf_glm2, plot = T) # Not a good fit.
gratia::appraise(cirricov_inf_glm2, method = "simulate")

anova(cirricov_inf_glm1, cirricov_inf_glm2, test="Chisq") # Nearly significant. 

# Pairwise comparisons on abundance in Cirripedia ####
summary(glht(cirricov_inf_glm2, linfct = mcp(inf.percent = "Tukey")))










# 4.2. Ochrophyta ==============================================================================
cover_ochro <- cover_group_length[cover_group_length$higher.group %in% "Ochrophyta",]

cover_ochro1 <- cover_ochro %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_ochro1)

hist(cover_ochro1$cover) # Gamma distribution

# 4.2.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
ochrocov_length_aov <- aov(cover ~ length_range, data = cover_ochro1)
summary(ochrocov_length_aov)

# Check ANOVA assumptions
plot(cirricov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_cirri1) # But, variances are actually NOT homogeneous

plot(ochrocov_length_aov, 2) # Normal distribution of the data NOT assumed
ochrocov_aov_residuals <- residuals(object = ochrocov_length_aov)
shapiro.test(x = ochrocov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
ochrocov_length_glm <- glm(cover ~ length_range, data = cover_ochro1,
                           family = "Gamma")
summary(ochrocov_length_glm)
Anova(ochrocov_length_glm)
AIC(ochrocov_length_glm) # AIC(glm) = 1103.188

# Check models
DHARMa::simulateResiduals(fittedModel = ochrocov_length_glm, plot = T) # A very good fit !
gratia::appraise(ochrocov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Ochrophyta ####
summary(glht(ochrocov_length_glm, linfct = mcp(length_range = "Tukey")))





# 4.2.2. Percent cover x (Infestation x Degradation x Position) ####
cover_ochro <- cover_group_length[cover_group_length$higher.group %in% "Ochrophyta",]

cover_ochro2 <- cover_ochro %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_ochro2)

hist(cover_ochro2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
ochrocov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_ochro2)
summary(ochrocov_inf_aov)

# Check ANOVA assumptions
plot(ochrocov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_ochro2) # But, variances are actually NOT homogeneous

plot(ochrocov_inf_aov, 2) # Normal distribution of the data NOT assumed
ochrocov_inf_aov_residuals <- residuals(object = ochrocov_inf_aov)
shapiro.test(x = ochrocov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
ochrocov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_ochro2,
                           family = "Gamma")
summary(ochrocov_inf_glm)
Anova(ochrocov_inf_glm)
AIC(ochrocov_inf_glm) # AIC(glm) = 1530.827

# Check models
DHARMa::simulateResiduals(fittedModel = ochrocov_inf_glm, plot = T) # Not a bad fit !
gratia::appraise(ochrocov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
ochrocov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.lvl:position + inf.percent:position, 
                        data = cover_ochro2,
                        family = "Gamma")
summary(ochrocov_inf_glm1)
Anova(ochrocov_inf_glm1)
AIC(ochrocov_inf_glm1) # AIC(glm) = 1521.468

# Check models
DHARMa::simulateResiduals(fittedModel = ochrocov_inf_glm1, plot = T) # Not a bad fit !
gratia::appraise(ochrocov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) w/ an interaction ####
ochrocov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position + inf.percent:position, 
                         data = cover_ochro2,
                         family = "Gamma")
summary(ochrocov_inf_glm2)
Anova(ochrocov_inf_glm2)
AIC(ochrocov_inf_glm2) # AIC(glm) = 1520.344

# Check models
DHARMa::simulateResiduals(fittedModel = ochrocov_inf_glm2, plot = T) # A very good fit.
gratia::appraise(ochrocov_inf_glm2, method = "simulate")

# Pairwise comparisons on percent cover in Ochrophyta ####
emmeans(ochrocov_inf_glm2, pairwise ~ position | inf.percent)





# 4.3. Rhodophyta ==============================================================================
cover_rhodo <- cover_group_length[cover_group_length$higher.group %in% "Rhodophyta",]

cover_rhodo1 <- cover_rhodo %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_rhodo1)

hist(cover_rhodo1$cover) # Gamma distribution

# 4.3.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
rhodocov_length_aov <- aov(cover ~ length_range, data = cover_rhodo1)
summary(rhodocov_length_aov)

# Check ANOVA assumptions
plot(rhodocov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_rhodo1) # But, variances are actually NOT homogeneous

plot(rhodocov_length_aov, 2) # Normal distribution of the data NOT assumed
rhodocov_aov_residuals <- residuals(object = rhodocov_length_aov)
shapiro.test(x = rhodocov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
rhodocov_length_glm <- glm(cover ~ length_range, data = cover_rhodo1,
                           family = "Gamma")
summary(rhodocov_length_glm)
Anova(rhodocov_length_glm)
AIC(rhodocov_length_glm) # AIC(glm) = 1247.73

# Check models
DHARMa::simulateResiduals(fittedModel = rhodocov_length_glm, plot = T) # A very good fit !
gratia::appraise(rhodocov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Rhodophyta ####
summary(glht(rhodocov_length_glm, linfct = mcp(length_range = "Tukey")))





# 4.3.2. Percent cover x (Infestation x Degradation x Position) ####
cover_rhodo <- cover_group_length[cover_group_length$higher.group %in% "Rhodophyta",]

cover_rhodo2 <- cover_rhodo %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_rhodo2)

hist(cover_rhodo2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
rhodocov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_rhodo2)
summary(rhodocov_inf_aov)

# Check ANOVA assumptions
plot(rhodocov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_rhodo2) # But, variances are actually NOT homogeneous

plot(rhodocov_inf_aov, 2) # Normal distribution of the data NOT assumed
rhodocov_inf_aov_residuals <- residuals(object = rhodocov_inf_aov)
shapiro.test(x = rhodocov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
rhodocov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_rhodo2,
                           family = "Gamma"(link = "log"))
summary(rhodocov_inf_glm)
Anova(rhodocov_inf_glm)
AIC(rhodocov_inf_glm) # AIC(glm) = 1791.977

# Check models
DHARMa::simulateResiduals(fittedModel = rhodocov_inf_glm, plot = T) # Not the best fit.
gratia::appraise(rhodocov_inf_glm, method = "simulate")

# Model #2: Additive GLM (Gamma) ####
rhodocov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_rhodo2,
                        family = "Gamma")
summary(rhodocov_inf_glm1)
Anova(rhodocov_inf_glm1)
AIC(rhodocov_inf_glm1) # AIC(glm) = 1781.962

# Check models
DHARMa::simulateResiduals(fittedModel = rhodocov_inf_glm1, plot = T) # Very good fit !
gratia::appraise(rhodocov_inf_glm1, method = "simulate")

# Pairwise comparisons on abundance in Rhodophyta ####
summary(glht(rhodocov_inf_glm1, linfct = mcp(inf.percent = "Tukey")))
summary(glht(rhodocov_inf_glm1, linfct = mcp(position = "Tukey")))





# 4.4. Sedentaria ==============================================================================
cover_seden <- cover_group_length[cover_group_length$higher.group %in% "Sedentaria",]

cover_seden1 <- cover_seden %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_seden1)

hist(cover_seden1$cover) # Gamma distribution

# 4.4.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
sedencov_length_aov <- aov(cover ~ length_range, data = cover_seden1)
summary(sedencov_length_aov)

# Check ANOVA assumptions
plot(sedencov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_seden1) # But, variances are actually homogeneous

plot(sedencov_length_aov, 2) # Normal distribution of the data NOT assumed
sedencov_aov_residuals <- residuals(object = sedencov_length_aov)
shapiro.test(x = sedencov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
sedencov_length_glm <- glm(cover ~ length_range, data = cover_seden1,
                           family = "Gamma")
summary(sedencov_length_glm)
Anova(sedencov_length_glm)
AIC(sedencov_length_glm) # AIC(glm) = 826.3935

# Check models
DHARMa::simulateResiduals(fittedModel = sedencov_length_glm, plot = T) # Not a bad fit.
gratia::appraise(sedencov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Sedentaria ####
summary(glht(sedencov_length_glm, linfct = mcp(length_range = "Tukey")))





# 4.4.2. Percent cover x (Infestation x Degradation x Position) ####
cover_seden <- cover_group_length[cover_group_length$higher.group %in% "Sedentaria",]

cover_seden2 <- cover_seden %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_seden2)

hist(cover_seden2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
sedencov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_seden2)
summary(sedencov_inf_aov)

# Check ANOVA assumptions
plot(sedencov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_seden2) # But, variances are actually homogeneous

plot(sedencov_inf_aov, 2) # Normal distribution of the data NOT assumed
sedencov_inf_aov_residuals <- residuals(object = sedencov_inf_aov)
shapiro.test(x = sedencov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
sedencov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_seden2,
                           family = "Gamma"(link = "log"))
summary(sedencov_inf_glm)
Anova(sedencov_inf_glm)
AIC(sedencov_inf_glm) # AIC(glm) = 1007.822

# Check models
DHARMa::simulateResiduals(fittedModel = sedencov_inf_glm, plot = T) # Not a good fit.
gratia::appraise(sedencov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
sedencov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.percent:position, 
                        data = cover_seden2,
                        family = "Gamma")
summary(sedencov_inf_glm1)
Anova(sedencov_inf_glm1)
AIC(sedencov_inf_glm1) # AIC(glm) = 987.3691

# Check models
DHARMa::simulateResiduals(fittedModel = sedencov_inf_glm1, plot = T) # Much better.
gratia::appraise(sedencov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) ####
sedencov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_seden2,
                         family = "Gamma")
summary(sedencov_inf_glm2)
Anova(sedencov_inf_glm2)
AIC(sedencov_inf_glm2) # AIC(glm) = 1003.066

# Check models
DHARMa::simulateResiduals(fittedModel = sedencov_inf_glm2, plot = T) # Not a good fit.
gratia::appraise(sedencov_inf_glm2, method = "simulate")

# Pairwise comparisons on abundance in Sedentaria ####
emmeans(sedencov_inf_glm1, pairwise ~ position | inf.percent)





###############################################################################################
# Section 5: Percent cover of epibiotic species -----------------------------------------------
###############################################################################################

# Summarize percent cover of significant epibiotic species
cover_long <- epibio[,c(1:15, 18,22)]
View(cover_long)

cover_group_sp <- cover_long[cover_long$epibiont %in% c("Spirobranchus kraussii", "Ralfsia verrucosa",
                                                        "Corallinales", "Tetraclita serrata",
                                                        "Hildenbrandia lecannellieri", "Juvenile barnacle",
                                                        "Spirorbis spp.", "Chthamalus dentatus"),]
View(cover_group_sp)





# 5.1. Tetraclita serrata (Cirripedia) ==================================================
cover_tester <- cover_group_sp[cover_group_sp$epibiont %in% "Tetraclita serrata",]

cover_tester1 <- cover_tester %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_tester1)

hist(cover_tester1$cover) # Gamma distribution

# 5.1.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
testercov_length_aov <- aov(cover ~ length_range, data = cover_tester1)
summary(testercov_length_aov)

# Check ANOVA assumptions
plot(testercov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_tester1) # But, variances are actually homogeneous

plot(testercov_length_aov, 2) # Normal distribution of the data NOT assumed
testercov_aov_residuals <- residuals(object = testercov_length_aov)
shapiro.test(x = testercov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
tetsercov_length_glm <- glm(cover ~ length_range, data = cover_tester1,
                           family = "Gamma")
summary(tetsercov_length_glm)
Anova(tetsercov_length_glm)
AIC(tetsercov_length_glm) # AIC(glm) = 927.7784

# Check models
DHARMa::simulateResiduals(fittedModel = tetsercov_length_glm, plot = T) # A very good fit !
gratia::appraise(tetsercov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Tetraclita serrata ####
summary(glht(tetsercov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.1.2. Percent cover x (Infestation x Degradation x Position) ####
cover_tester <- cover_group_sp[cover_group_sp$epibiont %in% "Tetraclita serrata",]

cover_tester2 <- cover_tester %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_tester2)

hist(cover_tester2$cover) # Gamma distribution


# Model #1: Regular ANOVA ####
testercov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_tester2)
summary(testercov_inf_aov)

# Check ANOVA assumptions
plot(testercov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_tester2) # But, variances are actually homogeneous

plot(testercov_inf_aov, 2) # Normal distribution of the data NOT assumed
testercov_inf_aov_residuals <- residuals(object = testercov_inf_aov)
shapiro.test(x = testercov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
tetsercov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_tester2,
                            family = "Gamma")
summary(tetsercov_inf_glm)
Anova(tetsercov_inf_glm) # Significant: inf.lvl:inf.percent and inf.percent:position
AIC(tetsercov_inf_glm) # AIC(glm) = 1195.346

# Check models
DHARMa::simulateResiduals(fittedModel = tetsercov_inf_glm, plot = T) # A very good fit !
gratia::appraise(tetsercov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interaction ####
tetsercov_inf_glm1 <- glm(cover ~ inf.lvl * inf.percent + position + inf.percent:position, data = cover_tester2,
                         family = "Gamma")
summary(tetsercov_inf_glm1)
Anova(tetsercov_inf_glm1)
AIC(tetsercov_inf_glm1) # AIC(glm) = 1186.37

# Check models
DHARMa::simulateResiduals(fittedModel = tetsercov_inf_glm1, plot = T) # A very good fit !
gratia::appraise(tetsercov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) w/ some interaction ####
tetsercov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_tester2,
                          family = "Gamma")
summary(tetsercov_inf_glm2)
Anova(tetsercov_inf_glm2)
AIC(tetsercov_inf_glm2) # AIC(glm) = 1194.284

# Check models
DHARMa::simulateResiduals(fittedModel = tetsercov_inf_glm2, plot = T) # A very good fit !
gratia::appraise(tetsercov_inf_glm1, method = "simulate")

# Pairwise comparisons on abundance in Tetraclita serrata ####
summary(glht(tetsercov_inf_glm2, linfct = mcp(inf.percent = "Tukey")))
summary(glht(tetsercov_inf_glm2, linfct = mcp(inf.lvl = "Tukey")))
summary(glht(tetsercov_inf_glm2, linfct = mcp(position = "Tukey")))

emmeans(tetsercov_inf_glm1, pairwise ~ inf.lvl | inf.percent)
emmeans(tetsercov_inf_glm1, pairwise ~ position | inf.percent)





# 5.2. Chthamalus dentatus (Cirripedia) ==================================================
cover_chtden <- cover_group_sp[cover_group_sp$epibiont %in% "Chthamalus dentatus",]

cover_chtden1 <- cover_chtden %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_chtden1)

hist(cover_chtden1$cover) # Gamma distribution

# 5.2.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
chtdencov_length_aov <- aov(cover ~ length_range, data = cover_chtden1)
summary(chtdencov_length_aov)

# Check ANOVA assumptions
plot(chtdencov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_chtden1) # But, variances are actually homogeneous

plot(chtdencov_length_aov, 2) # Normal distribution of the data NOT assumed
chtdencov_aov_residuals <- residuals(object = chtdencov_length_aov)
shapiro.test(x = chtdencov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
chtdencov_length_glm <- glm(cover ~ length_range, data = cover_chtden1,
                            family = "Gamma")
summary(chtdencov_length_glm)
Anova(chtdencov_length_glm)
AIC(chtdencov_length_glm) # AIC(glm) = 88.3736

# Check models
DHARMa::simulateResiduals(fittedModel = chtdencov_length_glm, plot = T) # A very good fit !
gratia::appraise(chtdencov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Chthamalus dentatus ####
summary(glht(chtdencov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.3. Juvenile barnacle (Cirripedia) ==================================================
cover_juvbar <- cover_group_sp[cover_group_sp$epibiont %in% "Juvenile barnacle",]

cover_juvbar1 <- cover_juvbar %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_juvbar1)

hist(cover_juvbar1$cover) # Gamma distribution

# 5.3.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
juvbarcov_length_aov <- aov(cover ~ length_range, data = cover_juvbar1)
summary(juvbarcov_length_aov)

# Check ANOVA assumptions
plot(juvbarcov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_juvbar1) # But, variances are actually NOT homogeneous

plot(juvbarcov_length_aov, 2) # Normal distribution of the data NOT assumed
juvbarcov_aov_residuals <- residuals(object = juvbarcov_length_aov)
shapiro.test(x = juvbarcov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
juvbarcov_length_glm <- glm(cover ~ length_range, data = cover_juvbar1,
                            family = "Gamma")
summary(juvbarcov_length_glm)
Anova(juvbarcov_length_glm)
AIC(juvbarcov_length_glm) # AIC(glm) = -39.14415

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarcov_length_glm, plot = T) # Not such a bad fit. 
gratia::appraise(juvbarcov_length_glm, method = "simulate")

# Pairwise comparisons on percent cover in Chthamalus dentatus ####
summary(glht(juvbarcov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.3.2. Percent cover x (Infestation x Degradation x Position) ####
cover_juvbar <- cover_group_sp[cover_group_sp$epibiont %in% "Juvenile barnacle",]

cover_juvbar2 <- cover_juvbar %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_juvbar2)

hist(cover_juvbar2$cover) # Gamma distribution, skewed to lower values.

# Model #1: Regular ANOVA ####
juvbarcov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_juvbar2)
summary(juvbarcov_inf_aov)

# Check ANOVA assumptions
plot(juvbarcov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_juvbar2) # But, variances are actually NOT homogeneous

plot(juvbarcov_inf_aov, 2) # Normal distribution of the data NOT assumed
juvbarcov_inf_aov_residuals <- residuals(object = juvbarcov_inf_aov)
shapiro.test(x = juvbarcov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) w/ identity link ####
juvbarcov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_juvbar2,
                            family = "Gamma"(link = "identity"))
summary(juvbarcov_inf_glm)
Anova(juvbarcov_inf_glm)
AIC(juvbarcov_inf_glm) # AIC(glm) = -440.1222

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarcov_inf_glm, plot = T) # Not a good fit. 
gratia::appraise(juvbarcov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
juvbarcov_inf_glm1 <- glm(cover ~ inf.lvl * inf.percent + position, data = cover_juvbar2,
                         family = "Gamma"(link = "identity"))
summary(juvbarcov_inf_glm1)
Anova(juvbarcov_inf_glm1)
AIC(juvbarcov_inf_glm1) # AIC(glm) = -440.1395

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarcov_inf_glm1, plot = T) # Not a good fit. 
gratia::appraise(juvbarcov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) ####
juvbarcov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_juvbar2,
                          family = "Gamma"(link = "identity"))
summary(juvbarcov_inf_glm2)
Anova(juvbarcov_inf_glm2)
AIC(juvbarcov_inf_glm2) # AIC(glm) = -334.4069

# Check models
DHARMa::simulateResiduals(fittedModel = juvbarcov_inf_glm2, plot = T) # Not a good fit at all.
gratia::appraise(juvbarcov_inf_glm2, method = "simulate")

# Pairwise comparisons on percent cover of juvenile barnacles ####
emmeans(juvbarcov_inf_glm1, pairwise ~ inf.lvl | inf.percent)








# 5.4. Ralfsia verrucosa (Ochrophyta) ==================================================
cover_ralver <- cover_group_sp[cover_group_sp$epibiont %in% "Ralfsia verrucosa",]

cover_ralver1 <- cover_ralver %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_ralver1)

hist(cover_ralver1$cover) # Gamma distribution

# 5.4.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
ralvercov_length_aov <- aov(cover ~ length_range, data = cover_ralver1)
summary(ralvercov_length_aov)

# Check ANOVA assumptions
plot(ralvercov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_ralver1) # But, variances are actually homogeneous

plot(ralvercov_length_aov, 2) # Normal distribution of the data NOT assumed
ralvercov_aov_residuals <- residuals(object = ralvercov_length_aov)
shapiro.test(x = ralvercov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
ralvercov_length_glm <- glm(cover ~ length_range, data = cover_ralver1,
                            family = "Gamma")
summary(ralvercov_length_glm)
Anova(ralvercov_length_glm)
AIC(ralvercov_length_glm) # AIC(glm) = 1106.441

# Check models
DHARMa::simulateResiduals(fittedModel = ralvercov_length_glm, plot = T) # A very good fit !
gratia::appraise(ralvercov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Chthamalus dentatus ####
summary(glht(ralvercov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.4.2. Percent cover x (Infestation x Degradation x Position) ####
cover_ralver <- cover_group_sp[cover_group_sp$epibiont %in% "Ralfsia verrucosa",]

cover_ralver2 <- cover_ralver %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_ralver2)

hist(cover_ralver2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
ralvercov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_ralver2)
summary(ralvercov_inf_aov)

# Check ANOVA assumptions
plot(ralvercov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_ralver2) # But, variances are actually homogeneous

plot(ralvercov_inf_aov, 2) # Normal distribution of the data NOT assumed
ralvercov_inf_aov_residuals <- residuals(object = ralvercov_inf_aov)
shapiro.test(x = ralvercov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
ralvercov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_ralver2,
                            family = "Gamma")
summary(ralvercov_inf_glm)
Anova(ralvercov_inf_glm)
AIC(ralvercov_inf_glm) # AIC(glm) = 1535.19

# Check models
DHARMa::simulateResiduals(fittedModel = ralvercov_inf_glm, plot = T) # A very good fit !
gratia::appraise(ralvercov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
ralvercov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.lvl:position, data = cover_ralver2,
                         family = "Gamma")
summary(ralvercov_inf_glm1)
Anova(ralvercov_inf_glm1)
AIC(ralvercov_inf_glm1) # AIC(glm) = 1519.471

# Check models
DHARMa::simulateResiduals(fittedModel = ralvercov_inf_glm1, plot = T) # A very good fit !
gratia::appraise(ralvercov_inf_glm1, method = "simulate")

# Pairwise comparisons on abundance in Ralfsia verrucosa ####
emmeans(ralvercov_inf_glm1, pairwise ~ position | inf.lvl)





# 5.5. Hildenbrandia lecannellieri (Rhodophyta) ==================================================
cover_hillec <- cover_group_sp[cover_group_sp$epibiont %in% "Hildenbrandia lecannellieri",]

cover_hillec1 <- cover_hillec %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_hillec1)

hist(cover_hillec1$cover) # Gamma distribution

# 5.5.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
hilleccov_length_aov <- aov(cover ~ length_range, data = cover_hillec1)
summary(hilleccov_length_aov)

# Check ANOVA assumptions
plot(hilleccov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_hillec1) # But, variances are actually homogeneous

plot(hilleccov_length_aov, 2) # Normal distribution of the data NOT assumed
hilleccov_aov_residuals <- residuals(object = hilleccov_length_aov)
shapiro.test(x = hilleccov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
hilleccov_length_glm <- glm(cover ~ length_range, data = cover_hillec1,
                            family = "Gamma")
summary(hilleccov_length_glm)
Anova(hilleccov_length_glm)
AIC(hilleccov_length_glm) # AIC(glm) = 533.8952

# Check models
DHARMa::simulateResiduals(fittedModel = hilleccov_length_glm, plot = T) # A very good fit !
gratia::appraise(ralvercov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Hildenbrandia lecannellieri ####
summary(glht(hilleccov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.5.2. Percent cover x (Infestation x Degradation x Position) ####
cover_hillec <- cover_group_sp[cover_group_sp$epibiont %in% "Hildenbrandia lecannellieri",]

cover_hillec2 <- cover_hillec %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_hillec2)

hist(cover_hillec2$cover) # Gamma distributio

# Model #1: Regular ANOVA ####
hilleccov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_hillec2)
summary(hilleccov_inf_aov)

# Check ANOVA assumptions
plot(hilleccov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_hillec2) # But, variances are actually homogeneous

plot(hilleccov_inf_aov, 2) # Normal distribution of the data NOT assumed
hilleccov_inf_aov_residuals <- residuals(object = hilleccov_inf_aov)
shapiro.test(x = hilleccov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
hilleccov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_hillec2,
                            family = "Gamma"(link = "log"))
summary(hilleccov_inf_glm)
Anova(hilleccov_inf_glm)
AIC(hilleccov_inf_glm) # AIC(glm) = 701.9109

# Check models
DHARMa::simulateResiduals(fittedModel = hilleccov_inf_glm, plot = T) # A very good fit !
gratia::appraise(hilleccov_inf_glm, method = "simulate")

# Model #2: GLM (Gamma) w/ length range ####
hilleccov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.percent:position, data = cover_hillec2,
                         family = "Gamma")
summary(hilleccov_inf_glm1)
Anova(hilleccov_inf_glm1)
AIC(hilleccov_inf_glm1) # AIC(glm) = 699.5997

# Check models
DHARMa::simulateResiduals(fittedModel = hilleccov_inf_glm1, plot = T) # A very good fit !
gratia::appraise(hilleccov_inf_glm1, method = "simulate")

# Pairwise comparisons on abundance in Hildenbrandia lecannellieri ####
emmeans(hilleccov_inf_glm1, pairwise ~ position | inf.percent)





# 5.6. Corallinales (Rhodophyta) ==================================================
cover_coral <- cover_group_sp[cover_group_sp$epibiont %in% "Corallinales",]

cover_coral1 <- cover_coral %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_coral1)

hist(cover_coral1$cover) # Gamma distribution

# 5.6.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
coralcov_length_aov <- aov(cover ~ length_range, data = cover_coral1)
summary(coralcov_length_aov)

# Check ANOVA assumptions
plot(coralcov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_coral1) # But, variances are actually homogeneous

plot(coralcov_length_aov, 2) # Normal distribution of the data NOT assumed
coralcov_aov_residuals <- residuals(object = coralcov_length_aov)
shapiro.test(x = coralcov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
coralcov_length_glm <- glm(cover ~ length_range, data = cover_coral1,
                            family = "Gamma")
summary(coralcov_length_glm)
Anova(coralcov_length_glm)
AIC(coralcov_length_glm) # AIC(glm) = 618.3854

# Check models
DHARMa::simulateResiduals(fittedModel = coralcov_length_glm, plot = T) # A very good fit !
gratia::appraise(coralcov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Corallinales ####
summary(glht(coralcov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.6.2. Percent cover x (Infestation x Degradation x Position) ####
cover_coral <- cover_group_sp[cover_group_sp$epibiont %in% "Corallinales",]

cover_coral2 <- cover_coral %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_coral2)

hist(cover_coral2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
coralcov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_coral2)
summary(coralcov_inf_aov)

# Check ANOVA assumptions
plot(coralcov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_coral2) # But, variances are actually homogeneous

plot(coralcov_inf_aov, 2) # Normal distribution of the data NOT assumed
coralcov_inf_aov_residuals <- residuals(object = coralcov_inf_aov)
shapiro.test(x = coralcov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
coralcov_inf_glm <- glm(cover ~ inf.lvl + inf.percent * position, data = cover_coral2,
                           family = "Gamma"(link = "log"))
summary(coralcov_inf_glm)
Anova(coralcov_inf_glm)
AIC(coralcov_inf_glm) # AIC(glm) = 781.1399

# Check models
DHARMa::simulateResiduals(fittedModel = coralcov_inf_glm, plot = T) # A very good fit !
gratia::appraise(coralcov_inf_glm, method = "simulate")

# Pairwise comparisons on abundance in Corallinales ####
emmeans(coralcov_inf_glm, pairwise ~ inf.lvl)
emmeans(coralcov_inf_glm, pairwise ~ position | inf.percent)




# 5.7. Spirobranchus kraussii (Sedentaria) ==================================================
cover_spikra <- cover_group_sp[cover_group_sp$epibiont %in% "Spirobranchus kraussii",]

cover_spikra1 <- cover_spikra %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_spikra1)

hist(cover_spikra1$cover) # Gamma distribution

# 5.7.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
spikracov_length_aov <- aov(cover ~ length_range, data = cover_spikra1)
summary(spikracov_length_aov)

# Check ANOVA assumptions
plot(spikracov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_spikra1) # But, variances are actually homogeneous

plot(spikracov_length_aov, 2) # Normal distribution of the data NOT assumed
spikracov_aov_residuals <- residuals(object = spikracov_length_aov)
shapiro.test(x = spikracov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
spikracov_length_glm <- glm(cover ~ length_range, data = cover_spikra1,
                           family = "Gamma")
summary(spikracov_length_glm)
Anova(spikracov_length_glm)
AIC(spikracov_length_glm) # AIC(glm) = 294.8484

# Check models
DHARMa::simulateResiduals(fittedModel = spikracov_length_glm, plot = T) # A very good fit !
gratia::appraise(spikracov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Spirobranchus kraussii ####
summary(glht(spikracov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.7.2. Percent cover x (Infestation x Degradation x Position) ####
cover_spikra <- cover_group_sp[cover_group_sp$epibiont %in% "Spirobranchus kraussii",]

cover_spikra2 <- cover_spikra %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_spikra2)

hist(cover_spikra2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
spikracov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_spikra2)
summary(spikracov_inf_aov)

# Check ANOVA assumptions
plot(spikracov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_spikra2) # But, variances are actually homogeneous

plot(spikracov_inf_aov, 2) # Normal distribution of the data NOT assumed
spikracov_inf_aov_residuals <- residuals(object = spikracov_inf_aov)
shapiro.test(x = spikracov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
spikracov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_spikra2,
                            family = "Gamma"(link = "log"))
summary(spikracov_inf_glm)
Anova(spikracov_inf_glm)
AIC(spikracov_inf_glm) # AIC(glm) = 314.0788

# Check models
DHARMa::simulateResiduals(fittedModel = spikracov_inf_glm, plot = T) # A very good fit !
gratia::appraise(spikracov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
spikracov_inf_glm1 <- glm(cover ~ inf.lvl + inf.percent + position + inf.lvl:inf.percent:position, data = cover_spikra2,
                         family = "Gamma"(link = "log"))
summary(spikracov_inf_glm1)
Anova(spikracov_inf_glm1)
AIC(spikracov_inf_glm1) # AIC(glm) = 314.0788

# Check models
DHARMa::simulateResiduals(fittedModel = spikracov_inf_glm1, plot = T) # A very good fit !
gratia::appraise(spikracov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) ####
spikracov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_spikra2,
                          family = "Gamma"(link = "log"))
summary(spikracov_inf_glm2)
Anova(spikracov_inf_glm2)
AIC(spikracov_inf_glm2) # AIC(glm) = 307.5953

# Check models
DHARMa::simulateResiduals(fittedModel = spikracov_inf_glm2, plot = T) # A very good fit !
gratia::appraise(spikracov_inf_glm2, method = "simulate")

# Pairwise comparisons on abundance in Spirobranchus kraussii ####
summary(glht(spikracov_inf_glm2, linfct = mcp(inf.percent = "Tukey")))
summary(glht(spikracov_inf_glm2, linfct = mcp(position = "Tukey")))





# 5.8. Spirorbis sp. (Sedentaria) ==================================================
cover_spiror <- cover_group_sp[cover_group_sp$epibiont %in% "Spirorbis spp.",]

cover_spiror1 <- cover_spiror %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_spiror1)

hist(cover_spiror1$cover) # Gamma distribution

# 5.7.1. Percent cover x Shell length ####

# Model #1: Regular ANOVA ####
spirorcov_length_aov <- aov(cover ~ length_range, data = cover_spiror1)
summary(spirorcov_length_aov)

# Check ANOVA assumptions
plot(spirorcov_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ length_range, data = cover_spiror1) # But, variances are actually homogeneous

plot(spirorcov_length_aov, 2) # Normal distribution of the data NOT assumed
spirorcov_aov_residuals <- residuals(object = spirorcov_length_aov)
shapiro.test(x = spirorcov_aov_residuals)  # Residuals not normally distributed

# Model #2: GLM (Gamma) w/ length range ####
spirorcov_length_glm <- glm(cover ~ length_range, data = cover_spiror1,
                            family = "Gamma")
summary(spirorcov_length_glm)
Anova(spirorcov_length_glm)
AIC(spirorcov_length_glm) # AIC(glm) = 96.41255

# Check models
DHARMa::simulateResiduals(fittedModel = spirorcov_length_glm, plot = T) # A good enough fit.
gratia::appraise(spirorcov_length_glm, method = "simulate")

# Pairwise comparisons on abundance in Spirorbis spp. ####
summary(glht(spirorcov_length_glm, linfct = mcp(length_range = "Tukey")))





# 5.7.2. Percent cover x (Infestation x Degradation x Position) ####
cover_spiror <- cover_group_sp[cover_group_sp$epibiont %in% "Spirorbis spp.",]

cover_spiror2 <- cover_spiror %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_spiror2)

hist(cover_spiror2$cover) # Gamma distribution

# Model #1: Regular ANOVA ####
spirorcov_inf_aov <- aov(cover ~ inf.lvl * inf.percent * position, data = cover_spiror2)
summary(spirorcov_inf_aov)

# Check ANOVA assumptions
plot(spirorcov_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(cover ~ inf.lvl * inf.percent * position, data = cover_spiror2) # But, variances are actually homogeneous

plot(spirorcov_inf_aov, 2) # Normal distribution of the data NOT assumed
spirorcov_inf_aov_residuals <- residuals(object = spirorcov_inf_aov)
shapiro.test(x = spirorcov_inf_aov_residuals)  # Residuals not normally distributed

# Model #2: Saturated GLM (Gamma) ####
spirorcov_inf_glm <- glm(cover ~ inf.lvl * inf.percent * position, data = cover_spiror2,
                            family = "Gamma")
summary(spirorcov_inf_glm)
Anova(spirorcov_inf_glm)
AIC(spirorcov_inf_glm) # AIC(glm) = -59.20631

# Check models
DHARMa::simulateResiduals(fittedModel = spirorcov_inf_glm, plot = T) # Not a good fit.
gratia::appraise(spirorcov_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
spirorcov_inf_glm1 <- glm(cover ~ inf.lvl * inf.percent + position, data = cover_spiror2,
                         family = "Gamma")
summary(spirorcov_inf_glm1)
Anova(spirorcov_inf_glm1)
AIC(spirorcov_inf_glm1) # AIC(glm) = -72.35347

# Check models
DHARMa::simulateResiduals(fittedModel = spirorcov_inf_glm1, plot = T) # Better fit.
gratia::appraise(spirorcov_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) ####
spirorcov_inf_glm2 <- glm(cover ~ inf.lvl + inf.percent + position, data = cover_spiror2,
                          family = "Gamma")
summary(spirorcov_inf_glm2)
Anova(spirorcov_inf_glm2)
AIC(spirorcov_inf_glm2) # AIC(glm) = -81.44618

# Check models
DHARMa::simulateResiduals(fittedModel = spirorcov_inf_glm2, plot = T) # Not so good.
gratia::appraise(spirorcov_inf_glm2, method = "simulate")

# Model #5: Additive GLM (Gamma) w/ some interactions ####
spirorcov_inf_glm3 <- glm(cover ~ inf.lvl + inf.percent + position + inf.lvl:position + inf.percent:position, 
                          data = cover_spiror2,
                          family = "Gamma")
summary(spirorcov_inf_glm3)
Anova(spirorcov_inf_glm3)
AIC(spirorcov_inf_glm3) # AIC(glm) = -74.49347

# Check models
DHARMa::simulateResiduals(fittedModel = spirorcov_inf_glm3, plot = T) # Not so good.
gratia::appraise(spirorcov_inf_glm3, method = "simulate")

# Pairwise comparisons on abundance in Spirorbis spp. ####
summary(glht(spirorcov_length_glm, linfct = mcp(length_range = "Tukey")))





###############################################################################################
# Section 6: Biomass of higher groups -------------------------------------------------------
###############################################################################################

# Summarize biomass of significant epibiotic higher groups
View(epibio)

biom_long <- epibio[,c(1:15,19,22)]
View(biom_long)

biom_long_length <- biom_long[biom_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Rhodophyta",
                                                            "Bryozoa", "Bivalvia", "Sedentaria"),]
View(biom_long_length)





# 6.1. Cirripedia ====

biom_cirri <- biom_long_length[biom_long_length$higher.group %in% "Cirripedia",]
View(biom_cirri)

# 6.1.1. Cirripedia - Biomass vs Shell length ====

biom_cirri1 <- biom_cirri %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_cirri1)

hist(biom_cirri1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
cirribiom_length_aov <- aov(biomass ~ length_range, data = biom_cirri1)
summary(cirribiom_length_aov)

# Check ANOVA assumptions
plot(cirribiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_cirri1) # But, variances are actually homogeneous

plot(cirribiom_length_aov, 2) # Normal distribution of the data NOT assumed
cirribiom_length_aov_residuals <- residuals(object = cirribiom_length_aov)
shapiro.test(x = cirribiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
cirribiom_length_glm <- glm(biomass ~ length_range, data = biom_cirri1,
                         family = "Gamma")
summary(cirribiom_length_glm)
Anova(cirribiom_length_glm)
AIC(cirribiom_length_glm) # AIC(glm) = 1149.727

# Check models
DHARMa::simulateResiduals(fittedModel = cirribiom_length_glm, plot = T) # Not a bad fit.
gratia::appraise(cirribiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Cirripedia ####
summary(glht(cirribiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 6.1.2. Cirripedia - Biomass vs (Infestation x Degradation x Position) ====

biom_cirri2 <- biom_cirri %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_cirri2)

hist(biom_cirri2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
cirribiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_cirri2)
summary(cirribiom_inf_aov)

# Check ANOVA assumptions
plot(cirribiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_cirri2) # But, variances are actually homogeneous

plot(cirribiom_inf_aov, 2) # Normal distribution of the data NOT assumed
cirribiom_inf_aov_residuals <- residuals(object = cirribiom_inf_aov)
shapiro.test(x = cirribiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
cirribiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_cirri2,
                            family = "Gamma")
summary(cirribiom_inf_glm)
Anova(cirribiom_inf_glm)
AIC(cirribiom_inf_glm) # AIC(glm) = 1675.999

# Check models
DHARMa::simulateResiduals(fittedModel = cirribiom_inf_glm, plot = T) # Not such a bad fit.
gratia::appraise(cirribiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
cirribiom_inf_glm1 <- glm(biomass ~ inf.lvl * inf.percent + position, data = biom_cirri2,
                         family = "Gamma")
summary(cirribiom_inf_glm1)
Anova(cirribiom_inf_glm1)
AIC(cirribiom_inf_glm1) # AIC(glm) = 1666.02

# Check models
DHARMa::simulateResiduals(fittedModel = cirribiom_inf_glm1, plot = T) # Not such a bad fit.
gratia::appraise(cirribiom_inf_glm1, method = "simulate")

# Model #4: Interactive GLM (Gamma) w/out position ####
cirribiom_inf_glm2 <- glm(biomass ~ inf.lvl * inf.percent, data = biom_cirri2,
                          family = "Gamma")
summary(cirribiom_inf_glm2)
Anova(cirribiom_inf_glm2)
AIC(cirribiom_inf_glm2) # AIC(glm) = 1664.878

# Check models
DHARMa::simulateResiduals(fittedModel = cirribiom_inf_glm2, plot = T) # Not such a bad fit.
gratia::appraise(cirribiom_inf_glm2, method = "simulate")

# Pairwise comparisons on biomass in Cirripedia ####
emmeans(cirribiom_inf_glm1, pairwise ~ inf.lvl | inf.percent)





# 6.2. Ochrophyta ====

biom_ochro <- biom_long_length[biom_long_length$higher.group %in% "Ochrophyta",]
View(biom_ochro)

# 6.2.1. Ochrophyta - Biomass vs Shell length ====

biom_ochro1 <- biom_ochro %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_ochro1)

hist(biom_ochro1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
ochrobiom_length_aov <- aov(biomass ~ length_range, data = biom_ochro1)
summary(ochrobiom_length_aov)

# Check ANOVA assumptions
plot(ochrobiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_ochro1) # But, variances are actually homogeneous

plot(ochrobiom_length_aov, 2) # Normal distribution of the data NOT assumed
ochrobiom_length_aov_residuals <- residuals(object = ochrobiom_length_aov)
shapiro.test(x = ochrobiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
ochrobiom_length_glm <- glm(biomass ~ length_range, data = biom_ochro1,
                            family = "Gamma")
summary(ochrobiom_length_glm)
Anova(ochrobiom_length_glm)
AIC(ochrobiom_length_glm) # AIC(glm) = 1320.333

# Check models
DHARMa::simulateResiduals(fittedModel = ochrobiom_length_glm, plot = T) # A very good fit !
gratia::appraise(ochrobiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Ochrophyta ####
summary(glht(ochrobiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 6.2.2. Ochrophyta - Biomass vs (Infestation x Degradation x Position) ====

biom_ochro2 <- biom_ochro %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_ochro2)

hist(biom_ochro2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
ochrobiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_ochro2)
summary(ochrobiom_inf_aov)

# Check ANOVA assumptions
plot(ochrobiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_ochro2) # But, variances are actually homogeneous

plot(ochrobiom_inf_aov, 2) # Normal distribution of the data NOT assumed
ochrobiom_inf_aov_residuals <- residuals(object = ochrobiom_inf_aov)
shapiro.test(x = ochrobiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
ochrobiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_ochro2,
                            family = "Gamma")
summary(ochrobiom_inf_glm)
Anova(ochrobiom_inf_glm)
AIC(ochrobiom_inf_glm) # AIC(glm) = 1905.265

# Check models
DHARMa::simulateResiduals(fittedModel = ochrobiom_inf_glm, plot = T) # A very good fit !
gratia::appraise(ochrobiom_inf_glm, method = "simulate")

# Model #2: Additive GLM (Gamma) ####
ochrobiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_ochro2,
                         family = "Gamma")
summary(ochrobiom_inf_glm1)
Anova(ochrobiom_inf_glm1)
AIC(ochrobiom_inf_glm1) # AIC(glm) = 1899.233

# Check models
DHARMa::simulateResiduals(fittedModel = ochrobiom_inf_glm1, plot = T) # A very good fit !
gratia::appraise(ochrobiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in Ochrophyta ####
summary(glht(ochrobiom_inf_glm1, linfct = mcp(position = "Tukey")))




# 6.3. Rhodophyta ====

biom_rhodo <- biom_long_length[biom_long_length$higher.group %in% "Rhodophyta",]
View(biom_rhodo)

# 6.3.1. Rhodophyta - Biomass vs Shell length ====

biom_rhodo1 <- biom_rhodo %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_rhodo1)

hist(biom_rhodo1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
rhodobiom_length_aov <- aov(biomass ~ length_range, data = biom_rhodo1)
summary(rhodobiom_length_aov)

# Check ANOVA assumptions
plot(rhodobiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_rhodo1) # But, variances are actually homogeneous

plot(rhodobiom_length_aov, 2) # Normal distribution of the data NOT assumed
rhodobiom_length_aov_residuals <- residuals(object = rhodobiom_length_aov)
shapiro.test(x = rhodobiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
rhodobiom_length_glm <- glm(biomass ~ length_range, data = biom_rhodo1,
                            family = "Gamma")
summary(rhodobiom_length_glm)
Anova(rhodobiom_length_glm)
AIC(rhodobiom_length_glm) # AIC(glm) = 1176.815

# Check models
DHARMa::simulateResiduals(fittedModel = rhodobiom_length_glm, plot = T) # Not a bad fit.
gratia::appraise(rhodobiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Ochrophyta ####
summary(glht(rhodobiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 6.3.2. Rhodophyta - Biomass vs (Infestation x Degradation x Position) ====

biom_rhodo2 <- biom_rhodo %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_rhodo2)

hist(biom_rhodo2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
rhodobiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_rhodo2)
summary(rhodobiom_inf_aov)

# Check ANOVA assumptions
plot(rhodobiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_rhodo2) # But, variances are actually homogeneous

plot(rhodobiom_inf_aov, 2) # Normal distribution of the data NOT assumed
rhodobiom_inf_aov_residuals <- residuals(object = rhodobiom_inf_aov)
shapiro.test(x = rhodobiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
rhodobiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_rhodo2,
                            family = "Gamma"(link = "log"))
summary(rhodobiom_inf_glm)
Anova(rhodobiom_inf_glm)
AIC(rhodobiom_inf_glm) # AIC(glm) = 1611.118

# Check models
DHARMa::simulateResiduals(fittedModel = rhodobiom_inf_glm, plot = T) # Not a good fit.
gratia::appraise(rhodobiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
rhodobiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_rhodo2,
                         family = "Gamma"(link = "identity"))
summary(rhodobiom_inf_glm1)
Anova(rhodobiom_inf_glm1)
AIC(rhodobiom_inf_glm1) # AIC(glm) = 1622.127

# Check models
DHARMa::simulateResiduals(fittedModel = rhodobiom_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(rhodobiom_inf_glm1, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
rhodobiom_inf_glm2 <- glm(biomass ~ inf.lvl + inf.percent + inf.lvl:position + inf.percent:position, 
                          data = biom_rhodo2,
                          family = "Gamma"(link = "log"))
summary(rhodobiom_inf_glm2)
Anova(rhodobiom_inf_glm2)
AIC(rhodobiom_inf_glm2) # AIC(glm) = 1619.004

# Check models
DHARMa::simulateResiduals(fittedModel = rhodobiom_inf_glm2, plot = T) # Not a good fit.
gratia::appraise(rhodobiom_inf_glm2, method = "simulate")

# Pairwise comparisons on biomass in Ochrophyta ####
summary(glht(rhodobiom_inf_glm1, linfct = mcp(inf.lvl = "Tukey")))





# 6.4. Bryozoa ====

biom_bryo <- biom_long_length[biom_long_length$higher.group %in% "Bryozoa",]
View(biom_bryo)

# 6.4.1. Bryozoa - Biomass vs Shell length ====

biom_bryo1 <- biom_bryo %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bryo1)

hist(biom_bryo1$biomass)

# Model #1: Regular ANOVA ####
bryobiom_length_aov <- aov(biomass ~ length_range, data = biom_bryo1)
summary(bryobiom_length_aov)

# Check ANOVA assumptions
plot(bryobiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_bryo1) # But, variances are actually homogeneous

plot(bryobiom_length_aov, 2) # Normal distribution of the data NOT assumed
bryobiom_length_aov_residuals <- residuals(object = bryobiom_length_aov)
shapiro.test(x = bryobiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
bryobiom_length_glm <- glm(biomass ~ length_range, data = biom_bryo1,
                            family = "Gamma")
summary(bryobiom_length_glm)
Anova(bryobiom_length_glm)
AIC(bryobiom_length_glm) # AIC(glm) = 739.1365

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_length_glm, plot = T) # Good fit !
gratia::appraise(bryobiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Bryozoa ####
summary(glht(bryobiom_length_glm, linfct = mcp(length_range = "Tukey")))






# 6.4.2. Bryozoa - Biomass vs (Infestation x Degradation x Position) ====

biom_bryo2 <- biom_bryo %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bryo2)

hist(biom_bryo2$biomass)

# Model #1: Regular ANOVA ####
bryobiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_bryo2)
summary(bryobiom_inf_aov)

# Check ANOVA assumptions
plot(bryobiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_bryo2) # But, variances are actually homogeneous

plot(bryobiom_inf_aov, 2) # Normal distribution of the data NOT assumed
bryobiom_inf_aov_residuals <- residuals(object = bryobiom_inf_aov)
shapiro.test(x = bryobiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
bryobiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_bryo2,
                           family = "Gamma")
summary(bryobiom_inf_glm)
Anova(bryobiom_inf_glm)
AIC(bryobiom_inf_glm) # AIC(glm) = 855.2084

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
bryobiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_bryo2,
                        family = "Gamma")
summary(bryobiom_inf_glm1)
Anova(bryobiom_inf_glm1)
AIC(bryobiom_inf_glm1) # AIC(glm) = 845.662

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) w/ some interactions ####
bryobiom_inf_glm2 <- glm(biomass ~ inf.lvl * inf.percent + position, data = biom_bryo2,
                         family = "Gamma")
summary(bryobiom_inf_glm2)
Anova(bryobiom_inf_glm2)
AIC(bryobiom_inf_glm2) # AIC(glm) = 851.2689

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm2, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm2, method = "simulate")

# Model #5: Additive GLM (Gamma) w/ some interactions ####
bryobiom_inf_glm3 <- glm(biomass ~ inf.lvl * inf.percent + position + inf.lvl:position, data = biom_bryo2,
                         family = "Gamma")
summary(bryobiom_inf_glm3)
Anova(bryobiom_inf_glm3)
AIC(bryobiom_inf_glm3) # AIC(glm) = 849.3644

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm3, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm3, method = "simulate")

# Model #6: Additive GLM (Gamma) w/ some interactions ####
bryobiom_inf_glm4 <- glm(biomass ~ inf.lvl + inf.percent + position + inf.lvl:position, data = biom_bryo2,
                         family = "Gamma")
summary(bryobiom_inf_glm4)
Anova(bryobiom_inf_glm4)
AIC(bryobiom_inf_glm4) # AIC(glm) = 843.3194

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm4, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm4, method = "simulate")

# Model #7: Additive GLM (Gamma) w/ some interactions ####
bryobiom_inf_glm5 <- glm(biomass ~ inf.lvl + inf.percent + position + inf.lvl:position + inf.percent:position, 
                         data = biom_bryo2,
                         family = "Gamma")
summary(bryobiom_inf_glm5)
Anova(bryobiom_inf_glm5)
AIC(bryobiom_inf_glm5) # AIC(glm) = 846.9642

# Check models
DHARMa::simulateResiduals(fittedModel = bryobiom_inf_glm5, plot = T) # Not a good fit.
gratia::appraise(bryobiom_inf_glm5, method = "simulate")

# Pairwise comparisons on biomass in Bryozoa ####
summary(glht(bryobiom_length_glm, linfct = mcp(length_range = "Tukey")))




# 6.5. Bivalvia ====

biom_bival <- biom_long_length[biom_long_length$higher.group %in% "Bivalvia",]
View(biom_bival)

# 6.5.1. Bivalvia - Biomass vs Shell length ====

biom_bival1 <- biom_bival %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bival1)

hist(biom_bival1$biomass)

# Model #1: Regular ANOVA ####
bivalbiom_length_aov <- aov(biomass ~ length_range, data = biom_bival1)
summary(bivalbiom_length_aov)

# Check ANOVA assumptions
plot(bivalbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_bival1) # But, variances are actually homogeneous

plot(bivalbiom_length_aov, 2) # Normal distribution of the data NOT assumed
bivalbiom_length_aov_residuals <- residuals(object = bivalbiom_length_aov)
shapiro.test(x = bivalbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
bivalbiom_length_glm <- glm(biomass ~ length_range, data = biom_bival1,
                           family = "Gamma")
summary(bivalbiom_length_glm)
Anova(bivalbiom_length_glm)
AIC(bivalbiom_length_glm) # AIC(glm) = 406.9936

# Check models
DHARMa::simulateResiduals(fittedModel = bivalbiom_length_glm, plot = T) # Good fit !
gratia::appraise(bivalbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Bivalvia ####
summary(glht(bivalbiom_length_glm, linfct = mcp(length_range = "Tukey")))






# 6.5.2. Bivalvia - Biomass vs (Infestation x Degradation x Position) ====

biom_bival2 <- biom_bival %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bival2)

hist(biom_bival2$biomass)

# Model #1: Regular ANOVA ####
bivalbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_bival2)
summary(bivalbiom_inf_aov)

# Check ANOVA assumptions
plot(bivalbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_bival2) # But, variances are actually homogeneous

plot(bivalbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
bivalbiom_inf_aov_residuals <- residuals(object = bivalbiom_inf_aov)
shapiro.test(x = bivalbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
bivalbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent + position, data = biom_bival2,
                            family = "Gamma"(link = "log"))
summary(bivalbiom_inf_glm)
Anova(bivalbiom_inf_glm)
AIC(bivalbiom_inf_glm) # AIC(glm) = 444.2247

# Check models
DHARMa::simulateResiduals(fittedModel = bivalbiom_inf_glm, plot = T) # Good fit !
gratia::appraise(bivalbiom_inf_glm, method = "simulate")

# Pairwise comparisons on biomass in Bivalvia ####
summary(glht(bivalbiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 6.6. Sedentaria ====

biom_seden <- biom_long_length[biom_long_length$higher.group %in% "Sedentaria",]
View(biom_seden)

# 6.6.1. Sedentaria - Biomass vs Shell length ====

biom_seden1 <- biom_seden %>%
  dplyr::group_by(quadrat, specimen, higher.group, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_seden1)

hist(biom_seden1$biomass)

# Model #1: Regular ANOVA ####
sedenbiom_length_aov <- aov(biomass ~ length_range, data = biom_seden1)
summary(sedenbiom_length_aov)

# Check ANOVA assumptions
plot(sedenbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_seden1) # But, variances are actually homogeneous

plot(sedenbiom_length_aov, 2) # Normal distribution of the data NOT assumed
sedenbiom_length_aov_residuals <- residuals(object = sedenbiom_length_aov)
shapiro.test(x = sedenbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
sedenbiom_length_glm <- glm(biomass ~ length_range, data = biom_seden1,
                            family = "Gamma")
summary(sedenbiom_length_glm)
Anova(sedenbiom_length_glm)
AIC(sedenbiom_length_glm) # AIC(glm) = 498.1261

# Check models
DHARMa::simulateResiduals(fittedModel = sedenbiom_length_glm, plot = T) # Not a bad fit.
gratia::appraise(sedenbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Sedentaria ####
summary(glht(sedenbiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 6.6.2. Sedentaria - Biomass vs (Infestation x Degradation x Position) ====

biom_seden2 <- biom_seden %>%
  dplyr::group_by(quadrat, specimen, valve, higher.group, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_seden2)

hist(biom_seden2$biomass)

# Model #1: Regular ANOVA ####
sedenbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_seden2)
summary(sedenbiom_inf_aov)

# Check ANOVA assumptions
plot(sedenbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_seden2) # But, variances are actually homogeneous

plot(sedenbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
sedenbiom_inf_aov_residuals <- residuals(object = sedenbiom_inf_aov)
shapiro.test(x = sedenbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
sedenbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_seden2,
                            family = "Gamma")
summary(sedenbiom_inf_glm)
Anova(sedenbiom_inf_glm)
AIC(sedenbiom_inf_glm) # AIC(glm) = 557.637

# Check models
DHARMa::simulateResiduals(fittedModel = sedenbiom_inf_glm, plot = T) # Not a bad fit.
gratia::appraise(sedenbiom_inf_glm, method = "simulate")

# Pairwise comparisons on biomass in Sedentaria ####
summary(glht(sedenbiom_length_glm, linfct = mcp(length_range = "Tukey")))



###############################################################################################
# Section 7: Biomass of epibiotic species -------------------------------------------------------
###############################################################################################

# Summarize biomass of significant epibiotic higher groups
View(epibio)

biom_long <- epibio[,c(1:15,19,22)]
View(biom_long)

biom_long_sp <- biom_long[biom_long$epibiont %in% c("Tetraclita serrata", "Chthamalus dentatus",
                                                    "Ralfsia verrucosa", "Hildenbrandia lecannellieri",
                                                    "Electra verticillata", "Bryozoa", "Perna perna",
                                                    "Spirobranchus kraussii", "Gelidium pristoides",
                                                    "Gunnarea gaimardi"),]
View(biom_long_sp)





# 7.1. Tetraclita serrata (Cirripedia) ====

# 7.1.1. TETSER - Biomass vs Shell length ====
biom_tetser <- biom_long_sp[biom_long_sp$epibiont %in% "Tetraclita serrata",]

biom_tetser1 <- biom_tetser %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_tetser1)

hist(biom_tetser1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
tetserbiom_length_aov <- aov(biomass ~ length_range, data = biom_tetser1)
summary(tetserbiom_length_aov)

# Check ANOVA assumptions
plot(tetserbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_tetser1) # But, variances are actually homogeneous

plot(tetserbiom_length_aov, 2) # Normal distribution of the data NOT assumed
tetserbiom_length_aov_residuals <- residuals(object = tetserbiom_length_aov)
shapiro.test(x = tetserbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
tetserbiom_length_glm <- glm(biomass ~ length_range, data = biom_tetser1,
                            family = "Gamma")
summary(tetserbiom_length_glm)
Anova(tetserbiom_length_glm)
AIC(tetserbiom_length_glm) # AIC(glm) = 1129.339

# Check models
DHARMa::simulateResiduals(fittedModel = tetserbiom_length_glm, plot = T) # A very good fit !
gratia::appraise(tetserbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in TETSER ####
summary(glht(tetserbiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 7.1.2. TETSER - Biomass vs (Infestation x Degradation x Position) ====
biom_tetser <- biom_long_sp[biom_long_sp$epibiont %in% "Tetraclita serrata",]

biom_tetser2 <- biom_tetser %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_tetser2)

hist(biom_tetser2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
tetserbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_tetser2)
summary(tetserbiom_inf_aov)

# Check ANOVA assumptions
plot(tetserbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_tetser2) # But, variances are actually homogeneous

plot(tetserbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
tetserbiom_inf_aov_residuals <- residuals(object = tetserbiom_inf_aov)
shapiro.test(x = tetserbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
tetserbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_tetser2,
                             family = "Gamma")
summary(tetserbiom_inf_glm)
Anova(tetserbiom_inf_glm)
AIC(tetserbiom_inf_glm) # AIC(glm) = 1623.086

# Check models
DHARMa::simulateResiduals(fittedModel = tetserbiom_inf_glm, plot = T) # Not such a bad fit.
gratia::appraise(tetserbiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
tetserbiom_inf_glm1 <- glm(biomass ~ inf.lvl * inf.percent + position, data = biom_tetser2,
                          family = "Gamma")
summary(tetserbiom_inf_glm1)
Anova(tetserbiom_inf_glm1)
AIC(tetserbiom_inf_glm1) # AIC(glm) = 1612.584

# Check models
DHARMa::simulateResiduals(fittedModel = tetserbiom_inf_glm1, plot = T) # Very good fit !
gratia::appraise(tetserbiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in TETSER ####
emmeans(tetserbiom_inf_glm1, pairwise ~ inf.lvl | inf.percent)





# 7.2. Chthamalus dentatus (Cirripedia) ====

# 7.2.1. CHTDEN - Biomass vs Shell length ====
biom_chtden <- biom_long_sp[biom_long_sp$epibiont %in% "Chthamalus dentatus",]

biom_chtden1 <- biom_chtden %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_chtden1) # Not a lot of data. 

hist(biom_chtden1$biomass) # Gamma distribution




# 7.3. Ralfsia verrucosa ====

# 7.3.1. RALVER - Biomass vs Shell length ====
biom_ralver <- biom_long_sp[biom_long_sp$epibiont %in% "Ralfsia verrucosa",]

biom_ralver1 <- biom_ralver %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_ralver1)

hist(biom_ralver1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
ralverbiom_length_aov <- aov(biomass ~ length_range, data = biom_ralver1)
summary(ralverbiom_length_aov)

# Check ANOVA assumptions
plot(ralverbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_ralver1) # But, variances are actually homogeneous

plot(ralverbiom_length_aov, 2) # Normal distribution of the data NOT assumed
ralverbiom_length_aov_residuals <- residuals(object = ralverbiom_length_aov)
shapiro.test(x = ralverbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
ralverbiom_length_glm <- glm(biomass ~ length_range, data = biom_ralver1,
                             family = "Gamma")
summary(ralverbiom_length_glm)
Anova(ralverbiom_length_glm)
AIC(ralverbiom_length_glm) # AIC(glm) = 1330.788

# Check models
DHARMa::simulateResiduals(fittedModel = ralverbiom_length_glm, plot = T) # A very good fit !
gratia::appraise(ralverbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in RALVER ####
summary(glht(ralverbiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 7.3.2. RALVER - Biomass vs (Infestation x Degradation x Position) ====
biom_ralver <- biom_long_sp[biom_long_sp$epibiont %in% "Ralfsia verrucosa",]

biom_ralver2 <- biom_ralver %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_ralver2)

hist(biom_ralver2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
ralverbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_ralver2)
summary(ralverbiom_inf_aov)

# Check ANOVA assumptions
plot(ralverbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_ralver2) # But, variances are actually homogeneous

plot(ralverbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
ralverbiom_inf_aov_residuals <- residuals(object = ralverbiom_inf_aov)
shapiro.test(x = ralverbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
ralverbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_ralver2,
                             family = "Gamma")
summary(ralverbiom_inf_glm)
Anova(ralverbiom_inf_glm)
AIC(ralverbiom_inf_glm) # AIC(glm) = 1917.122

# Check models
DHARMa::simulateResiduals(fittedModel = ralverbiom_inf_glm, plot = T) # A very good fit !
gratia::appraise(ralverbiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
ralverbiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_ralver2,
                          family = "Gamma")
summary(ralverbiom_inf_glm1)
Anova(ralverbiom_inf_glm1)
AIC(ralverbiom_inf_glm1) # AIC(glm) = 1910.195

# Check models
DHARMa::simulateResiduals(fittedModel = ralverbiom_inf_glm1, plot = T) # A very good fit !
gratia::appraise(ralverbiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in RALVER ####
summary(glht(ralverbiom_length_glm, linfct = mcp(length_range = "Tukey")))




# 7.4. Hildenbrandia lecannellieri ====

# 7.4.1. HILLEC - Biomass vs Shell length ====
biom_hillec <- biom_long_sp[biom_long_sp$epibiont %in% "Hildenbrandia lecannellieri",]

biom_hillec1 <- biom_hillec %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_hillec1)

hist(biom_hillec1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
hillecbiom_length_aov <- aov(biomass ~ length_range, data = biom_hillec1)
summary(hillecbiom_length_aov)

# Check ANOVA assumptions
plot(hillecbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_hillec1) # But, variances are actually homogeneous

plot(hillecbiom_length_aov, 2) # Normal distribution of the data NOT assumed
hillecbiom_length_aov_residuals <- residuals(object = hillecbiom_length_aov)
shapiro.test(x = hillecbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
hillecbiom_length_glm <- glm(biomass ~ length_range, data = biom_hillec1,
                             family = "Gamma")
summary(hillecbiom_length_glm)
Anova(hillecbiom_length_glm)
AIC(hillecbiom_length_glm) # AIC(glm) = 697.136

# Check models
DHARMa::simulateResiduals(fittedModel = hillecbiom_length_glm, plot = T) # A very good fit !
gratia::appraise(hillecbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in HILLEC ####
summary(glht(hillecbiom_length_glm, linfct = mcp(length_range = "Tukey")))






# 7.4.2. HILLEC - Biomass vs (Infestation x Degradation x Position) ====
biom_hillec <- biom_long_sp[biom_long_sp$epibiont %in% "Hildenbrandia lecannellieri",]

biom_hillec2 <- biom_hillec %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_hillec2)

hist(biom_hillec2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
hillecbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_hillec2)
summary(hillecbiom_inf_aov)

# Check ANOVA assumptions
plot(hillecbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_hillec2) # But, variances are actually homogeneous

plot(hillecbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
hillecbiom_inf_aov_residuals <- residuals(object = hillecbiom_inf_aov)
shapiro.test(x = hillecbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
hillecbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_hillec2,
                             family = "Gamma"(link = "log"))
summary(hillecbiom_inf_glm)
Anova(hillecbiom_inf_glm)
AIC(hillecbiom_inf_glm) # AIC(glm) = 988.0392

# Check models
DHARMa::simulateResiduals(fittedModel = hillecbiom_inf_glm, plot = T) # A very good fit !
gratia::appraise(hillecbiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
hillecbiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent * position, 
                           data = biom_hillec2,
                          family = "Gamma")
summary(hillecbiom_inf_glm1)
Anova(hillecbiom_inf_glm1)
AIC(hillecbiom_inf_glm1) # AIC(glm) = 980.1791

# Check models
DHARMa::simulateResiduals(fittedModel = hillecbiom_inf_glm1, plot = T) # A very good fit !
gratia::appraise(hillecbiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in HILLEC ####
emmeans(hillecbiom_inf_glm1, pairwise ~ inf.lvl)
emmeans(hillecbiom_inf_glm1, pairwise ~ position | inf.percent)




# 7.5. Gelidium pristoides ====

# 7.5.1. GELPRI - Biomass vs Shell length ====
biom_gelpri <- biom_long_sp[biom_long_sp$epibiont %in% "Gelidium pristoides",]

biom_gelpri1 <- biom_gelpri %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_gelpri1)

hist(biom_gelpri1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
gelpribiom_length_aov <- aov(biomass ~ length_range, data = biom_gelpri1)
summary(gelpribiom_length_aov)

# Check ANOVA assumptions
plot(gelpribiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_gelpri1) # But, variances are actually homogeneous

plot(gelpribiom_length_aov, 2) # Normal distribution of the data NOT assumed
gelpribiom_length_aov_residuals <- residuals(object = gelpribiom_length_aov)
shapiro.test(x = gelpribiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
gelpribiom_length_glm <- glm(biomass ~ length_range, data = biom_gelpri1,
                             family = "Gamma")
summary(gelpribiom_length_glm)
Anova(gelpribiom_length_glm)
AIC(gelpribiom_length_glm) # AIC(glm) = 409.3314

# Check models
DHARMa::simulateResiduals(fittedModel = gelpribiom_length_glm, plot = T) # A very good fit !
gratia::appraise(gelpribiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in HILLEC ####
summary(glht(gelpribiom_length_glm, linfct = mcp(length_range = "Tukey")))






# 7.5.2. GELPRI - Biomass vs (Infestation x Degradation x Position) ====
biom_gelpri <- biom_long_sp[biom_long_sp$epibiont %in% "Gelidium pristoides",]

biom_gelpri2 <- biom_gelpri %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_gelpri2)

hist(biom_gelpri2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
gelpribiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_gelpri2)
summary(gelpribiom_inf_aov)

# Check ANOVA assumptions
plot(gelpribiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_gelpri2) # Variances are actually NOT homogeneous

plot(gelpribiom_inf_aov, 2) # Normal distribution of the data NOT assumed
gelpribiom_inf_aov_residuals <- residuals(object = gelpribiom_inf_aov)
shapiro.test(x = gelpribiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
gelpribiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_gelpri2,
                             family = "Gamma")
summary(gelpribiom_inf_glm)
Anova(gelpribiom_inf_glm)
AIC(gelpribiom_inf_glm) # AIC(glm) = 443.637

# Check models
DHARMa::simulateResiduals(fittedModel = gelpribiom_inf_glm, plot = T) # Not such as bad fit.
gratia::appraise(gelpribiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) w/ some interactions ####
gelpribiom_inf_glm1 <- glm(biomass ~ inf.lvl * inf.percent + position, data = biom_gelpri2,
                          family = "Gamma")
summary(gelpribiom_inf_glm1)
Anova(gelpribiom_inf_glm1)
AIC(gelpribiom_inf_glm1) # AIC(glm) = 450.7665

# Check models
DHARMa::simulateResiduals(fittedModel = gelpribiom_inf_glm1, plot = T) # Better fit !
gratia::appraise(gelpribiom_inf_glm1, method = "simulate")

# Model #4: Additive GLM (Gamma) w/ some interactions ####
gelpribiom_inf_glm2 <- glm(biomass ~ inf.lvl * inf.percent + position + inf.lvl:position, data = biom_gelpri2,
                           family = "Gamma")
summary(gelpribiom_inf_glm2)
Anova(gelpribiom_inf_glm2)
AIC(gelpribiom_inf_glm2) # AIC(glm) = 446.0441

# Check models
DHARMa::simulateResiduals(fittedModel = gelpribiom_inf_glm2, plot = T) # Not a bad fit !
gratia::appraise(gelpribiom_inf_glm2, method = "simulate")

# Model #5: Additive GLM (Gamma) w/ some interactions ####
gelpribiom_inf_glm3 <- glm(biomass ~ inf.lvl * inf.percent + position + inf.lvl:position + inf.percent:position, 
                           data = biom_gelpri2,
                           family = "Gamma")
summary(gelpribiom_inf_glm3)
Anova(gelpribiom_inf_glm3)
AIC(gelpribiom_inf_glm3) # AIC(glm) = 443.637

# Check models
DHARMa::simulateResiduals(fittedModel = gelpribiom_inf_glm3, plot = T) # Not a good fit !
gratia::appraise(gelpribiom_inf_glm2, method = "simulate")

# Pairwise comparisons on biomass in HILLEC ####
emmeans(gelpribiom_inf_glm2, pairwise ~ inf.lvl | inf.percent)
emmeans(gelpribiom_inf_glm2, pairwise ~ position | inf.lvl)



# 7.6. Electra verticillata ====

# 7.6.1. ELEVER - Biomass vs Shell length ====
biom_elever <- biom_long_sp[biom_long_sp$epibiont %in% "Electra verticillata",]

biom_elever1 <- biom_elever %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_elever1)

hist(biom_elever1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
eleverbiom_length_aov <- aov(biomass ~ length_range, data = biom_elever1)
summary(eleverbiom_length_aov)

# Check ANOVA assumptions
plot(eleverbiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_elever1) # But, variances are actually homogeneous

plot(eleverbiom_length_aov, 2) # Normal distribution of the data NOT assumed
eleverbiom_length_aov_residuals <- residuals(object = eleverbiom_length_aov)
shapiro.test(x = eleverbiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
eleverbiom_length_glm <- glm(biomass ~ length_range, data = biom_elever1,
                             family = "Gamma")
summary(eleverbiom_length_glm)
Anova(eleverbiom_length_glm)
AIC(eleverbiom_length_glm) # AIC(glm) = 351.3309

# Check models
DHARMa::simulateResiduals(fittedModel = eleverbiom_length_glm, plot = T) # A very good fit !
gratia::appraise(eleverbiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in HILLEC ####
summary(glht(eleverbiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 7.6.2. ELEVER - Biomass vs (Infestation x Degradation x Position) ====
biom_elever <- biom_long_sp[biom_long_sp$epibiont %in% "Electra verticillata",]

biom_elever2 <- biom_elever %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_elever2)

hist(biom_elever2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
eleverbiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_elever2)
summary(eleverbiom_inf_aov)

# Check ANOVA assumptions
plot(eleverbiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_elever2) # But, variances are actually homogeneous

plot(eleverbiom_inf_aov, 2) # Normal distribution of the data NOT assumed
eleverbiom_inf_aov_residuals <- residuals(object = eleverbiom_inf_aov)
shapiro.test(x = eleverbiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
eleverbiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_elever2,
                             family = "Gamma")
summary(eleverbiom_inf_glm)
Anova(eleverbiom_inf_glm)
AIC(eleverbiom_inf_glm) # AIC(glm) = 411.8425

# Check models
DHARMa::simulateResiduals(fittedModel = eleverbiom_inf_glm, plot = T) # Not such a good fit.
gratia::appraise(eleverbiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
eleverbiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_elever2,
                          family = "Gamma")
summary(eleverbiom_inf_glm1)
Anova(eleverbiom_inf_glm1)
AIC(eleverbiom_inf_glm1) # AIC(glm) = 410.3575

# Check models
DHARMa::simulateResiduals(fittedModel = eleverbiom_inf_glm1, plot = T) # A very good fit !
gratia::appraise(eleverbiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in ELEVER ####
emmeans(eleverbiom_inf_glm1, pairwise ~ inf.lvl)
emmeans(eleverbiom_inf_glm1, pairwise ~ inf.percent)





# 7.7. Bryozoa ====

# 7.7.1. Bryozoa - Biomass vs Shell length ====
biom_bryozoa <- biom_long_sp[biom_long_sp$epibiont %in% "Bryozoa",]

biom_bryozoa1 <- biom_bryozoa %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bryozoa1)

hist(biom_bryozoa1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
bryozoabiom_length_aov <- aov(biomass ~ length_range, data = biom_bryozoa1)
summary(bryozoabiom_length_aov)

# Check ANOVA assumptions
plot(bryozoabiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_bryozoa1) # But, variances are actually homogeneous

plot(bryozoabiom_length_aov, 2) # Normal distribution of the data NOT assumed
bryozoabiom_length_aov_residuals <- residuals(object = bryozoabiom_length_aov)
shapiro.test(x = bryozoabiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
bryozoabiom_length_glm <- glm(biomass ~ length_range, data = biom_bryozoa1,
                             family = "Gamma")
summary(bryozoabiom_length_glm)
Anova(bryozoabiom_length_glm)
AIC(bryozoabiom_length_glm) # AIC(glm) = 401.1831

# Check models
DHARMa::simulateResiduals(fittedModel = bryozoabiom_length_glm, plot = T) # A good fit !
gratia::appraise(bryozoabiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Bryozoa ####
summary(glht(bryozoabiom_length_glm, linfct = mcp(length_range = "Tukey")))






# 7.7.2. Bryozoa - Biomass vs (Infestation x Degradation x Position) ====
biom_bryozoa <- biom_long_sp[biom_long_sp$epibiont %in% "Bryozoa",]

biom_bryozoa2 <- biom_bryozoa %>%
  dplyr::group_by(quadrat, specimen, valve, epibiont, inf.lvl, inf.percent, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_bryozoa2)

hist(biom_bryozoa2$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
bryozoabiom_inf_aov <- aov(biomass ~ inf.lvl * inf.percent * position, data = biom_bryozoa2)
summary(bryozoabiom_inf_aov)

# Check ANOVA assumptions
plot(bryozoabiom_inf_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ inf.lvl * inf.percent * position, data = biom_bryozoa2) # But, variances are actually homogeneous

plot(bryozoabiom_inf_aov, 2) # Normal distribution of the data NOT assumed
bryozoabiom_inf_aov_residuals <- residuals(object = bryozoabiom_inf_aov)
shapiro.test(x = bryozoabiom_inf_aov_residuals)  # Residuals NOT normally distributed

# Model #2: Saturated GLM (Gamma) ####
bryozoabiom_inf_glm <- glm(biomass ~ inf.lvl * inf.percent * position, data = biom_bryozoa2,
                              family = "Gamma")
summary(bryozoabiom_inf_glm)
Anova(bryozoabiom_inf_glm)
AIC(bryozoabiom_inf_glm) # AIC(glm) = 448.4989

# Check models
DHARMa::simulateResiduals(fittedModel = bryozoabiom_inf_glm, plot = T) # Not a good fit.
gratia::appraise(bryozoabiom_inf_glm, method = "simulate")

# Model #3: Additive GLM (Gamma) ####
bryozoabiom_inf_glm1 <- glm(biomass ~ inf.lvl + inf.percent + position, data = biom_bryozoa2,
                           family = "Gamma")
summary(bryozoabiom_inf_glm1)
Anova(bryozoabiom_inf_glm1)
AIC(bryozoabiom_inf_glm1) # AIC(glm) = 438.7751

# Check models
DHARMa::simulateResiduals(fittedModel = bryozoabiom_inf_glm1, plot = T) # Not a good fit.
gratia::appraise(bryozoabiom_inf_glm1, method = "simulate")

# Pairwise comparisons on biomass in Bryozoa ####
summary(glht(bryozoabiom_inf_glm1, linfct = mcp(inf.lvl = "Tukey")))




# 7.8. Perna perna ====

# 7.8.1. PERPER - Biomass vs Shell length ====
biom_perna <- biom_long_sp[biom_long_sp$epibiont %in% "Perna perna",]

biom_perna1 <- biom_perna %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_perna1)

hist(biom_perna1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
pernabiom_length_aov <- aov(biomass ~ length_range, data = biom_perna1)
summary(pernabiom_length_aov)

# Check ANOVA assumptions
plot(pernabiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_perna1) # But, variances are actually homogeneous

plot(pernabiom_length_aov, 2) # Normal distribution of the data NOT assumed
pernabiom_length_aov_residuals <- residuals(object = pernabiom_length_aov)
shapiro.test(x = pernabiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
pernabiom_length_glm <- glm(biomass ~ length_range, data = biom_perna1,
                              family = "Gamma")
summary(pernabiom_length_glm)
Anova(pernabiom_length_glm)
AIC(pernabiom_length_glm) # AIC(glm) = 406.9936

# Check models
DHARMa::simulateResiduals(fittedModel = pernabiom_length_glm, plot = T) # A good fit !
gratia::appraise(pernabiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Perna perna ####
summary(glht(pernabiom_length_glm, linfct = mcp(length_range = "Tukey")))





# 7.9. Spirobranchus kraussii ====

# 7.9.1. SPIKRA - Biomass vs Shell length ====
biom_spikra <- biom_long_sp[biom_long_sp$epibiont %in% "Spirobranchus kraussii",]

biom_spikra1 <- biom_spikra %>%
  dplyr::group_by(quadrat, specimen, epibiont, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_spikra1)

hist(biom_spikra1$biomass) # Gamma distribution

# Model #1: Regular ANOVA ####
spikrabiom_length_aov <- aov(biomass ~ length_range, data = biom_spikra1)
summary(spikrabiom_length_aov)

# Check ANOVA assumptions
plot(spikrabiom_length_aov, 1) # Homogeneity of variances NOT assumed
leveneTest(biomass ~ length_range, data = biom_spikra1) # But, variances are actually homogeneous

plot(spikrabiom_length_aov, 2) # Normal distribution of the data NOT assumed
spikrabiom_length_aov_residuals <- residuals(object = spikrabiom_length_aov)
shapiro.test(x = spikrabiom_length_aov_residuals)  # Residuals NOT normally distributed

# Model #2: GLM (Gamma) w/ shell length ####
spikrabiom_length_glm <- glm(biomass ~ length_range, data = biom_spikra1,
                            family = "Gamma")
summary(spikrabiom_length_glm)
Anova(spikrabiom_length_glm)
AIC(spikrabiom_length_glm) # AIC(glm) = 287.3831

# Check models
DHARMa::simulateResiduals(fittedModel = spikrabiom_length_glm, plot = T) # A good fit !
gratia::appraise(spikrabiom_length_glm, method = "simulate")

# Pairwise comparisons on biomass in Spirobranchus kraussii ####
summary(glht(spikrabiom_length_glm, linfct = mcp(length_range = "Tukey")))


