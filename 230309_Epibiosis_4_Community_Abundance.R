## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##              Script 4 - Community analysis on abundance
##
##
## Purpose of script: 
##
##
## Author: Alexia Dievart
##
## Date Created: 2023-01-18
## Dates Updated:
##
## Copyright (c) Alexia DIEVART 2023
## Email: alexia.dievart@hotmail.fr
##
## ---------------------------
##
## Notes: 
##   - 1 mussel = 2 valves
##   - Using Generalised Linear Latent Variable models (GLLVM) on abundance
##   - Using classic nMDS on biomass and percent cover
## ---------------------------





###############################################################################################
# Section: Session setup ----------------------------------------------------------------------
###############################################################################################

# Set working directory
setwd("D:/Little Bull/Etudes/phD - Rhodes/1.3_EPIBIOSIS/STATS")

# Load required packages
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(
  tidyverse,
  mvabund,
  MuMIn,
  gllvm
)





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



# Abundance ==============================================================================

# Select only environmental variables and abundance
abund_long <- epibio[,c(1:16,22)]
head(abund_long)
View(abund_long)

# Focus on higher groups ################################################################
abund_long_group <- abund_long %>%
  dplyr::group_by(quadrat, specimen, length_range, valve, inf.lvl, inf.percent,
                  higher.group, position, basibiont) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_long_group)

# Transform from long to wide format
abund_wide_group <- abund_long_group %>%
  pivot_wider(names_from = "higher.group", values_from = "abundance", values_fill = 0)
View(abund_wide_group)

# Store the columns containing group abundances into a new dataframe 
sp_abun_group <- abund_wide_group[, -c(1:8)]
head(sp_abun_group)

# Convert the group abundance data into an 'mvabund' object
sp_mva_group <- mvabund::mvabund(sp_abun_group)
head(sp_mva_group)

# Check that the object == 'mvabund'
mvabund::is.mvabund(sp_mva_group)





# Focus on species #################################################################
abund_long_sp <- abund_long %>%
  dplyr::group_by(quadrat, specimen, length_range, valve, inf.lvl, inf.percent,
                  epibiont, position, basibiont) %>%
  dplyr::summarise(
    abundance = sum(abundance)
  )
View(abund_long_sp)

# Transform from long to wide format
abund_wide_sp <- abund_long_sp %>%
  pivot_wider(names_from = "epibiont", values_from = "abundance", values_fill = 0)
View(abund_wide_sp)

# Store the columns containing species abundances into a new dataframe 
sp_abun_sp <- abund_wide_sp[, -c(1:8)]
head(sp_abun_sp)

# Convert the species abundance data into an 'mvabund' object
sp_mva_sp <- mvabund::mvabund(sp_abun_sp)
head(sp_mva_sp)

# Check that the object == 'mvabund'
mvabund::is.mvabund(sp_mva_sp)




###############################################################################################
# Section 2: Exploratory analyses -------------------------------------------------------------
###############################################################################################

# Abundance ==============================================================================

# Focus on higher group ################################################################

# Check the mean-variance relationship for the shell length range
meanvar.plot(sp_mva_group,
             col = as.factor(abund_wide_group$length_range),
             xlab = "Mean group abundance",
             ylab = "Variance")

# Add linear increase line
abline(a=0, 
       b=1, 
       untf = TRUE, 
       lty=3)
# Poisson model data should fall on this line (i.e., linear increase in variance with mean)

# Check the mean-variance relationship for the euendolithic infestation levels 
meanvar.plot(sp_mva_group,
             col = as.factor(abund_wide_group$inf.lvl),
             xlab = "Mean group abundance",
             ylab = "Variance")

# Check the mean-variance relationship for the shell degradation indexes 
meanvar.plot(sp_mva_group,
             col = as.factor(abund_wide_group$inf.percent),
             xlab = "Mean group abundance",
             ylab = "Variance")

# Check the mean-variance relationship for the position on the shell
meanvar.plot(sp_mva_pos_group,
             col = as.factor(abund_wide_pos_group$position),
             xlab = "Mean group abundance",
             ylab = "Variance")





#Focus on species ###################################################################

# Check the mean-variance relationship for the shell length range
meanvar.plot(sp_mva_sp,
             col = as.factor(abund_wide_sp$length_range),
             xlab = "Mean species abundance",
             ylab = "Variance")

# Add linear increase line
abline(a=0, 
       b=1, 
       untf = TRUE, 
       lty=3)
# Poisson model data should fall on this line (i.e., linear increase in variance with mean)

# Check the mean-variance relationship for the euendolithic infestation levels 
meanvar.plot(sp_mva_sp,
             col = as.factor(abund_wide_sp$inf.lvl),
             xlab = "Mean species abundance",
             ylab = "Variance")

# Check the mean-variance relationship for the shell degradation indexes 
meanvar.plot(sp_mva_sp,
             col = as.factor(abund_wide_sp$inf.percent),
             xlab = "Mean species abundance",
             ylab = "Variance")

# Check the mean-variance relationship for the position on the shell
meanvar.plot(sp_mva_sp,
             col = as.factor(abund_wide_sp$position),
             xlab = "Mean species abundance",
             ylab = "Variance")





###############################################################################################
# Section 3: Fit models  ----------------------------------------------------------------------
###############################################################################################

# Abundance  =============================================================================

# 3.1. Higher group - Abundance vs shell length ####
View(abund_wide_group)

# Model #1: log-link Poisson model on shell length ####
mod1_length <- mvabund::manyglm(
  sp_mva_group ~ abund_wide_group$length_range,
  family = "poisson"
)
mean(AIC(mod1_length)) # AIC(Poisson) = 1336.611

# Check model fit 
mvabund::plot.manyglm(mod1_length) # Not too much overdispersion, so not a bad fit.

# Model #2: simple negative binomial GLM on shell length ####
mod2_length <- manyglm(
  sp_mva_group ~ abund_wide_group$length_range,
  family = "negative.binomial"
)
mean(AIC(mod2_length)) # AIC(neg.bin) = 827.3046

# Check model fit 
mvabund::plot.manyglm(mod2_length) # Dispersion is much more constant, much better fit. 

# Multivariate hypothesis tests w/ negative binomial model ####
# - PIT-trap bootstrapping (Warton et al., 2017) used to calculate p-values 
# - cor.type  is set to "shrink" to account for potential correlations between
#   epibiotic species in their responses to predictor variables
#   - This process applies ridge regularisation (Warton 2008) to shrink 
#     the sample correlation matrix 
# - Wald's test statistics are calculated using a generalised estimating 
#   equation approach (GEE) to perform hypothesis testing 
mvabund::anova.manyglm(
  object = mod2_length,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "none",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)

# The results indicate that epibiotic community structure is significantly affected
# by shell length size (p = 0.001).

# Univariate hypothesis tests on shell length ####
# - These tests allow us to determine which epibiotic groups are driving the 
#   differences in epibiotic community structure
mvabund::anova.manyglm(
  object = mod2_length,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "adjusted",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)





# 3.2. Species - Abundance vs shell length ####
View(abund_wide_sp)

# Model #1: log-link Poisson model on shell length ####
mod1sp_length <- mvabund::manyglm(
  sp_mva_sp ~ abund_wide_sp$length_range,
  family = "poisson"
)
mean(AIC(mod1sp_length)) # AIC(Poisson) = 596.7904

# Check model fit 
mvabund::plot.manyglm(mod1sp_length) # Not too much overdispersion, so not a bad fit.

# Model #2: simple negative binomial GLM on shell length ####
mod2sp_length <- manyglm(
  sp_mva_sp ~ abund_wide_sp$length_range,
  family = "negative.binomial"
)
mean(AIC(mod2sp_length)) # AIC(neg.bin) = 375.448

# Check model fit 
mvabund::plot.manyglm(mod2sp_length) # Dispersion is much more constant, much better fit. 

# Multivariate hypothesis tests w/ negative binomial model ####
mvabund::anova.manyglm(
  object = mod2sp_length,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "none",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)

# The results indicate that epibiotic community structure is significantly affected
# by shell length size (p = 0.001).

# Univariate hypothesis tests on shell length ####
# - These tests allow us to determine which epibiotic groups are driving the 
#   differences in epibiotic community structure
mvabund::anova.manyglm(
  object = mod2sp_length,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "adjusted",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)





# 3.3. Higher group - Abundance vs (Infestation x Degradation) ####

# Model #1: Saturated log-link Poisson model on (infestation x degradation) ####
mod1_inf <- mvabund::manyglm(
  sp_mva_group ~ abund_wide_group$inf.lvl * abund_wide_group$inf.percent,
  family = "poisson"
)
mean(AIC(mod1_inf)) # AIC(Poisson) = 1370.842

# Check model fit 
mvabund::plot.manyglm(mod1_inf) # Not a bad fit. 

# Model #2: Additive log-link Poisson model on (infestation x degradation) ####
mod2_inf <- mvabund::manyglm(
  sp_mva_group ~ abund_wide_group$inf.lvl + abund_wide_group$inf.percent,
  family = "poisson"
)
mean(AIC(mod2_inf)) # AIC(Poisson) = 1356.934

# Check model fit 
mvabund::plot.manyglm(mod2_inf) # Much better dispersion.

# Model #3: Saturated negative binomial GLM on (infestation x degradation) ####
mod3_inf <- manyglm(
  sp_mva_group ~ abund_wide_group$inf.lvl * abund_wide_group$inf.percent,
  family = "negative.binomial"
)
mean(AIC(mod3_inf)) # AIC(neg.bin) = 859.7925

# Check model fit 
mvabund::plot.manyglm(mod3_inf) # Good fit too. 

# Model #4: Additive negative binomial GLM on (infestation x degradation) ####
mod4_inf <- manyglm(
  sp_mva_group ~ abund_wide_group$inf.lvl + abund_wide_group$inf.percent,
  family = "negative.binomial"
)
mean(AIC(mod4_inf)) # AIC(neg.bin) = 831.6194

# Check model fit 
mvabund::plot.manyglm(mod4_inf) # Very good fit too. 

# Multivariate hypothesis tests w/ negative binomial model ####
# - PIT-trap bootstrapping (Warton et al., 2017) used to calculate p-values 
# - cor.type  is set to "shrink" to account for potential correlations between
#   epibiotic species in their responses to predictor variables
#   - This process applies ridge regularisation (Warton 2008) to shrink 
#     the sample correlation matrix 
# - Wald's test statistics are calculated using a generalised estimating 
#   equation approach (GEE) to perform hypothesis testing 
mvabund::anova.manyglm(
  object = mod4_inf,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "none",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)

# The results indicate that epibiotic community structure is significantly affected
# by euendolithic infestation (p = 0.001) and shell degradation (p = 0.001).

# Univariate hypothesis tests on infestation and degradation ####
# - These tests allow us to determine which epibiotic groups are driving the 
#   differences in epibiotic community structure
mvabund::anova.manyglm(
  object = mod4_inf,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "adjusted",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)





# 3.4. Species - Abundance vs (Infestation x Degradation) ####

# Model #1: Saturated log-link Poisson model on (infestation x degradation) ####
mod1sp_inf <- mvabund::manyglm(
  sp_mva_sp ~ abund_wide_sp$inf.lvl * abund_wide_sp$inf.percent,
  family = "poisson"
)
mean(AIC(mod1sp_inf)) # AIC(Poisson) = 627.8448

# Check model fit 
mvabund::plot.manyglm(mod1sp_inf) # Not a bad fit.

# Model #2: Additive log-link Poisson model on (infestation x degradation) ####
mod2sp_inf <- mvabund::manyglm(
  sp_mva_sp ~ abund_wide_sp$inf.lvl + abund_wide_sp$inf.percent,
  family = "poisson"
)
mean(AIC(mod2sp_inf)) # AIC(Poisson) = 603.9483

# Check model fit 
mvabund::plot.manyglm(mod2sp_inf) # Much better.

# Model #3: Saturated negative binomial GLM on (infestation x degradation) ####
mod3sp_inf <- manyglm(
  sp_mva_sp ~ abund_wide_sp$inf.lvl * abund_wide_sp$inf.percent,
  family = "negative.binomial"
)
mean(AIC(mod3sp_inf)) # AIC(neg.bin) = 409.3665

# Check model fit 
mvabund::plot.manyglm(mod3sp_inf) # Good fit too. 

# Model #4: Additive negative binomial GLM on (infestation x degradation) ####
mod4sp_inf <- manyglm(
  sp_mva_sp ~ abund_wide_sp$inf.lvl + abund_wide_sp$inf.percent,
  family = "negative.binomial"
)
mean(AIC(mod4sp_inf)) # AIC(neg.bin) = 377.9602

# Check model fit 
mvabund::plot.manyglm(mod4sp_inf) # Good fit too. 

# Multivariate hypothesis tests w/ negative binomial model ####
# - PIT-trap bootstrapping (Warton et al., 2017) used to calculate p-values 
# - cor.type  is set to "shrink" to account for potential correlations between
#   epibiotic species in their responses to predictor variables
#   - This process applies ridge regularisation (Warton 2008) to shrink 
#     the sample correlation matrix 
# - Wald's test statistics are calculated using a generalised estimating 
#   equation approach (GEE) to perform hypothesis testing 
mvabund::anova.manyglm(
  object = mod4sp_inf,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "none",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)

# The results indicate that epibiotic community structure is significantly affected
# by euendolithic infestation (p = 0.001) and shell degradation (p = 0.001).

# Univariate hypothesis tests on infestation and degradation ####
# - These tests allow us to determine which epibiotic groups are driving the 
#   differences in epibiotic community structure
mvabund::anova.manyglm(
  object = mod4sp_inf,
  resamp = "pit.trap", 
  test = "wald",            # "LR" can only be used when cor.type = "I"
  p.uni = "adjusted",
  cor.type = "shrink",
  rep.seed = 2012, 
  nBoot = 999
)