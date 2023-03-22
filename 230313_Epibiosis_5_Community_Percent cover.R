## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##              Script 5 - Community analysis on percent cover
##
##
## Purpose of script: 
##
##
## Author: Alexia Dievart
##
## Date Created: 2023-03-10
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

# Load packages
library(pacman)
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
               vegan)

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
colors_position1 <- c("darkgrey", "brown3")
colors_degrad <- c("chocolate4", "chocolate", "darkgoldenrod3", "darkgoldenrod1", "bisque2", "cornsilk2")





###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Load raw data ===============================================================================
epibio <- read.csv("epibiosis.csv", h=T, dec=",", sep=";")
head(epibio)

# Curate raw data set =========================================================================

# Convert variables into factor
factor_epibio <- c("location", "quadrat", "specimen", "mussel.species", "sex", "valve", "inf.lvl", "inf.percent",
                   "higher.group", "epibiont", "position", "basibiont")
epibio[,factor_epibio] <- lapply(epibio[,factor_epibio], factor)

# Recode inf.percent into shell degradation index
epibio$inf.percent <- recode_factor(epibio$inf.percent, "0" = "1", "0-25" = "2", "25-50" = "3", "50-75" = "4",
                                    "75-100" = "5", "100" = "6")

# Create a new "shell length range" variable in the raw data set
length_range_bins <- c(30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
length_range_labels <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)
epibio$length_range <- cut(epibio$shell.length, breaks = length_range_bins,
                           labels = length_range_labels)

# Exclude epibionts recorded on Mytilus galloprovincialis
epibio <- epibio[! epibio$mussel.species == "Mytilus galloprovincialis",]

# Exclude epibionts w/ cover.percent = 0 or is NA
epibio <- epibio[! epibio$cover.percent == 0,]
epibio <- epibio[! is.na(epibio$cover.percent),]
View(epibio)

# Exclude mobile epibionts 
epibio <- epibio[! epibio$mobility == "mobile",]

# Exclude IIary epibionts and epibionts in byssus
epibio <- epibio[! epibio$position %in% c("II epibiosis", "byssus"),]
View(epibio)

# Isolate percent cover =======================================================================
cover_long <- epibio[,c(1:15, 18, 22)]
View(cover_long)

# Focus on higher groups ####
cover_long_group <- cover_long %>%
  dplyr::group_by(quadrat, specimen, valve, length_range, inf.lvl, inf.percent,
                  higher.group, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_long_group)

# Transform from long to wide format
cover_wide_group <- cover_long_group %>%
  pivot_wider(names_from = "higher.group", values_from = "cover", values_fill = 0)
View(cover_wide_group)

# Focus on epibiotic species ####
cover_long_sp <- cover_long %>%
  dplyr::group_by(quadrat, specimen, valve, length_range, inf.lvl, inf.percent,
                  epibiont, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_long_sp)

# Transform from long to wide format
cover_wide_sp <- cover_long_sp %>%
  pivot_wider(names_from = "epibiont", values_from = "cover", values_fill = 0)
View(cover_wide_sp)

# Store percent cover and environmental variables into separate data frames ===================
env_group <- cover_wide_group[, c(1:7)]
cover_group <- cover_wide_group[, c(8:18)]
View(cover_group)

env_sp <- cover_wide_sp[, c(1:7)]
cover_sp <- cover_wide_sp[, c(8:43)]
View(cover_sp)






###############################################################################################
# Section 1: Community analysis on higher groups ------------------------------------------------
###############################################################################################

# Calculate the distance matrix ===============================================================
# Compute distance matrix using Bray-Curtis on transformed abundances
range(cover_group)
range(cover_group^0.5) # Range < 10
range(cover_group^0.25) 
# If working on community data, down-weigth the abundant species to a range between 0 and 10
cover_dist_group <- vegdist(cover_group^0.5, method ='bray')
cover_dist_group

# Calculate nMDS ==============================================================================

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.05 --> Good fit

cover_bray_group <- metaMDS(cover_group^0.5, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to not output all of these random starts
cover_bray_group # Stress = 0.0735 shows this nMDS is fair.

# Exploring the results of the nMDS ===========================================================
plot(cover_bray_group, type = "t", xlim = c(-1.0, 5.0), ylim = c(-1.0, 1.0))
stressplot(cover_bray_group)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(cover_bray_group)
plot(cover_bray_group, type = "t", main = "goodness of fit")
points(cover_bray_group, display = "sites", cex = gof*100) # sites refers to lines

# What parameters are significant ?
cover_fit_group <- envfit(cover_bray_group, env_group, permu = 999) # Caution: envfit does not allow missing values
cover_fit_group
# Quadrat, infestation levels and position have a significant effect on the goodness of fit of the nMDS model.

library(RVAideMemoire)
cover_fit_group_adj <- cover_fit_group
pvals_group_adj <- p.adjust(cover_fit_group$factors$pvals, method = 'bonferroni')
cover_fit_group_adj$factors$pvals <- pvals_group_adj
cover_fit_group_adj
# W/ a Bonferroni correction, only position had still a significant effect on the goodness of fit of the nMDS model.

# Pairwise comparisons ==========================================================================

pairwise.factorfit(cover_bray_group, env_group$length_range, p.method = "bonferroni")
# - 50 vs 70 is nearly significant.
# - 70 vs 80, and 70 vs 100 is significant.

pairwise.factorfit(cover_bray_group, env_group$inf.lvl, p.method = "bonferroni")
# No significant difference between infestation levels

pairwise.factorfit(cover_bray_group, env_group$inf.percent, p.method = "bonferroni")
# - 1 vs 4 is significant.

pairwise.factorfit(cover_bray_group, env_group$position, p.method = "bonferroni")
# - infested vs non-infested is significant.

plot(cover_bray_group, display = "sites", xlim = c(-1.0, 5.0), ylim = c(-1.0, 1.0))
plot(cover_fit_group_adj, p.max = 0.05) 

# Plot nMDS w/ infestation levels and degradation =================================================

ordiplot(cover_bray_group, type = "none", xlim = c(-0, 5), ylim = c(-0.5, 0.5))
points(cover_bray_group, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
orditorp(cover_bray_group, display = "species", col = "black", air = 0.4, cex = 0.8)
legend(0, - 0.2, title = "SD index",
       legend = levels(env_group$inf.percent),
       pch = c(16, 17, 18, 8, 5, 3), 
       col = c("black", 'black', "black", "black", "black", "black"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(1, - 0.2, title = "Euendolithic infestation",
       legend = levels(env_group$inf.lvl),
       pch = 16,
       col = c("coral4", "brown3", "indianred1", "pink3", "darkgrey"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-0.3, 1.5, "Stress = 0.0735", bty = "n", cex = 1)

# Run PERMANOVA ==================================================================================
cover_group_pmv <- adonis2(cover_group^0.5 ~ length_range + inf.lvl + inf.percent + position + inf.lvl:position +
                             inf.percent:position, data = env_group,
                           permutations = 999, method = 'bray')
cover_group_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities differ significantly between :
#   - Shell length sizes
#   - Infestation levels
#   - Shell degradation index
#   - Position
#   - Position x Inf.percent

# Plot permuted F-values
densityplot(permustats(cover_group_pmv))

# Average distance to median and significance ==================================================

# For shell length sizes ####
cover_disp_group_length <- betadisper(cover_dist_group, group = env_group$length_range)
boxplot(cover_disp_group_length)
# - Less dispersion for shell length ranges 60, 70, 80, 90, 100 compared to the others
cover_disp_group_length

# Is there a significant difference in dispersion between shell length sizes ?
anova(cover_disp_group_length)
permutest(cover_disp_group_length)
# Significant difference in dispersion between shell length sizes
col_length <- hcl(h = seq(180, 270, length = 9), c = 70, l = 65)

plot(cover_disp_group_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = T, cex = 0.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
round(100 * cover_disp_group_length$eig / sum(cover_disp_group_length$eig), 2)

# For infestation levels ####
cover_disp_group_inf <- betadisper(cover_dist_group, group = env_group$inf.lvl)
boxplot(cover_disp_group_inf)
# - Less dispersion for infestation levels B and C
cover_disp_group_inf

plot(cover_disp_group_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
round(100 * cover_disp_group_inf$eig / sum(cover_disp_group_inf$eig), 2)

# For degradation index ####
cover_disp_group_deg <- betadisper(cover_dist_group, group = env_group$inf.percent)
boxplot(cover_disp_group_deg)
# - Less dispersion for degradation index 3 and 4
cover_disp_group_deg

plot(cover_disp_group_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = T, cex = 1.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
round(100 * cover_disp_group_deg$eig / sum(cover_disp_group_deg$eig), 2)

# For position ####
cover_disp_group_pos <- betadisper(cover_dist_group, group = env_group$position)
boxplot(cover_disp_group_pos)
# - Less dispersion for infested portion of the shell.
cover_disp_group_pos

plot(cover_disp_group_pos, hull = F, ellipse = T, segments = T, seg.col = colors_position1, seg.lwd = 1,
     pch = c(16, 16),
     col = colors_position1, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
round(100 * cover_disp_group_pos$eig / sum(cover_disp_group_pos$eig), 2)

## SIMPER analyses ===================================================================================
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)

# SIMPER on shell length ####
cover_simp_length <- simper(cover_group^0.5, group = env_group$length_range)
cover_simp_length
summary(cover_simp_length)

# SIMPER on infestation
cover_simp_inf <- simper(cover_group^0.5, group = env_group$inf.lvl)
cover_simp_inf
summary(cover_simp_inf)

# SIMPER on degradation
cover_simp_deg <- simper(cover_group^0.5, group = env_group$inf.percent)
cover_simp_deg
summary(cover_simp_deg)

# SIMPER on position
cover_simp_pos <- simper(cover_group^0.5, group = env_group$position)
cover_simp_pos
summary(cover_simp_pos)





###############################################################################################
# Section 2: Community analysis on epibiotic species ------------------------------------------
###############################################################################################

# Calculate the distance matrix ===============================================================
# Compute distance matrix using Bray-Curtis on transformed abundances
range(cover_sp)
range(cover_sp^0.5) # Range < 10
range(cover_sp^0.25) 
# If working on community data, down-weigth the abundant species to a range between 0 and 10
cover_dist_sp <- vegdist(cover_sp^0.5, method ='bray')
cover_dist_group

# Calculate nMDS ==============================================================================

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.05 --> Good fit

cover_bray_sp <- metaMDS(cover_sp^0.5, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to not output all of these random starts
cover_bray_sp # Stress = 0.003 shows this nMDS is fair.

# Exploring the results of the nMDS ===========================================================
plot(cover_bray_sp, type = "t")
stressplot(cover_bray_sp)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(cover_bray_sp)
plot(cover_bray_sp, type = "t", main = "goodness of fit")
points(cover_bray_sp, display = "sites", cex = gof*100)

cover_fit_sp <- envfit(cover_bray_sp, env_sp, permu = 999) # Caution: envfit does not allow missing values
cover_fit_sp
# Degradation has a significant effect on the goodness of fit of the nMDS model.
# Quadrat and infestation levels have nearly a significant effect on the goodness of fit of the nMDS model.

library(RVAideMemoire)
cover_fit_sp_adj <- cover_fit_sp
pvals_sp_adj <- p.adjust(cover_fit_sp$factors$pvals, method = 'bonferroni')
cover_fit_sp_adj$factors$pvals <- pvals_sp_adj
cover_fit_sp_adj
# W/ a Bonferroni correction, it is not significant anymore.

# Pairwise comparisons

pairwise.factorfit(cover_bray_sp, env_sp$length_range, p.method = "bonferroni")
# - 60 vs 110 are significantly different.

pairwise.factorfit(cover_bray_sp, env_sp$inf.lvl, p.method = "bonferroni")
# - Levels B and C are significantly different from level E

pairwise.factorfit(cover_bray_sp, env_sp$inf.percent, p.method = "bonferroni")
# - No significant differences

pairwise.factorfit(cover_bray_sp, env_sp$position, p.method = "bonferroni")
# - No significant differences

plot(cover_bray_sp, display = "sites")
plot(cover_fit_sp_adj, p.max = 0.05) 

# Plot nMDS w/ infestation levels and sites =============================================
ordiplot(cover_bray_sp, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(cover_bray_sp, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
legend(0, - 0.2, title = "SD index",
       legend = levels(env_group$inf.percent),
       pch = c(16, 17, 18, 8, 5, 3), 
       col = c("black", 'black', "black", "black", "black", "black"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(1, - 0.2, title = "Euendolithic infestation",
       legend = levels(env_group$inf.lvl),
       pch = 16,
       col = c("coral4", "brown3", "indianred1", "pink3", "darkgrey"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-0.3, 1.5, "Stress = 0.003", bty = "n", cex = 1)

# Run PERMANOVA ============================================================================
cover_sp_pmv <- adonis2(cover_sp^0.5 ~ length_range + inf.lvl + inf.percent + position + inf.lvl:position +
                             inf.percent:position, data = env_sp,
                           permutations = 999, method = 'bray')
cover_sp_pmv
# Communities differ between :
#   - Shell length sizes
#   - Infestation
#   - Degradation
#   - Position
#   - Degradation x Position (nearly significant)

# Plot permuted F-values
densityplot(permustats(cover_sp_pmv))

# Average distance to median and significance ==================================================

# For shell length sizes ####
cover_disp_sp_length <- betadisper(cover_dist_sp, group = env_sp$length_range)
boxplot(cover_disp_sp_length)
# - Less dispersion for shell length ranges 60, 70, 80, 90, 100 compared to the others
cover_disp_sp_length

# Is there a significant difference in dispersion between shell length sizes ?
anova(cover_disp_sp_length)
permutest(cover_disp_sp_length)
# Significant difference in dispersion between shell length sizes
col_length <- hcl(h = seq(180, 270, length = 9), c = 70, l = 65)

plot(cover_disp_sp_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = T, cex = 0.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
round(100 * cover_disp_sp_length$eig / sum(cover_disp_sp_length$eig), 2)

# For infestation levels ####
cover_disp_sp_inf <- betadisper(cover_dist_sp, group = env_sp$inf.lvl)
boxplot(cover_disp_sp_inf)
# - Less dispersion for infestation levels B and C
cover_disp_sp_inf

plot(cover_disp_sp_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
round(100 * cover_disp_sp_inf$eig / sum(cover_disp_sp_inf$eig), 2)

# For degradation index ####
cover_disp_sp_deg <- betadisper(cover_dist_sp, group = env_sp$inf.percent)
boxplot(cover_disp_sp_deg)
# - Less dispersion for degradation index 3 and 4
cover_disp_sp_deg

plot(cover_disp_sp_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = T, cex = 1.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
round(100 * cover_disp_sp_deg$eig / sum(cover_disp_sp_deg$eig), 2)

# For position ####
cover_disp_sp_pos <- betadisper(cover_dist_sp, group = env_sp$position)
boxplot(cover_disp_sp_pos)
# - Less dispersion for infested portion of the shell.
cover_disp_sp_pos

plot(cover_disp_sp_pos, hull = F, ellipse = T, segments = T, seg.col = colors_position1, seg.lwd = 1,
     pch = c(16, 16),
     col = colors_position1, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)

round(100 * cover_disp_sp_deg$eig / sum(cover_disp_sp_deg$eig), 2)

## SIMPER analyses ===================================================================================
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)

# SIMPER on shell length ####
cover_simp_length1 <- simper(cover_sp^0.5, group = env_sp$length_range)
cover_simp_length1
summary(cover_simp_length1)

summary_cover_simp_length1 <- summary(cover_simp_length1)
View(summary_cover_simp_length1$'90_80')
View(summary_cover_simp_length1$'90_70')
View(summary_cover_simp_length1$'90_40')
View(summary_cover_simp_length1$'90_60')
View(summary_cover_simp_length1$'90_50')
View(summary_cover_simp_length1$'90_110')
View(summary_cover_simp_length1$'90_100')
View(summary_cover_simp_length1$'90_120')
View(summary_cover_simp_length1$'80_70')
View(summary_cover_simp_length1$'80_40')
View(summary_cover_simp_length1$'80_60')
View(summary_cover_simp_length1$'80_50')
View(summary_cover_simp_length1$'80_110')
View(summary_cover_simp_length1$'80_100')
View(summary_cover_simp_length1$'80_120')
View(summary_cover_simp_length1$'70_40')
View(summary_cover_simp_length1$'70_60')
View(summary_cover_simp_length1$'70_50')
View(summary_cover_simp_length1$'70_110')
View(summary_cover_simp_length1$'70_100')
View(summary_cover_simp_length1$'70_120')
View(summary_cover_simp_length1$'40_60')
View(summary_cover_simp_length1$'40_50')
View(summary_cover_simp_length1$'40_110')
View(summary_cover_simp_length1$'40_100')
View(summary_cover_simp_length1$'40_120')
View(summary_cover_simp_length1$'60_50')
View(summary_cover_simp_length1$'60_110')
View(summary_cover_simp_length1$'60_100')
View(summary_cover_simp_length1$'60_120')
View(summary_cover_simp_length1$'50_110')
View(summary_cover_simp_length1$'50_100')
View(summary_cover_simp_length1$'50_120')
View(summary_cover_simp_length1$'110_100')
View(summary_cover_simp_length1$'110_120')
View(summary_cover_simp_length1$'100_120')


# SIMPER on infestation
cover_simp_inf1 <- simper(cover_sp^0.5, group = env_sp$inf.lvl)
cover_simp_inf1
summary(cover_simp_inf1)

# SIMPER on degradation
cover_simp_deg1 <- simper(cover_sp^0.5, group = env_sp$inf.percent)
cover_simp_deg1
summary(cover_simp_deg1)

# SIMPER on degradation
cover_simp_pos1 <- simper(cover_sp^0.5, group = env_sp$position)
cover_simp_pos1
summary(cover_simp_pos1)
