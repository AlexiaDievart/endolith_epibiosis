## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##              Script 6 - Community analysis on biomass
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
colors_position2 <- c("cadetblue", "darkgrey", "brown3")
colors_degrad <- c("chocolate4", "chocolate", "darkgoldenrod3", "darkgoldenrod1", "bisque2", "cornsilk2")






###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Load raw data ===============================================================================
epibio <- read.csv("epibiosis.csv", h = T, dec = ",", sep = ";")
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
View(epibio)

# Exclude epibionts w/ biomass = 0 or is NA
epibio <- epibio[! epibio$biomass == 0,]
epibio <- epibio[! is.na(epibio$biomass),]
View(epibio)

# Exclude epibionts in byssus
epibio <- epibio[! epibio$position %in% "byssus",]
View(epibio)

# Isolate biomass =======================================================================
biom_long <- epibio[,c(1:15, 19, 22)]
View(biom_long)

# Focus on higher groups ####
biom_long_group <- biom_long %>%
  dplyr::group_by(quadrat, specimen, valve, length_range, inf.lvl, inf.percent,
                  higher.group, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_long_group)

# Transform from long to wide format
biom_wide_group <- biom_long_group %>%
  pivot_wider(names_from = "higher.group", values_from = "biomass", values_fill = 0)
View(biom_wide_group)

# Focus on epibiotic species ####
biom_long_sp <- biom_long %>%
  dplyr::group_by(quadrat, specimen, valve, length_range, inf.lvl, inf.percent,
                  epibiont, position) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_long_sp)

# Transform from long to wide format
biom_wide_sp <- biom_long_sp %>%
  pivot_wider(names_from = "epibiont", values_from = "biomass", values_fill = 0)
View(biom_wide_sp)

# Store percent cover and environmental variables into separate data frames ===================
env_group <- biom_wide_group[, c(1:7)]
View(env_group)
biom_group <- biom_wide_group[, c(8:24)]
View(biom_group)

env_sp <- biom_wide_sp[, c(1:7)]
biom_sp <- biom_wide_sp[, c(8:43)]
View(biom_sp)





###############################################################################################
# Section 1: Community analysis on higher groups ------------------------------------------------
###############################################################################################

# Calculate the distance matrix ===============================================================
# Compute distance matrix using Bray-Curtis on biomass
range(biom_group)
range(biom_group^0.5) 
range(biom_group^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
biom_dist_group <- vegdist(biom_group^0.25, method = 'bray')
biom_dist_group




# Calculate nMDS ==============================================================================

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.05 --> Good fit

biom_bray_group <- metaMDS(biom_group^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to not output all of these random starts
biom_bray_group # Stress = 0.004 is a very good fit.




# Exploring the results of the nMDS ===========================================================

plot(biom_bray_group, type = "t", xlim = c(-3, 2), ylim = c(-1.0, 1.0))
stressplot(biom_bray_group)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.

gof <- goodness(biom_bray_group)
plot(biom_bray_group, type = "t", main = "goodness of fit", xlim = c(-3, 2), ylim = c(-1.0, 1.0))
points(biom_bray_group, display = "sites", cex = gof*100) # sites refers to lines

# What parameters are significant ?
biom_fit_group <- envfit(biom_bray_group, env_group, permu = 999) # Caution: envfit does not allow missing values
biom_fit_group
# Valve have a nearly significant effect on the goodness of fit of the nMDS model.

library(RVAideMemoire)
biom_fit_group_adj <- biom_fit_group
pvals_group_adj <- p.adjust(biom_fit_group$factors$pvals, method = 'bonferroni')
biom_fit_group_adj$factors$pvals <- pvals_group_adj
biom_fit_group_adj
# W/ a Bonferroni correction, there is no significant effect on the goodness of fit of the nMDS model.





# Pairwise comparisons ==========================================================================

pairwise.factorfit(biom_bray_group, env_group$length_range, p.method = "bonferroni")
# - 60 vs 80/110 is significant
# - 70 vs 80 is nearly significant

pairwise.factorfit(biom_bray_group, env_group$inf.lvl, p.method = "bonferroni")
# - B vs C, C vs E are significant

pairwise.factorfit(biom_bray_group, env_group$inf.percent, p.method = "bonferroni")
# - 3 vs 6 is significant

pairwise.factorfit(biom_bray_group, env_group$position, p.method = "bonferroni")
# - II epibiosis vs non-infested is nearly significant

plot(biom_bray_group, display = "sites", xlim = c(-3, 2), ylim = c(-1.0, 1.0))
plot(biom_fit_group_adj, p.max = 0.05) 





# Plot nMDS w/ infestation levels and degradation =================================================

ordiplot(biom_bray_group, type = "none", xlim = c(-3, 2), ylim = c(-1.5, 1.5))
points(biom_bray_group, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
orditorp(biom_bray_group, display = "species", col = "black", air = 0.4, cex = 0.8)
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
legend(-0.3, 1.5, "Stress = 0.0044", bty = "n", cex = 1)





# Run PERMANOVA ==================================================================================
biom_group_pmv <- adonis2(biom_group^0.25 ~ length_range + inf.lvl + inf.percent + position + inf.lvl:position +
                             inf.percent:position, data = env_group,
                           permutations = 999, method = 'bray')
biom_group_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities differ significantly between :
#   - Shell length sizes
#   - Infestation levels
#   - Shell degradation index
#   - Position
#   - Infestation x Position

# Plot permuted F-values
densityplot(permustats(biom_group_pmv))





# Average distance to median and significance ==================================================

# For shell length sizes ####
biom_disp_group_length <- betadisper(biom_dist_group, group = env_group$length_range)
boxplot(biom_disp_group_length)
# - Less dispersion for shell length ranges 80 and 90 compared to the others
biom_disp_group_length

# Is there a significant difference in dispersion between shell length sizes ?
anova(biom_disp_group_length)
permutest(biom_disp_group_length)
# Significant difference in dispersion between shell length sizes

col_length <- hcl(h = seq(180, 270, length = 9), c = 70, l = 65)

round(100 * biom_disp_group_length$eig / sum(biom_disp_group_length$eig), 2)

plot(biom_disp_group_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = T, cex = 0.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.7, 0.7),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)





# For infestation levels ####
biom_disp_group_inf <- betadisper(biom_dist_group, group = env_group$inf.lvl)
boxplot(biom_disp_group_inf)
# - Less dispersion for infestation levels C and D compared to others
biom_disp_group_inf

round(100 * biom_disp_group_inf$eig / sum(biom_disp_group_inf$eig), 2)

plot(biom_disp_group_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)





# For degradation index ####
biom_disp_group_deg <- betadisper(biom_dist_group, group = env_group$inf.percent)
boxplot(biom_disp_group_deg)
# - More dispersion for shell degradation index 2
biom_disp_group_deg

round(100 * biom_disp_group_deg$eig / sum(biom_disp_group_deg$eig), 2)

plot(biom_disp_group_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = T, cex = 1.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)





# For position ####
biom_disp_group_pos <- betadisper(biom_dist_group, group = env_group$position)
boxplot(biom_disp_group_pos)
# - Less dispersion for II epibiosis and infested portion of the shell.
biom_disp_group_pos

round(100 * biom_disp_group_pos$eig / sum(biom_disp_group_pos$eig), 2)

plot(biom_disp_group_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position2, seg.lwd = 1,
     pch = c(16, 16, 16),
     col = colors_position2, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)





# SIMPER analyses ===================================================================================

# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)

# SIMPER on shell length ####
biom_simp_length <- simper(biom_group^0.25, group = env_group$length_range)
biom_simp_length
summary(biom_simp_length)




# SIMPER on infestation
biom_simp_inf <- simper(biom_group^0.25, group = env_group$inf.lvl)
biom_simp_inf
summary(biom_simp_inf)





# SIMPER on degradation
biom_simp_deg <- simper(biom_group^0.25, group = env_group$inf.percent)
biom_simp_deg
summary(biom_simp_deg)





# SIMPER on position
biom_simp_pos <- simper(biom_group^0.25, group = env_group$position)
biom_simp_pos
summary(biom_simp_pos)





###############################################################################################
# Section 2: Community analysis on epibiotic species ------------------------------------------
###############################################################################################

# Calculate the distance matrix ===============================================================

# Compute distance matrix using Bray-Curtis on transformed abundances
range(biom_sp)
range(biom_sp^0.5) 
range(biom_sp^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
biom_dist_sp <- vegdist(biom_sp^0.25, method ='bray')
biom_dist_sp





# Calculate nMDS ==============================================================================

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.05 --> Good fit

biom_bray_sp <- metaMDS(biom_sp^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to not output all of these random starts
biom_bray_sp # Stress = 0.007 shows this nMDS is a good fit.





# Exploring the results of the nMDS ===========================================================

plot(biom_bray_sp, type = "t", xlim = c(-150, 250), ylim = c(-25, 100))
stressplot(biom_bray_sp)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.

gof <- goodness(biom_bray_sp)
plot(biom_bray_sp, type = "t", main = "goodness of fit", xlim = c(-150, 250), ylim = c(-25, 100))
points(biom_bray_sp, display = "sites", cex = gof*100)

biom_fit_sp <- envfit(biom_bray_sp, env_sp, permu = 999) # Caution: envfit does not allow missing values
biom_fit_sp
# Infestation has a nearly significant effect on the goodness of fit of the nMDS model.

library(RVAideMemoire)
biom_fit_sp_adj <- biom_fit_sp
pvals_sp_adj <- p.adjust(biom_fit_sp$factors$pvals, method = 'bonferroni')
biom_fit_sp_adj$factors$pvals <- pvals_sp_adj
biom_fit_sp_adj
# W/ a Bonferroni correction, it is not significant anymore.





# Pairwise comparisons ####

pairwise.factorfit(biom_bray_sp, env_sp$length_range, p.method = "bonferroni")
# - Nothing is significant.

pairwise.factorfit(biom_bray_sp, env_sp$inf.lvl, p.method = "bonferroni")
# - Nothing is significant.

pairwise.factorfit(biom_bray_sp, env_sp$inf.percent, p.method = "bonferroni")
# - Nothing is significant.

pairwise.factorfit(biom_bray_sp, env_sp$position, p.method = "bonferroni")
# - Nothing is significant.

plot(biom_bray_sp, display = "sites")
plot(biom_fit_sp_adj, p.max = 0.05) 





# Plot nMDS w/ infestation levels and sites =============================================
ordiplot(biom_bray_sp, type = "none", xlim = c(-50, 50), ylim = c(-35, 30))
points(biom_bray_sp, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_sp$inf.percent],
       col = colors_inflvl[env_sp$inf.lvl])
legend(0, - 0.2, title = "SD index",
       legend = levels(env_sp$inf.percent),
       pch = c(16, 17, 18, 8, 5, 3), 
       col = c("black", 'black', "black", "black", "black", "black"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(1, - 0.2, title = "Euendolithic infestation",
       legend = levels(env_sp$inf.lvl),
       pch = 16,
       col = c("coral4", "brown3", "indianred1", "pink3", "darkgrey"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-0.3, 1.5, "Stress = 0.007", bty = "n", cex = 1)





# Run PERMANOVA ============================================================================
biom_sp_pmv <- adonis2(biom_sp^0.25 ~ length_range + inf.lvl + inf.percent + position + inf.lvl:position +
                          inf.percent:position, data = env_sp,
                        permutations = 999, method = 'bray')
biom_sp_pmv
# Communities differ between :
#   - Shell length sizes
#   - Infestation
#   - Degradation
#   - Position
#   - Infestation x Position

# Plot permuted F-values
densityplot(permustats(biom_sp_pmv))






# Average distance to median and significance ==================================================

# For shell length sizes ####
biom_disp_sp_length <- betadisper(biom_dist_sp, group = env_sp$length_range)
boxplot(biom_disp_sp_length)
# - Less dispersion for shell length ranges 70, 80, 90, 100 compared to the others
biom_disp_sp_length

# Is there a significant difference in dispersion between shell length sizes ?
anova(biom_disp_sp_length)
permutest(biom_disp_sp_length)
# Significant difference in dispersion between shell length sizes
col_length <- hcl(h = seq(180, 270, length = 9), c = 70, l = 65)

round(100 * biom_disp_sp_length$eig / sum(biom_disp_sp_length$eig), 2)

plot(biom_disp_sp_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = T, cex = 0.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)





# For infestation levels ####
biom_disp_sp_inf <- betadisper(biom_dist_sp, group = env_sp$inf.lvl)
boxplot(biom_disp_sp_inf)
# - Less dispersion for infestation levels C and D
biom_disp_sp_inf

round(100 * biom_disp_sp_inf$eig / sum(biom_disp_sp_inf$eig), 2)

plot(biom_disp_sp_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)




# For degradation index ####
biom_disp_sp_deg <- betadisper(biom_dist_sp, group = env_sp$inf.percent)
boxplot(biom_disp_sp_deg)
# - More dispersion for shell degradation 2
biom_disp_sp_deg

round(100 * biom_disp_sp_deg$eig / sum(biom_disp_sp_deg$eig), 2)

plot(biom_disp_sp_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = T, cex = 1.5, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)




# For position ####
biom_disp_sp_pos <- betadisper(biom_dist_sp, group = env_sp$position)
boxplot(biom_disp_sp_pos)
# - Less dispersion for infested portion of the shell.
biom_disp_sp_pos

round(100 * biom_disp_sp_pos$eig / sum(biom_disp_sp_pos$eig), 2)

plot(biom_disp_sp_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position2, seg.lwd = 1,
     pch = c(16, 16, 16),
     col = colors_position2, label = T, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)





## SIMPER analyses ===================================================================================

# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)

# SIMPER on shell length ####
biom_simp_length1 <- simper(biom_sp^0.25, group = env_sp$length_range)
biom_simp_length1
summary(biom_simp_length1)

summary_biom_simp_length1 <- summary(biom_simp_length1)
View(summary_biom_simp_length1$'90_80')
View(summary_biom_simp_length1$'90_70')
View(summary_biom_simp_length1$'90_60')
View(summary_biom_simp_length1$'90_50')
View(summary_biom_simp_length1$'90_40')
View(summary_biom_simp_length1$'90_110')
View(summary_biom_simp_length1$'90_100')
View(summary_biom_simp_length1$'90_120')
View(summary_biom_simp_length1$'80_70')
View(summary_biom_simp_length1$'80_60')
View(summary_biom_simp_length1$'80_50')
View(summary_biom_simp_length1$'80_40')
View(summary_biom_simp_length1$'80_110')
View(summary_biom_simp_length1$'80_100')
View(summary_biom_simp_length1$'80_120')
View(summary_biom_simp_length1$'70_60')
View(summary_biom_simp_length1$'70_50')
View(summary_biom_simp_length1$'70_40')
View(summary_biom_simp_length1$'70_110')
View(summary_biom_simp_length1$'70_100')
View(summary_biom_simp_length1$'70_120')
View(summary_biom_simp_length1$'60_50')
View(summary_biom_simp_length1$'60_40')
View(summary_biom_simp_length1$'60_110')
View(summary_biom_simp_length1$'60_100')
View(summary_biom_simp_length1$'60_120')
View(summary_biom_simp_length1$'50_40')
View(summary_biom_simp_length1$'50_110')
View(summary_biom_simp_length1$'50_100')
View(summary_biom_simp_length1$'50_120')
View(summary_biom_simp_length1$'40_110')
View(summary_biom_simp_length1$'40_100')
View(summary_biom_simp_length1$'40_120')
View(summary_biom_simp_length1$'110_100')
View(summary_biom_simp_length1$'110_120')
View(summary_biom_simp_length1$'100_120')

# SIMPER on infestation
biom_simp_inf1 <- simper(biom_sp^0.25, group = env_sp$inf.lvl)
biom_simp_inf1
summary(biom_simp_inf1)

# SIMPER on degradation
biom_simp_deg1 <- simper(biom_sp^0.25, group = env_sp$inf.percent)
biom_simp_deg1
summary(biom_simp_deg1)

# SIMPER on degradation
biom_simp_pos1 <- simper(biom_sp^0.25, group = env_sp$position)
biom_simp_pos1
summary(biom_simp_pos1)
