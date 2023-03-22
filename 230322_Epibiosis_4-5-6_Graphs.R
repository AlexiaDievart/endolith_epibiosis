## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##              Script 4-5-6 - Graphs on community analyses
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
# Section: Load and curate data ---------------------------------------------------------------
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



# Summarize abundance of significant epibiotic higher groups =====================================
View(abund_long)

# Abundance of higher groups per shell length ####






###############################################################################################
# Section 1: Abundance of epibiotic communities -----------------------------------------------
###############################################################################################

# Figure 5.6 - Abundance vs Shell length ======================================================
figure4 <- ggarrange(barplot_abund_group_length2, barplot_abund_sp_length2,
                     nrow = 2, common.legend = F,
                     labels = c("A", "B"), 
                     hjust = c(-1.5, -1.5))
figure4


# BARPLOT - Abundance x Shell length for higher groups ####
abund_group_length <- abund_long[abund_long$higher.group %in% c("Bryozoa", "Cirripedia", "Ochrophyta", "Sedentaria",
                                                                "Bivalvia", "Rhodophyta"),]

abund_group_length1 <- abund_group_length %>%
  dplyr::group_by(quadrat, specimen, length_range, higher.group) %>%
  dplyr::summarise(
    tot_ab = sum(abundance)
  )
View(abund_group_length1)

abund_group_length2 <- abund_group_length %>%
  dplyr::group_by(length_range, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_length2)

barplot_abund_group_length2_text <- data.frame(
  label = c("A", "AB", "AB", "A", "ABC", "BC", "C", "C", "ABC", # Cirripedia
            "AB", "A", "BC", "BCD", "BCD", "CD", "CD", "D", "ABCD"), # Sedentaria
  higher.group = c("Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia",
                   "Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia",
                   "Sedentaria", "Sedentaria", "Sedentaria", "Sedentaria", "Sedentaria",
                   "Sedentaria", "Sedentaria", "Sedentaria", "Sedentaria"),
  length_range = c("40", "50", "60", "70", "80", "90", "100", "110", "120"),
  y = c(6.5, 8.5, 8.5, 6, 7, 8, 9.5, 7.5, 9, # Cirripedia
        3, 3, 6, 11.5, 13.5, 15.5, 16, 22.5, 16) # Sedentaria
)

barplot_abund_group_length2 <- ggplot(abund_group_length2, aes(x = length_range, y = mean_ab, 
                                                                color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_group_length2_text, aes(x = length_range, y = y, label = label), 
            color = "black", size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 23), breaks = seq(from = 0, to = 20, by = 5)) +
  facet_wrap(~ higher.group, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total abundance per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_abund_group_length2

# BARPLOT - Abundance x Shell length for epibiotic species ####
abund_sp_length <- abund_long[abund_long$epibiont %in% c("Tetraclita serrata", "Juvenile barnacle",
                                                         "Ralfsia verrucosa", "Spirorbis spp.",
                                                         "Corallinales"),]
View(abund_sp_length)

abund_sp_length1 <- abund_sp_length %>%
  dplyr::group_by(quadrat, specimen, length_range, epibiont) %>%
  dplyr::summarise(
    tot_ab = sum(abundance)
  )
View(abund_sp_length1)

abund_sp_length2 <- abund_sp_length %>%
  dplyr::group_by(length_range, epibiont) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_sp_length2)

barplot_abund_sp_length2_text <- data.frame(
  label = c("AB", "ABC", "AB", "A", "ABC", "BD", "CD", "ABC", "ABC", # Juvenile barnacle
            "AC", "A", "AC", "BC", "BC", "BC", "B", "B", "ABC"), # Spirorbis spp.
  epibiont = c("Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle",
               "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle",
               "Juvenile barnacle",
               "Spirorbis spp.", "Spirorbis spp.", "Spirorbis spp.", "Spirorbis spp.", "Spirorbis spp.",
               "Spirorbis spp.", "Spirorbis spp.", "Spirorbis spp.", "Spirorbis spp."),
  length_range = c("40", "50", "60", "70", "80", "90", "100", "110", "120"),
  y = c(7, 9.5, 8.5, 7.5, 8.5, 9.5, 12, 9, 11.5, # Juvenile barnacle
        3.5, 3.5, 7, 14, 17.5, 20, 21, 11, 13) 
)

barplot_abund_sp_length2 <- ggplot(abund_sp_length2, aes(x = length_range, y = mean_ab, 
                                                               color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_sp_length2_text, aes(x = length_range, y = y, label = label), 
            color = "black", size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 22), breaks = seq(from = 0, to = 20, by = 5)) +
  facet_wrap(~ epibiont, ncol = 5) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range (mm)", y = "Mean total abundance per mussel") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))
barplot_abund_sp_length2





# Figure 5.7 - Abundance vs (Infestation + Degradation) ==============================================
figure5 <- ggarrange(ggarrange(barplot_abund_group_inf3, barplot_abund_group_inf4, barplot_abund_group_inf5,
                               ncol = 3, common.legend = F, widths = c(1, 1, 1),
                               labels = c("A", "B", "C"), 
                               hjust = c(-1.5, -2.5, -2.5)),
                     ggarrange(barplot_abund_sp_inf3, barplot_abund_sp_deg3,
                               ncol = 2, common.legend = F, widths = c(1,1),
                               labels = c("E", "F"),
                               hjust = c(-1.5, -3)), 
                     nrow = 2)

figure5

# BARPLOT - Abundance x (Infestation + Degradation) for higher groups ####
abund_group_inf <- abund_long[abund_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Rhodophyta"),]

abund_group_inf1 <- abund_group_inf %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent, higher.group) %>%
  dplyr::summarise(
    tot_ab = sum(abundance)
  )
View(abund_group_inf1)

abund_group_inf2 <- abund_group_inf %>%
  dplyr::group_by(inf.lvl, inf.percent, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_inf2)

abund_group_inf2 <- abund_group_inf %>%
  dplyr::group_by(inf.lvl, inf.percent, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_inf2)

# BARPLOT - Abundance x Infestation for Cirripedia and Ochrophyta ####
abund_group_inf3 <- abund_group_inf %>%
  dplyr::filter(higher.group %in% c("Cirripedia", "Ochrophyta")) %>%
  dplyr::group_by(inf.lvl, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_inf3)

barplot_abund_group_inf3_text <- data.frame(
  label = c("A", "AB", "AB", "AB", "B"), # Ochrophyta 
  higher.group = c("Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta"),
  inf.lvl = c("A", "B", "C", "D", "E"),
  y = c(2.5, 3.5, 3.5, 3, 4.5) 
)

barplot_abund_group_inf3 <- ggplot(abund_group_inf3, aes(x = inf.lvl, y = mean_ab,
                                               color = inf.lvl, fill = inf.lvl)) +
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_group_inf3_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3.5) +
  facet_wrap(~ higher.group, ncol = 2) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total abundance per mussel valve") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = c(0.52, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 11, face = "bold")) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_abund_group_inf3

# BARPLOT - Abundance x Degradation for Cirripedia and Ochrophyta ####
abund_group_inf4 <- abund_group_inf %>%
  dplyr::filter(higher.group %in% c("Cirripedia", "Ochrophyta")) %>%
  dplyr::group_by(inf.percent, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_inf4)

barplot_abund_group_inf4_text <- data.frame(
  label = c("AB", "A", "AB", "AB", "B", "AB"), # Ochrophyta 
  higher.group = c("Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta"),
  inf.percent = c("1", "2", "3", "4", "5", "6"),
  y = c(1.5, 3.5, 3.5, 3.8, 2.5, 3)  # Ochrophyta
)

barplot_abund_group_inf4 <- ggplot(abund_group_inf4, aes(x = inf.percent, y = mean_ab,
                                                         color = inf.percent, fill = inf.percent)) +
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_group_inf4_text, aes(x = inf.percent, y = y, label = label), 
            color = "black", size = 3.5) +
  facet_wrap(~ higher.group, ncol = 2) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean total abundance per mussel valve") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = c(0.52, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 11, face = "bold")) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_abund_group_inf4

# BARPLOT - Abundance x (Infestation x Degradation) for Rhodophyta ####
abund_group_inf5 <- abund_group_inf %>%
  dplyr::filter(higher.group %in% c("Rhodophyta")) %>%
  dplyr::group_by(inf.lvl, inf.percent, higher.group) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_group_inf5)

library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill = colors_degrad))

barplot_abund_group_inf5_text <- data.frame(
  label = c("", "A", "A", "A", "B"),
  inf.lvl = c("A", "B", "C", "D", "E"),
  inf.percent = c("4"),
  y = c(2, 4, 3, 6, 24)
)

barplot_abund_group_inf5 <- ggplot(abund_group_inf5, aes(x = inf.lvl, y = mean_ab, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_group_inf5_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3.5) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  scale_fill_manual(values = colors_inflvl) +
  scale_color_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(from = 0, to = 25, by = 5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Rhodophyta", color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean total abundance per mussel valve") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"))
barplot_abund_group_inf5





# BARPLOT - Abundance x (Infestation + Degradation) for epibiotic species ####
abund_sp_inf<- abund_long[abund_long$epibiont %in% c("Tetraclita serrata", "Juvenile barnacle",
                                                         "Ralfsia verrucosa", "Hildenbrandia lecannellieri",
                                                         "Corallinales"),]
View(abund_sp_inf)

abund_sp_inf1 <- abund_sp_inf %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent, epibiont) %>%
  dplyr::summarise(
    tot_ab = sum(abundance)
  )
View(abund_sp_inf1)

# BARPLOT - Abundance x Infestation for epibiotic species ####
abund_sp_inf3 <- abund_sp_inf %>%
  dplyr::group_by(inf.lvl, epibiont) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_sp_inf3)

barplot_abund_sp_inf3_text <- data.frame(
  label = c("A", "AB", "AB", "B", "B"), # RALVER
  epibiont = c("Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa",
               "Ralfsia verrucosa"),
  inf.lvl = c("A", "B", "C", "D", "E"),
  y = c(2.5, 3, 3, 3, 4) 
)

barplot_abund_sp_inf3 <- ggplot(abund_sp_inf3, aes(x = inf.lvl, y = mean_ab,
                                                         color = inf.lvl, fill = inf.lvl)) +
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_sp_inf3_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3.5) +
  facet_wrap(~ epibiont, ncol = 5) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean total abundance per mussel valve") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"))
barplot_abund_sp_inf3

# BARPLOT - Abundance x Degradation for epibiotic species ####
abund_sp_deg3 <- abund_sp_inf %>%
  dplyr::group_by(inf.percent, epibiont) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_sp_deg3)

barplot_abund_sp_deg3_text <- data.frame(
  label = c("AB", "A", "AB", "AB", "B", "AB"), # RALVER
  epibiont = c("Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa",
               "Ralfsia verrucosa", "Ralfsia verrucosa"),
  inf.percent = c("1", "2", "3", "4", "5", "6"),
  y = c(2, 3.5, 3.5, 3.5, 2.2, 3)
)

barplot_abund_sp_deg3 <- ggplot(abund_sp_deg3, aes(x = inf.percent, y = mean_ab,
                                                   color = inf.percent, fill = inf.percent)) +
  geom_errorbar(aes(ymin = mean_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_abund_sp_deg3_text, aes(x = inf.percent, y = y, label = label), 
            color = "black", size = 3.5) +
  facet_wrap(~ epibiont, ncol = 5) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean total abundance per mussel valve") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"))
barplot_abund_sp_deg3



















# Barplot of abundance vs position for each higher group ####
abund_group3 <- abund_long %>%
  dplyr::group_by(position, higher.group) %>%
  dplyr::summarise(
    mean_length = mean(shell.length),
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_group3)

barplot_abund_group3 <- ggplot(abund_group3, aes(x = position, y = mean_ab,
                                                 color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap2(~ higher.group, ncol = 3) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 12, face = "bold"))
barplot_abund_group3



# Summarize abundance of significant epibiotic species =====================================
abund_long1 <- epibio[,1:16]
abund_long1 <- abund_long1[abund_long1$epibiont %in% c("Chthamalus dentatus", "Ralfsia verrucosa",
                                                      "Tetraclita serrata", "Hildenbrandia lecannellieri",
                                                      "Corallinales"),]
View(abund_long1)

abund_juvbar <- abund_long1[abund_long1$epibiont %in% "Juvenile barnacle",]

abund_sp <- abund_long1 %>%
  dplyr::group_by(inf.lvl, inf.percent, epibiont) %>%
  dplyr::summarise(
    mean_length = mean(shell.length),
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_sp)

# Because each species drives community differences in a certain way, there will be several plots

# Barplot of CHTDEN vs position on the shell ####
abund_chtden <- abund_long1 %>%
  dplyr::group_by(position, epibiont) %>%
  dplyr::filter(epibiont == "Chthamalus dentatus") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_chtden)

barplot_abund_chtden <- ggplot(abund_chtden, aes(x = position, y = mean_ab,
                                                 color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 1.4,
           y = 5,
           label = "Chthamalus dentatus", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0,5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0, 0.92),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_chtden

abund_juvbar <- abund_long1 %>%
  dplyr::group_by(position, epibiont) %>%
  dplyr::filter(epibiont == "Juvenile barnacle") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_juvbar)

barplot_abund_juvbar <- ggplot(abund_juvbar, aes(x = position, y = mean_ab,
                                                 color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  scale_y_continuous(limits = c(0,6), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0, 0.92),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_juvbar

abund_chtden1 <- abund_long1 %>%
  dplyr::group_by(inf.percent, epibiont) %>%
  dplyr::filter(epibiont == "Chthamalus dentatus") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_chtden1)

abund_chtden1[nrow(abund_chtden1) + 1,] <- list("1", "Chthamalus dentatus", 0, 0, 0, 0, 0)

barplot_abund_chtden1 <- ggplot(abund_chtden1, aes(x = inf.percent, y = mean_ab,
                                                 color = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 1.4,
           y = 5,
           label = "Chthamalus dentatus", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0,5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0, 0.92),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_chtden1

abund_chtden2 <- abund_long1 %>%
  dplyr::group_by(inf.lvl, epibiont) %>%
  dplyr::filter(epibiont == "Chthamalus dentatus") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_chtden2)

abund_chtden2[nrow(abund_chtden2) + 1,] <- list("A", "Chthamalus dentatus", 0, 0, 0, 0, 0)

barplot_abund_chtden2 <- ggplot(abund_chtden2, aes(x = inf.lvl, y = mean_ab,
                                                   color = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 1.4,
           y = 5,
           label = "Chthamalus dentatus", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0,5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0, 0.92),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_chtden2

# Barplot of TETSER and corallinales vs infestation levels ####
abund_tetser_coral <- abund_long1 %>%
  dplyr::group_by(inf.lvl, epibiont) %>%
  dplyr::filter(epibiont %in% c("Tetraclita serrata", "Corallinales")) %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_tetser_coral)

barplot_abund_tetser_coral <- ggplot(abund_tetser_coral, aes(x = inf.lvl, y = mean_ab,
                                                 color = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap2(~ epibiont, ncol = 2) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 10, face = "bold"))
barplot_abund_tetser_coral

abund_tetser_coral1 <- abund_long1 %>%
  dplyr::group_by(inf.percent, epibiont) %>%
  dplyr::filter(epibiont %in% c("Tetraclita serrata", "Corallinales")) %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_tetser_coral1)

barplot_abund_tetser_coral1 <- ggplot(abund_tetser_coral1, aes(x = inf.percent, y = mean_ab,
                                                             color = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap2(~ epibiont, ncol = 2) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 10, face = "bold"))
barplot_abund_tetser_coral1

abund_tetser_coral2 <- abund_long1 %>%
  dplyr::group_by(position, epibiont) %>%
  dplyr::filter(epibiont %in% c("Tetraclita serrata", "Corallinales")) %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_tetser_coral2)

barplot_abund_tetser_coral2 <- ggplot(abund_tetser_coral2, aes(x = position, y = mean_ab,
                                                               color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  facet_wrap2(~ epibiont, ncol = 2) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'),
        strip.text = element_text(size = 10, face = "bold"))
barplot_abund_tetser_coral2


# Barplot of HILLEC vs shell degradation indexes ####
abund_hillec <- abund_long1 %>%
  dplyr::group_by(inf.percent, epibiont) %>%
  dplyr::filter(epibiont %in% "Hildenbrandia lecannellieri") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_hillec)
glimpse(abund_hillec)

abund_hillec[nrow(abund_hillec) + 1,] <- list("1", "Hildenbrandia lecannellieri", 0, 0, 0, 0, 0)
abund_hillec[nrow(abund_hillec) + 1,] <- list("6", "Hildenbrandia lecannellieri", 0, 0, 0, 0, 0)

barplot_abund_hillec <- ggplot(abund_hillec, aes(x = inf.percent, y = mean_ab,
                                                             color = inf.percent, fill = inf.percent)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 3,
           y = 5,
           label = "Hildenbrandia lecannellieri", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  labs(color = "Shell degradation indexes", fill = "Shell degradation indexes",
       x = "Shell degradation indexes", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_hillec

abund_hillec1 <- abund_long1 %>%
  dplyr::group_by(inf.lvl, epibiont) %>%
  dplyr::filter(epibiont %in% "Hildenbrandia lecannellieri") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_hillec1)

barplot_abund_hillec1 <- ggplot(abund_hillec1, aes(x = inf.lvl, y = mean_ab,
                                                 color = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 3,
           y = 5,
           label = "Hildenbrandia lecannellieri", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_hillec1

abund_hillec2 <- abund_long1 %>%
  dplyr::group_by(position, epibiont) %>%
  dplyr::filter(epibiont %in% "Hildenbrandia lecannellieri") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_hillec2)

barplot_abund_hillec2 <- ggplot(abund_hillec2, aes(x = position, y = mean_ab,
                                                   color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = 3,
           y = 5,
           label = "Hildenbrandia lecannellieri", size = 4, fontface = "bold.italic") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(from = 0, to = 4, by = 1)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_hillec2

# Barplot of RALVER vs shell degradation indexes and infestation levels ####
abund_ralver <- abund_long1 %>%
  dplyr::group_by(inf.percent, inf.lvl, epibiont) %>%
  dplyr::filter(epibiont %in% "Ralfsia verrucosa") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_ralver)

abund_ralver[nrow(abund_ralver) + 1,] <- list("1", "A", "Ralfsia verrucosa", 0, 0, 0, 0, 0)

barplot_abund_ralver <- ggplot(abund_ralver, aes(x = inf.lvl, y = mean_ab,
                                                 color = inf.lvl, fill = inf.lvl)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_text(data = data.frame(x = 2.9, y = 6.8, inf.percent = "1", inf.lvl = "A", label = "Ralfsia \n verrucosa"), 
            aes(x = x, y = y, label = label), size = 4, fontface = "bold.italic", color = "black") +
  facet_wrap2(~ inf.percent, strip = strip, ncol = 6) +
  scale_y_continuous(limits = c(0, 7), breaks = seq(from = 0, to = 7, by = 1)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_ralver

abund_ralver1 <- abund_long1 %>%
  dplyr::group_by(position, epibiont) %>%
  dplyr::filter(epibiont %in% "Ralfsia verrucosa") %>%
  dplyr::summarise(
    tot_ab = sum(abundance),
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab = sd_ab / sqrt(n)
  )
View(abund_ralver1)

barplot_abund_ralver1 <- ggplot(abund_ralver1, aes(x = position, y = mean_ab,
                                                 color = position, fill = position)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  scale_y_continuous(limits = c(0, 7), breaks = seq(from = 0, to = 7, by = 1)) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "Mean abundance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent'))
barplot_abund_ralver1





###############################################################################################
# Section: Epibiotic community percent cover --------------------------------------------------
###############################################################################################

# Figure 5.8 - nMDS on percent cover (higher groups and species) =====================================
par(mfrow = c(1,2), oma=c(1,1,1,1))





# Plot nMDS for higher groups w/ infestation levels and degradation ####
ordiplot(cover_bray_group, type = "none", xlim = c(-0.1, 5), ylim = c(-0.5, 0.5))
points(cover_bray_group, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
orditorp(cover_bray_group, display = "species", col = "black", air = 0.4, cex = 0.8)
legend(-0.6, - 1.1, title = "SD index",
       legend = levels(env_group$inf.percent),
       pch = c(16, 17, 18, 8, 5, 3), 
       col = c("black", 'black', "black", "black", "black", "black"),
       cex = 0.7, box.lty = 0, bg = "transparent", pt.cex = 1.5,
       y.intersp = 0.5)
legend(-0.5, - 1.1, title = "Euendolithic infestation",
       legend = levels(env_group$inf.lvl),
       pch = 16,
       col = c("coral4", "brown3", "indianred1", "pink3", "darkgrey"),
       cex = 0.7, box.lty = 0, bg = "transparent", pt.cex = 1.5,
       y.intersp = 0.5)
legend(-0.95, 3, "Stress = 0.0735", bty = "n", cex = 1)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

# Plot nMDS for epibiotic species w/ infestation levels and degradation ====================================
ordiplot(cover_bray_sp, type = "none", xlim = c(0, 1), ylim = c(-1.0, 1.0))
points(cover_bray_sp, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
orditorp(cover_bray_sp, display = "species", col = "black", air = 0.4, cex = 0.8)
legend(-0.85, 1.2, "Stress = 0.003", bty = "n", cex = 1)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")





# Figure 5.9 - Dispersion plots on percent cover (higher groups and species) =================================
par(mfrow = c(2,4), oma = c(1,1,1,1))

# Higher group - Dispersion for shell length
plot(cover_disp_group_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = F, cex = 0.2, lwd = 2,
     xlim = c(-0.5, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
ordilabel(scores(cover_disp_group_length, "centroids"), col = "black", border = col_length, fill = "white",
          cex = 0.7)

mtext("A", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Higher group - Dispersion for infestation levels ####
plot(cover_disp_group_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = F, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
ordilabel(scores(cover_disp_group_inf, "centroids"), col = "black", border = colors_inflvl, fill = "white",
          cex = 0.7)

mtext("B", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Higher group - Dispersion for degradation ####
plot(cover_disp_group_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16),
     col = colors_degrad, label = F, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
ordilabel(scores(cover_disp_group_deg, "centroids"), col = "black", border = colors_degrad, fill = "white",
          cex = 0.7)

mtext("C", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Higher group - Dispersion for position ####
plot(cover_disp_group_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position1, seg.lwd = 1,
     pch = c(16, 16),
     col = colors_position1, label = F, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
ordilabel(scores(cover_disp_group_pos, "centroids"), col = "black", border = colors_position1, fill = "white",
          cex = 0.7)

mtext("D", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Epibiotic species - Dispersion for shell length ####
plot(cover_disp_sp_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = F, cex = 0.5, lwd = 2,
     xlim = c(-1, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
ordilabel(scores(cover_disp_sp_length, "centroids"), col = "black", border = col_length, fill = "white",
          cex = 0.7)

mtext("E", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Epibiotic species - Dispersion for infestation levels ####
plot(cover_disp_sp_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = F, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.81 %)", ylab = "Dimension 2 (17.76 %)",
     main = ""
)
ordilabel(scores(cover_disp_sp_inf, "centroids"), col = "black", border = colors_inflvl, fill = "white",
          cex = 0.7)

mtext("F", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Epibiotic species - Dispersion for degradation index ####
plot(cover_disp_sp_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16),
     col = colors_degrad, label = F, cex = 1, lwd = 2,
     xlim = c(-0.9, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
ordilabel(scores(cover_disp_sp_deg, "centroids"), col = "black", border = colors_degrad, fill = "white",
          cex = 0.7)

mtext("G", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")

# Epibiotic species - Dispersion for position ####
plot(cover_disp_sp_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position1, seg.lwd = 1,
     pch = c(16, 16),
     col = colors_position1, label = F, cex = 1, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (19.16 %)", ylab = "Dimension 2 (13.88 %)",
     main = ""
)
ordilabel(scores(cover_disp_sp_pos, "centroids"), col = "black", border = colors_position1, fill = "white",
          cex = 0.7)

mtext("H", side = 3, line = -1, cex = 1.5, adj = -0.17, col = "grey30")





# Figure 5.10 - Percent cover vs Shell length =================================
figure10 <- ggarrange(barplot_cover_group_length2, barplot_cover_sp_length2,
                     nrow = 2, common.legend = F,
                     labels = c("A", "B"), 
                     hjust = c(-1.5, -1.5))
figure10








# BARPLOT - Percent cover x Shell length for higher groups ####

# Summarize percent cover of significant epibiotic higher groups ####
cover_long <- epibio[,c(1:15, 18,22)]
View(cover_long)

cover_group_length <- cover_long[cover_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Sedentaria",
                                                                "Rhodophyta"),]
View(cover_group_length)

cover_group_length1 <- cover_group_length %>%
  dplyr::group_by(quadrat, specimen, length_range, higher.group) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_group_length1)

cover_group_length2 <- cover_group_length %>%
  dplyr::group_by(length_range, higher.group) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_length2)

barplot_cover_group_length2_text <- data.frame(
  label = c("A", "AB", "AB", "AB", "AB", "B", "AB", "AB", ""), 
  higher.group = c("Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta",
                   "Ochrophyta", "Ochrophyta", "Ochrophyta", "Ochrophyta"),
  length_range = c("40", "50", "60", "70", "80", "90", "100", "110", "120"),
  y = c(30, 28, 26, 24, 30, 23.5, 27, 12, 1)
)

barplot_cover_group_length2 <- ggplot(cover_group_length2, aes(x = length_range, y = mean_cover, 
                                                               color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_group_length2_text, aes(x = length_range, y = y, label = label), 
            color = "black", size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 32), breaks = seq(from = 0, to = 30, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ higher.group, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_group_length2





# BARPLOT - Percent cover x Shell length for epibiotic species ####

# Summarize percent cover of significant epibiotic species
cover_long <- epibio[,c(1:15, 18,22)]
View(cover_long)

cover_sp <- cover_long[cover_long$epibiont %in% c("Spirobranchus kraussii", "Ralfsia verrucosa",
                                                        "Corallinales", "Tetraclita serrata",
                                                        "Hildenbrandia lecannellieri", "Juvenile barnacle",
                                                        "Spirorbis spp.", "Chthamalus dentatus"),]
View(cover_sp)

cover_sp_length1 <- cover_sp %>%
  dplyr::group_by(quadrat, specimen, length_range, epibiont) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_sp_length1)

cover_sp_length2 <- cover_sp %>%
  dplyr::group_by(length_range, epibiont) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_length2)

barplot_cover_sp_length2_text <- data.frame(
  label = c("AB", "A", "A", "A", "A", "A", "A", "B", "AB", # Juvenile barnacle
            "A", "AB", "AB", "AB", "AB", "B", "AB", "AB", "", # Ralfsia verrucosa
            "AB", "AB", "AB", "AB", "A", "AB", "B", "AB", "AB"), 
  epibiont = c("Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle",
               "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle", "Juvenile barnacle",
               "Juvenile barnacle",
               "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa",
               "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa", "Ralfsia verrucosa",
               "Ralfsia verrucosa",
               "Corallinales", "Corallinales", "Corallinales", "Corallinales", "Corallinales",
               "Corallinales", "Corallinales", "Corallinales", "Corallinales"),
  length_range = c("40", "50", "60", "70", "80", "90", "100", "110", "120"),
  y = c(2, 2, 2, 2, 2, 3, 2, 17, 3, # Juvenile barnacle
        30, 46, 43, 46, 32, 25, 28, 13, 1, # Ralfsia verrucosa
        20, 10, 4, 8, 12, 10, 7, 7, 5) # Corallinales
)

barplot_cover_sp_length2 <- ggplot(cover_sp_length2, aes(x = length_range, y = mean_cover, 
                                                         color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_sp_length2_text, aes(x = length_range, y = y, label = label), 
            color = "black", size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 46), breaks = seq(from = 0, to = 50, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ epibiont, ncol = 8) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range (mm)", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))
barplot_cover_sp_length2





# Figure 5.11 - Percent cover vs (Infestation + Degradation) =================================
figure11 <- ggarrange(ggarrange(barplot_cover_group_tot1, barplot_cover_group_sep1, barplot_cover_group_sep2,
                                nrow = 1, ncol = 3, common.legend = F,
                                labels = c("A", "B", "C"),
                                hjust = c(-2, -2, -2)),
                      ggarrange(barplot_cover_group_ochro_int1, barplot_cover_group_seden_int1,
                                nrow = 1, ncol = 2, common.legend = F,
                                labels = c("D", "E"),
                                hjust = c(-2, -3)),
                      nrow = 2
                      )
figure11





# BARPLOT - Percent cover x Infestation + Degradation for higher groups ####

# Summarize percent cover of significant epibiotic higher groups
cover_long <- epibio[,c(1:15, 18,22)]
View(cover_long)

cover_long <- cover_long[cover_long$position %in% c("infested", "non infested"),]

# BARPLOT - Infestation ####
cover_group_tot <- cover_long[cover_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Sedentaria",
                                                             "Rhodophyta"),]
View(cover_group_tot)

cover_group_tot1 <- cover_group_tot %>%
  dplyr::group_by(higher.group, inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_tot1)

barplot_cover_group_tot1 <- ggplot(cover_group_tot1, aes(x = inf.lvl, y = mean_cover, 
                                                         color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(from = 0, to = 50, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ higher.group, ncol = 4) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_group_tot1

# BARPLOT - Degradation for Cirripedia and Rhodophyta ####
cover_group_sep <- cover_long[cover_long$higher.group %in% c("Cirripedia", "Rhodophyta"),]
View(cover_group_sep)

cover_group_sep1 <- cover_group_sep %>%
  dplyr::group_by(higher.group, inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_sep1)

barplot_cover_group_sep1_text <- data.frame(
  label = c("AB", "AB", "A", "AB", "B", "AB"), 
  higher.group = c("Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia", "Cirripedia"),
  inf.percent = c("1", "2", "3", "4", "5", "6"),
  y = c(2, 5, 19, 12, 10, 14)
)

barplot_cover_group_sep1 <- ggplot(cover_group_sep1, aes(x = inf.percent, y = mean_cover, 
                                                         color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_group_sep1_text, aes(x = inf.percent, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(from = 0, to = 50, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ higher.group, ncol = 4) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_group_sep1





# BARPLOT - Position for Cirripedia and Rhodophyta ####
cover_group_sep <- cover_long[cover_long$higher.group %in% c("Cirripedia", "Rhodophyta"),]
View(cover_group_sep)

cover_group_sep2 <- cover_group_sep %>%
  dplyr::group_by(higher.group, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_sep2)

barplot_cover_group_sep2_text <- data.frame(
  label = c("A", "B"), 
  higher.group = c("Rhodophyta", "Rhodophyta"),
  position = c("infested", "non infested"),
  y = c(14, 29)
)

barplot_cover_group_sep2 <- ggplot(cover_group_sep2, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_group_sep2_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(from = 0, to = 50, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ higher.group, ncol = 4) +
  labs(color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_group_sep2





# BARPLOT - Position x Degradation for Ochrophyta ####
cover_group_ochro_int <- cover_long[cover_long$higher.group %in% "Ochrophyta",]
View(cover_group_ochro_int)

cover_group_ochro_int1 <- cover_group_ochro_int %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_ochro_int1)

barplot_cover_group_ochro_int1_text <- data.frame(
  label = c("A", "B"), 
  higher.group = c("Ochrophyta", "Ochrophyta"),
  inf.percent = "4",
  position = c("infested", "non infested"),
  y = c(8, 35)
)

library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill = colors_degrad))

barplot_cover_group_ochro_int1 <- ggplot(cover_group_ochro_int1, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_group_ochro_int1_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 55), breaks = seq(from = 0, to = 60, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Ochrophyta", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
barplot_cover_group_ochro_int1





# BARPLOT - Position x Degradation for Sedentaria ####
cover_group_seden_int <- cover_long[cover_long$higher.group %in% "Sedentaria",]
View(cover_group_seden_int)

cover_group_seden_int1 <- cover_group_seden_int %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_group_seden_int1)

barplot_cover_group_seden_int1 <- ggplot(cover_group_seden_int1, aes(x = position, y = mean_cover, 
                                                                     color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(from = 0, to = 60, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Sedentaria", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
barplot_cover_group_seden_int1





# Figure 5.12 - Percent cover vs (Infestation + Degradation) =================================
figure12 <- ggarrange(ggarrange(barplot_cover_sp_inf1, barplot_cover_sp_deg1, barplot_cover_sp_pos1,
                                nrow = 1, ncol = 3, common.legend = F,
                                labels = c("A", "B", "C"),
                                hjust = c(-2, -2, -2)),
                      ggarrange(barplot_cover_sp_juvbar1, barplot_cover_sp_ralver1,
                                nrow = 1, ncol = 2, common.legend = F,
                                labels = c("D", "E"),
                                hjust = c(-2, -3)),
                      ggarrange(barplot_cover_sp_hillec1, barplot_cover_sp_coral1,
                                nrow = 1, ncol = 2, common.legend = F,
                                labels = c("F", "G"),
                                hjust = c(-2, -3)),
                      nrow = 3
)
figure12





# BARPLOT - Infestation ####
cover_sp_inf <- cover_long[cover_long$epibiont %in% c("Tetraclita serrata", "Hildenbrandia lecannellieri",
                                                      "Corallinales", "Spirobranchus kraussii",
                                                      "Spirorbis spp."),]
View(cover_sp_inf)

cover_sp_inf1 <- cover_sp_inf %>%
  dplyr::group_by(epibiont, inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_inf1)


barplot_cover_sp_inf1 <- ggplot(cover_sp_inf1, aes(x = inf.lvl, y = mean_cover, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(from = 0, to = 40, by = 10),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ epibiont, ncol = 5, labeller = label_wrap_gen(width=10)) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_cover_sp_inf1

# BARPLOT - Degradation ####
cover_sp_deg <- cover_long[cover_long$epibiont %in% c("Tetraclita serrata", "Ralfsia verrucosa",
                                                      "Spirobranchus kraussii", "Spirorbis spp."),]
View(cover_sp_deg)

cover_sp_deg1 <- cover_sp_deg %>%
  dplyr::group_by(epibiont, inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_deg1)


barplot_cover_sp_deg1 <- ggplot(cover_sp_deg1, aes(x = inf.percent, y = mean_cover, 
                                                   color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(from = 0, to = 50, by = 10),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ epibiont, ncol = 4, labeller = label_wrap_gen(width=10)) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_cover_sp_deg1





# BARPLOT - Position ####
cover_sp_pos <- cover_long[cover_long$epibiont %in% c("Tetraclita serrata", "Juvenile barnacle",
                                                      "Spirobranchus kraussii", "Spirorbis spp."),]
View(cover_sp_pos)

cover_sp_pos1 <- cover_sp_pos %>%
  dplyr::group_by(epibiont, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_pos1)

barplot_cover_sp_pos1_text <- data.frame(
  label = c("A", "B"), 
  epibiont = c("Spirobranchus kraussii", "Spirobranchus kraussii"),
  position = c("infested", "non infested"),
  y = c(6, 4)
)


barplot_cover_sp_pos1 <- ggplot(cover_sp_pos1, aes(x = position, y = mean_cover, 
                                                   color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_sp_pos1_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 22), breaks = seq(from = 0, to = 20, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap(~ epibiont, ncol = 4, labeller = label_wrap_gen(width=10)) +
  labs(color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_cover_sp_pos1





# BARPLOT - JUVBAR - Infestation x Degradation ####
cover_sp_juvbar <- cover_long[cover_long$epibiont %in% "Juvenile barnacle",]
View(cover_sp_juvbar)

cover_sp_juvbar1 <- cover_sp_juvbar %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_juvbar1)

barplot_cover_sp_juvbar1_text <- data.frame(
  label = c("", "A", "B", "B", "AB"),
  inf.percent = "5",
  inf.lvl = c("A", "B", "C", "D", "E"),
  y = c(1, 2, 2, 2, 6)
)

barplot_cover_sp_juvbar1 <- ggplot(cover_sp_juvbar1, aes(x = inf.lvl, y = mean_cover, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_sp_juvbar1_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(from = 0, to = 25, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Juvenile barnacle", color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_cover_sp_juvbar1





# BARPLOT - RALVER - Infestation x Position ####
cover_sp_ralver <- cover_long[cover_long$epibiont %in% "Ralfsia verrucosa",]
View(cover_sp_ralver)

cover_sp_ralver1 <- cover_sp_ralver %>%
  dplyr::group_by(inf.lvl, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_ralver1)

barplot_cover_sp_ralver1_text <- data.frame(
  label = c("A", "B"),
  inf.lvl = "B",
  position = c("infested", "non infested"),
  y = c(7, 55)
)

strip1 <- strip_themed(background_x = elem_list_rect(fill = colors_inflvl))

barplot_cover_sp_ralver1 <- ggplot(cover_sp_ralver1, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_sp_ralver1_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 55), breaks = seq(from = 0, to = 55, by = 10),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap2(~ inf.lvl, ncol = 5, strip = strip1) +
  labs(title = "Ralfsia verrucosa", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_cover_sp_ralver1





# BARPLOT - HILLEC - Degradation x Position ####
cover_sp_hillec <- cover_long[cover_long$epibiont %in% "Hildenbrandia lecannellieri",]
View(cover_sp_hillec)

cover_sp_hillec1 <- cover_sp_hillec %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_hillec1)

cover_sp_hillec1[nrow(cover_sp_hillec1) + 1,] <- list("1", "infested", NA, NA)
cover_sp_hillec1[nrow(cover_sp_hillec1) + 1,] <- list("6", "infested", NA, NA)

barplot_cover_sp_hillec1 <- ggplot(cover_sp_hillec1, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(from = 0, to = 45, by = 10),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Hildenbrandia lecannellieri", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_cover_sp_hillec1





# BARPLOT - Corallinales - Degradation x Position ####
cover_sp_coral <- cover_long[cover_long$epibiont %in% "Corallinales",]
View(cover_sp_coral)

cover_sp_coral1 <- cover_sp_coral %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_coral1)

barplot_cover_sp_coral1_text <- data.frame(
  label = c("A", "B", "A", "B"),
  inf.percent = c("2", "2", "5", "5"),
  position = c("infested", "non infested",
               "infested", "non infested"),
  y = c(13, 4, 8, 4)
)

barplot_cover_sp_coral1 <- ggplot(cover_sp_coral1, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover_sp_coral1_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(from = 0, to = 15, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Corallinales", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean percent cover per mussel valve") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_cover_sp_coral1



























# BARPLOT - Percent cover x (Infestation x Degradation x Position) for epibiotic species
cover_group_sp <- cover_long[cover_long$epibiont %in% c("Spirobranchus kraussii", "Ralfsia verrucosa",
                                                        "Corallinales", "Tetraclita serrata",
                                                        "Hildenbrandia lecannellieri", "Juvenile barnacle",
                                                        "Spirorbis spp.", "Chthamalus dentatus"),]
View(cover_group_sp)

cover_juvbar3 <- cover_juvbar %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_juvbar3)

barplot_cover_juvbar_inf <- ggplot(cover_juvbar3, aes(x = inf.lvl, y = mean_cover, 
                                                         color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  #  geom_text(data = barplot_cover_group_length2_text, aes(x = length_range, y = y, label = label), 
  #            color = "black", size = 3) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(from = 0, to = 30, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ inf.percent, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_juvbar_inf

cover_sp_ralver <- cover_long[cover_long$epibiont %in% c("Ralfsia verrucosa"),]
View(cover_sp_ralver)

cover_sp_ralver1 <- cover_sp_ralver %>%
  dplyr::group_by(inf.lvl, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_ralver1)

barplot_cover_ralver_inf <- ggplot(cover_sp_ralver1, aes(x = position, y = mean_cover, 
                                                      color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  #  geom_text(data = barplot_cover_group_length2_text, aes(x = length_range, y = y, label = label), 
  #            color = "black", size = 3) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(from = 0, to = 30, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ inf.lvl, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_ralver_inf

cover_sp_coral <- cover_long[cover_long$epibiont %in% c("Corallinales"),]
View(cover_sp_coral)

cover_sp_coral1 <- cover_sp_coral %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_coral1)

cover_sp_coral2 <- cover_sp_coral %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover.percent),
    sd_cover = sd(cover.percent)
  )
View(cover_sp_coral2)

barplot_cover_coral_inf2 <- ggplot(cover_sp_coral2, aes(x = position, y = mean_cover, 
                                                         color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  #  geom_text(data = barplot_cover_group_length2_text, aes(x = length_range, y = y, label = label), 
  #            color = "black", size = 3) +
  scale_colour_manual(values = colors_position1) +
  scale_fill_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(from = 0, to = 30, by = 5),
                     labels = scales::percent_format(scale = 1)) +
  facet_wrap(~ inf.percent, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_cover_coral_inf2




# Figure 5.13 - nMDS on biomass =================================
par(mfrow = c(1,2), oma=c(1,1,1,1))

# Plot nMDS for higher groups w/ infestation levels and degradation ####
ordiplot(biom_bray_group, type = "none", xlim = c(-3, 2), ylim = c(-1.5, 1.5))
points(biom_bray_group, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_group$inf.percent],
       col = colors_inflvl[env_group$inf.lvl])
orditorp(biom_bray_group, display = "species", col = "black", air = 0.4, cex = 0.8)
legend(-3.5, -1.4, title = "SD index",
       legend = levels(env_group$inf.percent),
       pch = c(16, 17, 18, 8, 5, 3), 
       col = c("black", 'black', "black", "black", "black", "black"),
       cex = 0.7, box.lty = 0, bg = "transparent", pt.cex = 1.5,
       y.intersp = 0.5)
legend(-3.5, -1.4, title = "Euendolithic infestation",
       legend = levels(env_group$inf.lvl),
       pch = 16,
       col = c("coral4", "brown3", "indianred1", "pink3", "darkgrey"),
       cex = 0.7, box.lty = 0, bg = "transparent", pt.cex = 1.5,
       y.intersp = 0.5)
legend(-3.85, 3.3, "Stress = 0.0044", bty = "n", cex = 1)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")

# Plot nMDS for epibiotic species w/ infestation levels and degradation ====================================
ordiplot(biom_bray_sp, type = "none", xlim = c(-50, 50), ylim = c(-35, 30))
points(biom_bray_sp, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5, 3)[env_sp$inf.percent],
       col = colors_inflvl[env_sp$inf.lvl])
orditorp(biom_bray_sp, display = "species", col = "black", air = 0.4, cex = 0.8)
legend(-68, 64, "Stress = 0.007", bty = "n", cex = 1)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")





# Figure 5.14 - Dispersion plots on biomass =================================
par(mfrow = c(2,4), oma = c(1,1,1,1))

# Higher group - Dispersion for shell length
plot(biom_disp_group_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = F, cex = 0.5, lwd = 2,
     xlim = c(-0.8, 0.8), ylim = c(-0.7, 0.7),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)
ordilabel(scores(biom_disp_group_length, "centroids"), col = "black", border = col_length, fill = "white",
          cex = 0.7)

mtext("A", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Higher group - Dispersion for infestation levels ####
plot(biom_disp_group_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = F, cex = 1, lwd = 2,
     xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)
ordilabel(scores(biom_disp_group_inf, "centroids"), col = "black", border = colors_inflvl, fill = "white",
          cex = 0.7)

mtext("B", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Higher group - Dispersion for degradation ####
plot(biom_disp_group_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = F, cex = 1.5, lwd = 2,
     xlim = c(-1, 1), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)
ordilabel(scores(biom_disp_group_deg, "centroids"), col = "black", border = colors_degrad, fill = "white",
          cex = 0.7)

mtext("C", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Higher group - Dispersion for position ####
plot(biom_disp_group_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position2, seg.lwd = 1,
     pch = c(16, 16, 16),
     col = colors_position2, label = F, cex = 1, lwd = 2,
     xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (26.07 %)", ylab = "Dimension 2 (21.56 %)",
     main = ""
)
ordilabel(scores(biom_disp_group_pos, "centroids"), col = "black", border = colors_position2, fill = "white",
          cex = 0.7)

mtext("D", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Epibiotic species - Dispersion for shell length ####
plot(biom_disp_sp_length, hull = F, ellipse = T, segments = F, seg.col = col_length, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16, 16, 16, 16, 16),
     col = col_length, label = F, cex = 0.5, lwd = 2,
     xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)
ordilabel(scores(biom_disp_sp_length, "centroids"), col = "black", border = col_length, fill = "white",
          cex = 0.7)

mtext("E", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Epibiotic species - Dispersion for infestation levels ####
plot(biom_disp_sp_inf, hull = F, ellipse = T, segments = F, seg.col = colors_inflvl, seg.lwd = 1,
     pch = c(16, 16, 16, 16, 16),
     col = colors_inflvl, label = F, cex = 1, lwd = 2,
     xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)
ordilabel(scores(biom_disp_sp_inf, "centroids"), col = "black", border = colors_inflvl, fill = "white",
          cex = 0.7)

mtext("F", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Epibiotic species - Dispersion for degradation index ####
plot(biom_disp_sp_deg, hull = F, ellipse = T, segments = F, seg.col = colors_degrad, seg.lwd = 1,
     pch = c(16, 17, 18, 8, 5, 3),
     col = colors_degrad, label = F, cex = 1.5, lwd = 2,
     xlim = c(-0.8, 0.8), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)
ordilabel(scores(biom_disp_sp_deg, "centroids"), col = "black", border = colors_degrad, fill = "white",
          cex = 0.7)

mtext("G", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")

# Epibiotic species - Dispersion for position ####
plot(biom_disp_sp_pos, hull = F, ellipse = T, segments = F, seg.col = colors_position2, seg.lwd = 1,
     pch = c(16, 16, 16),
     col = colors_position2, label = F, cex = 1, lwd = 2,
     xlim = c(-0.7, 0.7), ylim = c(-0.6, 0.6),
     xlab = "Dimension 1 (24.59 %)", ylab = "Dimension 2 (18.11 %)",
     main = ""
)
ordilabel(scores(biom_disp_sp_pos, "centroids"), col = "black", border = colors_position2, fill = "white",
          cex = 0.7)

mtext("H", side = 3, line = -1, cex = 1.5, adj = -0.15, col = "grey30")





# Figure 5.15 - Biomass vs Shell length =================================
figure15 <- ggarrange(barplot_biom_long_length3, barplot_biom_long_sp3,
                      nrow = 2, common.legend = F,
                      labels = c("A", "B"), 
                      hjust = c(-1.5, -1.5))
figure15








# BARPLOT - Biomass x Shell length for higher groups ####

# Summarize biomass of significant epibiotic higher groups ####
biom_long_length <- biom_long[biom_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Rhodophyta",
                                                            "Bryozoa", "Bivalvia", "Sedentaria"),]
View(biom_long_length)

biom_long_length1 <- biom_long[biom_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Rhodophyta",
                                                             "Bryozoa"),]

biom_long_length2 <- biom_long_length1 %>%
  dplyr::group_by(quadrat, specimen, length_range, higher.group) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_long_length2)

biom_long_length3 <- biom_long_length1 %>%
  dplyr::group_by(length_range, higher.group) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_long_length3)

barplot_biom_long_length3 <- ggplot(biom_long_length3, aes(x = length_range, y = mean_biom, 
                                                               color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 1200), breaks = seq(from = 0, to = 1200, by = 200)) +
  facet_wrap(~ higher.group, ncol = 6) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean biomass per mussel (mg)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none", 
        strip.text = element_text(size = 11, face = "bold"))
barplot_biom_long_length3





# BARPLOT - Biomass x Shell length for epibiotic species ####

# Summarize biomass of significant epibiotic species
biom_long_sp <- biom_long[biom_long$epibiont %in% c("Tetraclita serrata", "Chthamalus dentatus",
                                                    "Ralfsia verrucosa", "Hildenbrandia lecannellieri",
                                                    "Electra verticillata", "Bryozoa", "Perna perna",
                                                    "Spirobranchus kraussii", "Gelidium pristoides",
                                                    "Gunnarea gaimardi"),]
View(biom_long_sp)

biom_long_sp1 <- biom_long[biom_long$epibiont %in% c("Tetraclita serrata", "Ralfsia verrucosa", 
                                                     "Hildenbrandia lecannellieri",
                                                     "Electra verticillata", "Bryozoa"),]

biom_long_sp2 <- biom_long_sp1 %>%
  dplyr::group_by(quadrat, specimen, length_range, epibiont) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_long_sp2)

biom_long_sp3 <- biom_long_sp1 %>%
  dplyr::group_by(length_range, epibiont) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_long_sp3)

barplot_biom_long_sp3 <- ggplot(biom_long_sp3, aes(x = length_range, y = mean_biom, 
                                                         color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(from = 0, to = 150, by = 20)) +
  facet_wrap(~ epibiont, ncol = 8) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range (mm)", y = "Mean biomass per mussel (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))
barplot_biom_long_sp3






# Figure 5.16 - Biomass vs (Infestation x Degradation x Position) =================================
figure16 <- ggarrange(barplot_biom_group_inf1, barplot_biom_group_deg,
                      barplot_biom_group_pos1, barplot_biom_group_cirri1,
                      nrow = 2, ncol = 2, common.legend = F,
                      labels = c("A", "B", "C", "D"), 
                      hjust = c(-1.5, -2, -1.5, -2))
figure16





# BARPLOT - Biomass of higher groups vs Infestation (except Cirripedia) ####
biom_group_inf <- biom_long[biom_long$higher.group %in% c("Ochrophyta", "Rhodophyta", "Bryozoa"),]
View(biom_group_inf)

biom_group_inf1 <- biom_group_inf %>%
  dplyr::group_by(higher.group, inf.lvl) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_group_inf1)

barplot_biom_group_inf1_text <- data.frame(
  label = c("A", "AB", "BC", "C", "ABC"),
  higher.group = c("Rhodophyta", "Rhodophyta", "Rhodophyta", "Rhodophyta", "Rhodophyta"),
  inf.lvl = c("A", "B", "C", "D", "E"),
  y = c(2100, 400, 780, 1100, 150)
)

barplot_biom_group_inf1 <- ggplot(biom_group_inf1, aes(x = inf.lvl, y = mean_biom, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_biom_group_inf1_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 2100), breaks = seq(from = 0, to = 2000, by = 250)) +
  facet_wrap(~ higher.group, ncol = 3) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))
barplot_biom_group_inf1






# BARPLOT - Biomass of higher groups vs Degradation (except Cirripedia) ####
biom_group_inf <- biom_long[biom_long$higher.group %in% c("Ochrophyta", "Rhodophyta", "Bryozoa"),]
View(biom_group_inf)

biom_group_deg <- biom_group_inf %>%
  dplyr::group_by(higher.group, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_group_deg)

barplot_biom_group_deg <- ggplot(biom_group_deg, aes(x = inf.percent, y = mean_biom, 
                                                       color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_color_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  scale_y_continuous(limits = c(0, 2050), breaks = seq(from = 0, to = 2000, by = 250)) +
  facet_wrap(~ higher.group, ncol = 3) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))
barplot_biom_group_deg






# BARPLOT - Biomass of higher groups vs Position ####
biom_group_pos <- biom_long[biom_long$higher.group %in% c("Cirripedia", "Ochrophyta", "Rhodophyta", "Bryozoa"),]
View(biom_group_pos)

biom_group_pos1 <- biom_group_pos %>%
  dplyr::group_by(higher.group, position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_group_pos1)

barplot_biom_group_pos1 <- ggplot(biom_group_pos1, aes(x = position, y = mean_biom, 
                                                     color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_color_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 1000), breaks = seq(from = 0, to = 2000, by = 250)) +
  scale_x_discrete(labels = c("II epibio", "Inf", "Non-inf")) +
  facet_wrap(~ higher.group, ncol = 4) +
  labs(color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"))
barplot_biom_group_pos1






# BARPLOT - Biomass of higher groups vs Position ####
biom_group_cirri <- biom_long[biom_long$higher.group %in% "Cirripedia",]
View(biom_group_cirri)

biom_group_cirri1 <- biom_group_cirri %>%
  dplyr::group_by(higher.group, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_group_cirri1)

biom_group_cirri1[nrow(biom_group_cirri1) + 1,] <- list("Cirripedia", "A", "1", NA, NA)

library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill = colors_degrad))

barplot_biom_group_cirri1 <- ggplot(biom_group_cirri1, aes(x = inf.lvl, y = mean_biom, 
                                                       color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 130), breaks = seq(from = 0, to = 25, by = 25)) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = " Cirripedia", color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
barplot_biom_group_cirri1






# Figure 5.17 - Biomass vs (Infestation x Degradation x Position) =================================
figure17 <- ggarrange(ggarrange(barplot_biom_sp_inf1, barplot_biom_sp_deg1, barplot_biom_sp_pos1,
                                nrow = 1, ncol = 3, common.legend = F,
                                labels = c("A", "B", "C"),
                                hjust = c(-2, -2, -2)),
                      ggarrange(barplot_tetser_sp_biom1, barplot_gelpri_sp_biom1,
                                nrow = 1, ncol = 2, common.legend = F,
                                labels = c("D", "E"),
                                hjust = c(-2, -3)),
                      ggarrange(barplot_hillec_sp_biom1, barplot_gelpri_sp_biom2,
                                nrow = 1, ncol = 2, common.legend = F,
                                labels = c("F", "G"),
                                hjust = c(-2, -3)),
                      nrow = 3
)
figure17


# BARPLOT - Species biomass vs Infestation ####
biom_sp_inf <- biom_long[biom_long$epibiont %in% c("Hildenbrandia lecannellieri",
                                                   "Ralfsia verrucosa", "Electra verticillata",
                                                   "Bryozoa"),]

biom_sp_inf1 <- biom_sp_inf %>%
  dplyr::group_by(epibiont, inf.lvl) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_sp_inf1)

barplot_biom_sp_inf1_text <- data.frame(
  label = c("AB", " AB", "A", "B", "AB",
            "", "A", "AB", "B", ""),
  epibiont = c("Hildenbrandia lecannellieri", "Hildenbrandia lecannellieri", "Hildenbrandia lecannellieri",
               "Hildenbrandia lecannellieri", "Hildenbrandia lecannellieri",
               "Electra verticillata", "Electra verticillata", "Electra verticillata", "Electra verticillata",
               "Electra verticillata"),
  inf.lvl = c("A", "B", "C", "D", "E"),
  y = c(40, 150, 220, 90, 40,
        25, 55, 30, 30, 10)
)

barplot_biom_sp_inf1 <- ggplot(biom_sp_inf1, aes(x = inf.lvl, y = mean_biom, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_biom_sp_inf1_text, aes(x = inf.lvl, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 220), breaks = seq(from = 0, to = 200, by = 50)) +
  facet_wrap(~ epibiont, ncol = 4, labeller = label_wrap_gen(width=10)) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_biom_sp_inf1





# BARPLOT - Species biomass vs Degradation ####
biom_sp_deg <- biom_long[biom_long$epibiont %in% c("Ralfsia verrucosa", "Electra verticillata",
                                                   "Bryozoa"),]

biom_sp_deg1 <- biom_sp_deg %>%
  dplyr::group_by(epibiont, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_sp_deg1)

barplot_biom_sp_deg1 <- ggplot(biom_sp_deg1, aes(x = inf.percent, y = mean_biom, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_degrad) +
  scale_fill_manual(values = colors_degrad) +
  scale_y_continuous(limits = c(0, 125), breaks = seq(from = 0, to = 125, by = 25)) +
  facet_wrap(~ epibiont, ncol = 3, labeller = label_wrap_gen(width=10)) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_biom_sp_deg1






# BARPLOT - Species biomass vs Position ####
biom_sp_pos <- biom_long[biom_long$epibiont %in% c("Ralfsia verrucosa", "Electra verticillata",
                                                   "Bryozoa", "Tetraclita serrata"),]

biom_sp_pos1 <- biom_sp_pos %>%
  dplyr::group_by(epibiont, position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_sp_pos1)

barplot_biom_sp_pos1 <- ggplot(biom_sp_pos1, aes(x = position, y = mean_biom, 
                                                 color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(from = 0, to = 125, by = 25)) +
  scale_x_discrete(labels = c("II epibio", " Inf", "Non-inf")) +
  facet_wrap(~ epibiont, ncol = 4, labeller = label_wrap_gen(width=10)) +
  labs(color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 8, face = "bold"))
barplot_biom_sp_pos1






# BARPLOT - TETSER - Infestation x Degradation ####
tetser_sp_biom <- biom_long[biom_long$epibiont %in% "Tetraclita serrata",]
View(tetser_sp_biom)

tetser_sp_biom1 <- tetser_sp_biom %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(tetser_sp_biom1)

tetser_sp_biom1[nrow(tetser_sp_biom1) + 1,] <- list("A", "1", NA, NA)

barplot_tetser_sp_biom1 <- ggplot(tetser_sp_biom1, aes(x = inf.lvl, y = mean_biom, 
                                                         color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 125), breaks = seq(from = 0, to = 150, by = 25)) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Tetraclita serrata", color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_tetser_sp_biom1






# BARPLOT - GELPRI - Infestation x Degradation ####
gelpri_sp_biom <- biom_long[biom_long$epibiont %in% "Gelidium pristoides",]
View(gelpri_sp_biom)

gelpri_sp_biom1 <- gelpri_sp_biom %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(gelpri_sp_biom1)

gelpri_sp_biom1[nrow(gelpri_sp_biom1) + 1,] <- list("A", "1", NA, NA)

barplot_gelpri_sp_biom1 <- ggplot(gelpri_sp_biom1, aes(x = inf.lvl, y = mean_biom, 
                                                       color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  scale_y_continuous(limits = c(0, 3000), breaks = seq(from = 0, to = 3000, by = 500)) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Gelidium pristoides", color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_gelpri_sp_biom1






# BARPLOT - HILLEC - Infestation x Position ####
hillec_sp_biom <- biom_long[biom_long$epibiont %in% "Hildenbrandia lecannellieri",]
View(hillec_sp_biom)

hillec_sp_biom1 <- hillec_sp_biom %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(hillec_sp_biom1)

hillec_sp_biom1[nrow(hillec_sp_biom1) + 1,] <- list("1", "infested", NA, NA)

barplot_hillec_sp_biom1 <- ggplot(hillec_sp_biom1, aes(x = position, y = mean_biom, 
                                                       color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_colour_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 175), breaks = seq(from = 0, to = 250, by = 50)) +
  scale_x_discrete(labels = c("II epibio", "Inf", "Non-inf")) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  labs(title = "Hildenbrandia lecannellieri", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_hillec_sp_biom1





# BARPLOT - GELPRI - Infestation x Position ####
View(gelpri_sp_biom)

gelpri_sp_biom2 <- gelpri_sp_biom %>%
  dplyr::group_by(inf.lvl, position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(gelpri_sp_biom2)

library(ggh4x)
strip1 <- strip_themed(background_x = elem_list_rect(fill = colors_inflvl))


barplot_gelpri_sp_biom2_text <- data.frame(
  label = c("B", "A", "A"),
  inf.lvl = c("D"),
  position = c("II epibiosis", "infested", "non infested"),
  y = c(550, 2700, 2550)
)

barplot_gelpri_sp_biom2 <- ggplot(gelpri_sp_biom2, aes(x = position, y = mean_biom, 
                                                       color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_gelpri_sp_biom2_text, aes(x = position, y = y, label = label), 
            color = "black", size = 3) +
  scale_colour_manual(values = colors_position2) +
  scale_fill_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 2700), breaks = seq(from = 0, to = 2500, by = 500)) +
  scale_x_discrete(labels = c("II epibio", "Inf", "Non-inf")) +
  facet_wrap2(~ inf.lvl, ncol = 5, strip = strip1) +
  labs(title = "Gelidium pristoides", color = "Position on the mussel shell", fill = "Position on the mussel shell",
       x = "Position on the mussel shell", y = "\n Mean biomass per mussel valve (mg)") +
  theme(axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none", 
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
barplot_gelpri_sp_biom2



