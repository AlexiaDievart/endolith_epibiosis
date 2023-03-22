## ---------------------------
##
## Script name: Life on collapsing shells - the relationship between epibiosis 
##              and euendolithic infestation in the brown mussel Perna perna
##
##
## Purpose of script: 
##        - Shell parameters
##        - Community descriptors 1
##        - Community descriptors 2
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

library(ggpubr)
remove.packages("ggplot2")
install.packages("ggplot2")
library(ggplot2)

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

# Load full dataset (long format) =============================================================
epibio <- read.csv("epibiosis.csv", h=T, dec=",", sep=";")
dplyr::glimpse(epibio)
factor_epibio <- c("location", "quadrat", "specimen", "mussel.species", "sex", "valve", "inf.lvl", "inf.percent",
                   "higher.group", "epibiont", "position", "basibiont")
epibio[,factor_epibio] <- lapply(epibio[,factor_epibio], factor)
epibio$inf.percent <- recode_factor(epibio$inf.percent, "0" = "1", "0-25" = "2", "25-50" = "3", "50-75" = "4",
                                    "75-100" = "5", "100" = "6")
epibio <- epibio[! epibio$mussel.species == "Mytilus galloprovincialis",]
View(epibio)

# Community variables for abundance ===========================================================
abund_long <- epibio[,1:16]
View(abund_long)

## Per epibiotic species
abund_long_species <- abund_long[, -12]
View(abund_long_species)

## Calculate the mean shell length for each euendolithic infestation level
length_inf <- epibio %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl) %>%
  dplyr::summarise(
    shell_length = mean(shell.length)
  )
View(length_inf)
# For each mussel, each valve can have different euendolithic infestation levels.

length_inf2 <- length_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length),
    n = n(),
    se_length  = sd_length / sqrt(n)
  )
head(length_inf2)

## Calculate the total number of valves and the number of valves w/ epibionts for each infestation level and shell degradation index
prev_epibio_wide <- abund_long %>%
  dplyr::group_by(quadrat, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    n = n_distinct(specimen),
    mean_length = mean(shell.length),
    n.clean = sum(abundance == "0"),
    n.inf = n - n.clean
  )
View(prev_epibio_wide)

prev_epibio_long <- prev_epibio_wide %>%
  pivot_longer(names_to = "epibiosis", values_to = "N", cols = -c(quadrat, valve, inf.lvl, inf.percent, n, mean_length))
View(prev_epibio_long)

prev_epibio_wide1 <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell.length),
    epibiosis = sum(abundance)
  )
View(prev_epibio_wide1)

prev_epibio_wide1$epibiosis[prev_epibio_wide1$epibiosis != 0] <- 1
prev_epibio_wide1$epibiosis <- as.factor(prev_epibio_wide1$epibiosis)


###############################################################################################
# Section: Graphical analyses -----------------------------------------------------------------
###############################################################################################

# 1. EPIBIOSIS and SHELL PARAMETERS ==============================================================

figure1 <- ggarrange(hist_infxdeg1, barplot_length_inf3, barplot_length_inf4,
                     prev_dens, prev_hist1, prev_hist2,
                          ncol = 3, nrow = 2, common.legend = F,
                          labels = c("A", "B", "C", "D", "E", "F"), 
                          hjust = c(-0.5, -1.5, -1.5, -0.5, -1.5, -1.5),
                          font.label = list(size = 20, color = "grey30"))
figure1

# HIST - Degradation x Infestation ####

# Summarize infestation and degradation for each valve of each specimen
infxdeg1 <- infxdeg %>%
  dplyr::group_by(inf.lvl, inf.percent) %>%
  dplyr::summarise(
    N = sum(n)
  )
View(infxdeg1)

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


# BARPLOT - Mean shell length x Infestation ####

# Summarize shell length for each infestation level
length_infxdeg <- epibio_perna %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    shell_length = mean(shell.length)
  )
View(length_infxdeg)

length_infxdeg3 <- length_infxdeg %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length)
  )
View(length_infxdeg3)

barplot_length_inf3 <- ggplot(length_infxdeg3, aes(x = inf.lvl, y = mean_length, 
                                                   color = inf.lvl, fill = inf.lvl)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = c(70, 80, 95, 105, 100),
           label = c("A", "A", "A", "B", "B"),
           size = 4) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100)) +
  scale_color_manual(values = colors_inflvl) +
  scale_fill_manual(values = colors_inflvl) +
  labs(color = "Infestation level", fill = "Infestation level",
       x = "Infestation levels", y = "\n Mean shell length (in mm)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_length_inf3

# BARPLOT - Mean shell length x Degradation ####

# Summarize shell length by shell degradation index
length_infxdeg4 <- length_infxdeg %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell_length),
    sd_length  = sd(shell_length)
  )
View(length_infxdeg4)

barplot_length_inf4 <- ggplot(length_infxdeg4, aes(x = inf.percent, y = mean_length, 
                                                   color = inf.percent, fill = inf.percent)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_length - sd_length, ymax = mean_length + sd_length), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  annotate("text", 
           x = c(1, 2, 3, 4, 5, 6),
           y = c(55, 70, 85, 95, 105, 100),
           label = c("A", "B", "C", "D", "E", "E"),
           size = 4) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100)) +
  scale_color_manual(values = col_deg) +
  scale_fill_manual(values = col_deg) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "SD Index", fill = "SD Index",
       x = "Shell degradation index", y = "\n Mean shell length (in mm)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold")) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_length_inf4

# DENSITY PLOT - Epibiosis x Shell length ####
prev_epibio_wide1 <- abund_long %>%
  dplyr::group_by(quadrat, specimen, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    mean_length = mean(shell.length),
    epibiosis = sum(abundance)
  )
View(prev_epibio_wide1)

prev_epibio_wide1$epibiosis[prev_epibio_wide1$epibiosis != 0] <- 1
prev_epibio_wide1$epibiosis <- as.factor(prev_epibio_wide1$epibiosis)

prev_dens <- ggplot(data = prev_epibio_wide1, aes(x = mean_length, y = after_stat(count), group = epibiosis)) +
  geom_density(aes(fill = epibiosis), position = position_fill(), size = 0) +
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

# HIST - Epibiosis x Infestation ####

# Summarize the number of mussels with and without epibionts per infestation levels
prev_epibio_wide <- epibio_perna %>%
  dplyr::group_by(quadrat, valve, inf.lvl, inf.percent) %>%
  dplyr::summarise(
    n = n_distinct(specimen),
    mean_length = mean(shell.length),
    n.clean = sum(abundance == "0"),
    n.inf = n - n.clean
  )
View(prev_epibio_wide)

prev_epibio_long <- prev_epibio_wide %>%
  pivot_longer(names_to = "epibiosis", values_to = "N", cols = -c(quadrat, valve, inf.lvl, inf.percent, n, mean_length))
View(prev_epibio_long)

prev_hist1 <- ggplot(prev_epibio_long, aes(x = inf.lvl, y = N, fill = epibiosis)) +
  geom_col(position = "fill") +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 1.05,
           label = c("AB", "AB", "A", "B", "AB"),
           size = 4) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = colors_epibiosis) +
  labs(y = "\n Proportion of mussel valves", x = "Infestation levels",
       fill = "Epibiosis") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), 
        legend.position = "none")
prev_hist1

# HIST - Epibiosis x Degradation ####
prev_hist2 <- ggplot(prev_epibio_long, aes(x = inf.percent, y = N, fill = epibiosis)) +
  geom_col(position = "fill") +
  annotate("text", 
           x = c(1, 2, 3, 4, 5, 6),
           y = 1.05,
           label = c("A", "AB", "B", "AB", "B", "AB"),
           size = 4) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = colors_epibiosis) +
  labs(y = "\n Proportion of mussel valves", x = "Shell degradation index",
       fill = "Epibiosis") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), 
        legend.position = "none")
prev_hist2





# 2. EPIBIOTIC COMMUNITY DESCRIPTORS 1 - Abundance, percent cover and biomass ====

figure2 <- ggarrange(ggarrange(barplot_abund_length2, barplot_abund_inf2, barplot_abund_deg2, barplot_abund_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("A", "B", "C", "D"), 
                               hjust = c(-1.3, -2.5, -2.5, -2.5)),
                     ggarrange(barplot_cover_length3, barplot_cover_inf3, barplot_cover_deg4,
                               ncol = 3, common.legend = F, widths = c(2, 1, 2),
                               labels = c("E", "F", "G"),
                               hjust = c(-1.5, -3, -2.5)), 
                     ggarrange(barplot_biom_length2, barplot_biom_inf2, barplot_biom_deg2, barplot_biom_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("H", "I", "J", "K"), 
                               hjust = c(-1.5, -9, -3.5, -2.5),
                               vjust = c(0.5, 1, 1, 1)),
                     nrow = 3)
figure2


# BARPLOT - Mean total abundance x Shell length ####
abund_epibio_length2 <- abund_epibio_length %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_length2)

barplot_abund_length2 <- ggplot(abund_epibio_length2, aes(x = length_range, y = mean_abund, 
                                                          color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_abund - 0, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(12, 15, 17, 25, 35, 40, 45, 60, 40),
           label = c("A", "ABC", "ABCD", "BDE", "BEF", "BFG", "BG", "BG", "BG"),
           size = 4) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(breaks = seq(from = 0, to = 60, by = 10)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total abundance per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none")
barplot_abund_length2

# BARPLOT - Mean total abundance x Infestation ####
abund_epibio_inf2 <- abund_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_inf2)

barplot_abund_inf2 <- ggplot(abund_epibio_inf2, aes(x = inf.lvl, y = mean_abund, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(13, 14, 22, 20, 19),
           label = c("AB", "B", "A", "A", "AB"),
           size = 4) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  labs(color = "Infestation level", fill = "Infestation level",
       x = "Infestation levels", y = "\n Mean total abundance per valve") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = c(0.02, 0.99),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_abund_inf2

# BARPLOT - Mean total abundance x Degradation ####
abund_epibio_deg2 <- abund_epibio_inf %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_deg2)

barplot_abund_deg2 <- ggplot(abund_epibio_deg2, aes(x = inf.percent, y = mean_abund, 
                                                    color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = mean_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6),
           y = c(8, 14, 15, 19, 22, 18),
           label = c("A", "AB", "AB", "B", "B", "AB"),
           size = 4) +
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "SD Index", fill = "SD Index",
       x = "Shell degradation index", y = "\n Mean total abundance per valve") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.02, 0.99),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2))
barplot_abund_deg2

# BARPLOT - Mean total abundance x Position ####
abund_epibio_pos2 <- abund_epibio_pos %>%
  dplyr::group_by(position) %>%
  dplyr::summarise(
    mean_abund = mean(abundance),
    sd_abund = sd(abundance),
  )
View(abund_epibio_pos2)

barplot_abund_pos2 <- ggplot(abund_epibio_pos2, aes(x = position, y = mean_abund, 
                                                    color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  annotate("text", 
           x = c(1, 2, 3),
           y = c(19, 13, 11),
           label = c("A", "B", "C"), size = 4) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(from = 0, to = 20, by = 5)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean total abundance per valve") +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.22, 0.99),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2))
barplot_abund_pos2


# BARPLOT - Cover x Shell length ####
cover_epibio_length1 <- cover_epibio_length %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  dplyr::summarise(
    cover = mean(cover)
  )
View(cover_epibio_length1)

cover_epibio_length3 <- cover_epibio_length1 %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_length3)

barplot_cover_length3 <- ggplot(cover_epibio_length3, aes(x = length_range, y = mean_cover, 
                                                          color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(13, 15, 27, 34, 27, 32, 25, 27, 7),
           label = c("A", "AB", "BC", "C", "C", "C", "C", "BC", "AC"),
           size = 4) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(from = 0, to = 35, by = 5), 
                     labels = scales::percent_format(scale = 1)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean percent cover per mussel") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 11),
        legend.position = "none")
barplot_cover_length3

# BARPLOT - Cover x Infestation ####
cover_epibio_inf3 <- cover_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_inf3)

barplot_cover_inf3 <- ggplot(cover_epibio_inf3, aes(x = inf.lvl, y = mean_cover, 
                                                    color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(from = 0, to = 35, by = 5), 
                     labels = scales::percent_format(scale = 1)) +
  labs(color = "Infestation level", fill = "Infestation level",
       x = "Infestation levels", y = "\n Mean percent cover per valve") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_cover_inf3

# BARPLOT - Cover x (Degradation x Position) ####
cover_epibio_deg3 <- cover_long1 %>%
  dplyr::group_by(quadrat, specimen, valve, inf.percent, position) %>%
  dplyr::summarise(
    cover = sum(cover.percent)
  )
View(cover_epibio_deg3)

cover_epibio_deg4 <- cover_epibio_deg3 %>%
  dplyr::group_by(inf.percent, position) %>%
  dplyr::summarise(
    mean_cover = mean(cover),
    sd_cover = sd(cover)
  )
View(cover_epibio_deg4)

library(ggh4x)
strip <- strip_themed(background_x = elem_list_rect(fill = col_deg))

barplot_cover4_text <- data.frame(
  label = c("", "BC", "BC", "A", "B", "A", "B", "AC", "B", "B", "BC", "BC"),
  inf.percent = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6"),
  position     = c("infested", "non infested", "infested", "non infested", "infested", "non infested", "infested", "non infested", "infested", "non infested", "infested", "non infested"),
  y     = c(0, 9, 9, 42, 12, 41, 16, 31, 17, 22, 21, 12)
)

barplot_cover_deg4 <- ggplot(cover_epibio_deg4, aes(x = position, y = mean_cover, 
                                            color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_cover, ymax = mean_cover + sd_cover), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(data = barplot_cover4_text, aes(x = position, y = y, label = label), color = "black", size = 4) +
  facet_wrap2(~ inf.percent, ncol = 6, strip = strip) +
  scale_fill_manual(values = colors_position1) +
  scale_color_manual(values = colors_position1) +
  scale_y_continuous(limits = c(0, 42), breaks = seq(from = 0, to = 40, by = 5), 
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels = c("Inf", "Non-inf")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean percent cover per valve") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"))
barplot_cover_deg4

# BARPLOT - Biomass x Shell length ####
biom_epibio_length <- biom_long1 %>%
  dplyr::group_by(quadrat, specimen, length_range) %>%
  dplyr::summarise(
    biomass = sum(biomass)
  )
View(biom_epibio_length)

biom_epibio_length2 <- biom_epibio_length %>%
  dplyr::group_by(length_range) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_length2)

barplot_biom_length2 <- ggplot(biom_epibio_length2, aes(x = length_range, y = mean_biom, 
                                                        color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(60, 80, 130, 400, 650, 1020, 500, 620, 100),
           label = c("A", "AB", "B", "C", "CD", "D", "CD", "CD", "ABCD"),
           size = 4) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1000, by = 250)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean total biomass per mussel (mg)") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 11, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1),
        legend.position = "none")
barplot_biom_length2

# BARPLOT - Biomass x Infestation ####
biom_epibio_inf2 <- biom_epibio_inf %>%
  dplyr::group_by(inf.lvl) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_inf2)

barplot_biom_inf2 <- ggplot(biom_epibio_inf2, aes(x = inf.lvl, y = mean_biom, 
                                                  color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(420, 180, 340, 460, 120),
           label = c("ABC", "A", "B", "C", "ABC"),
           size = 4) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(limits = c(0, 460), breaks = seq(from = 0, to = 450, by = 100)) +
  labs(color = "Infestation level", fill = "Infestation level",
       x = "Infestation levels", y = "\n Mean total biomass per valve (mg)") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_biom_inf2

# BARPLOT - Biomass x Degradation ####
biom_epibio_deg2 <- biom_epibio_inf %>%
  dplyr::group_by(inf.percent) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_deg2)

barplot_biom_deg2 <- ggplot(biom_epibio_deg2, aes(x = inf.percent, y = mean_biom, 
                                                  color = inf.percent, fill = inf.percent)) +
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6),
           y = c(40, 510, 470, 500, 290, 140),
           label = c("A", "BC", "C", "CD", "AB", "ABD"),
           size = 4) +
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  scale_y_continuous(limits = c(0, 510), breaks = seq(from = 0, to = 500, by = 100)) +
  labs(color = "SD Index", fill = "SD Index",
       x = "Shell degradation index", y = "\n Mean total biomass per valve (mg)") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_biom_deg2

# BARPLOT - Biomass x Position ####
biom_epibio_pos2 <- biom_epibio_pos %>%
  dplyr::group_by(position) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass)
  )
View(biom_epibio_pos2)

barplot_biom_pos2 <- ggplot(biom_epibio_pos2, aes(x = position, y = mean_biom, 
                                                  color = position, fill = position)) + 
  geom_errorbar(aes(ymin = mean_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 400), breaks = seq(from = 0, to = 400, by = 100)) +
  annotate("text", 
           x = c(1, 2, 3),
           y = c(200, 190, 390),
           label = c("A", "B", "A"), size = 4) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean total biomass per valve (mg)") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_biom_pos2





# 3. EPIBIOTIC COMMUNITY DESCRIPTORS 2 - Species richness, diversity and evenness ====

figure3 <- ggarrange(ggarrange(barplot_rich_length2, barplot_rich_inf2, barplot_rich_deg2, barplot_rich_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("A", "B", "C", "D"), 
                               hjust = c(-1.2, -2.4, -2.4, -2.4)),
                     ggarrange(barplot_shannon_length2, barplot_shannon_inf2, barplot_shannon_deg2, barplot_shannon_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("E", "F", "G", "H"),
                               hjust = c(-1.5, -3, -2.2, -2.2),
                               vjust = c(0.8, 0.8, 0.8, 0.8)), 
                     ggarrange(barplot_simpson_length2, barplot_simpson_inf2, barplot_simpson_deg2, barplot_simpson_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("I", "J", "K", "L"), 
                               hjust = c(-4.5, -3, -2.2, -2.2), 
                               vjust = c(0.8, 0.8, 0.8, 0.8)),
                     ggarrange(barplot_pielou_length2, barplot_pielou_inf2, barplot_pielou_deg2, barplot_pielou_pos2,
                               ncol = 4, common.legend = F, widths = c(2, 1, 1, 1),
                               labels = c("M", "N", "O", "P"),
                               hjust = c(-1.8, -2.8, -2.8, -3.1),
                               vjust = c(0.8, 0.8, 0.8, 0.8)), 
                     nrow = 4)
figure3

# BARPLOT - Mean species richness x Shell length ####
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
           y = c(3, 3, 4, 5.7, 6.2, 6.5, 6.7, 8, 6.5),
           label = c("A", "A", "A", "B", "BC", "CD", "CD", "D", "BCD"),
           size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(breaks = seq(from = 0, to = 8, by = 1)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean species richness (S)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_rich_length2

# BARPLOT - Mean species richness x Infestation ####
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
           size = 3) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(breaks = seq(from = 0, to = 3, by = 1)) +
  labs(color = "Infestation level", fill = "Infestation level",
       x = "Infestation levels", y = "\n Mean species richness (S)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_rich_inf2

# BARPLOT - Mean species richness x Degradation ####
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
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  scale_y_continuous(limits = c(0, 4.2), breaks = seq(from = 0, to = 4, by = 1)) +
  labs(color = "SD Index", fill = "SD Index",
       x = "Shell degradation index", y = "\n Mean species richness (S)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_rich_deg2

# BARPLOT - Mean species richness x Position ####
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
           y = c(2.8, 3.2, 2.7),
           label = c("A", "B", "A"),
           size = 3) +
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 3.2), breaks = seq(from = 0, to = 3, by = 1)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean species richness (S)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_rich_pos2

# BARPLOT - Mean Shannon x Shell length ####
shannon_length2 <- abund_div_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    shannon_mean = mean(shannon_length),
    shannon_sd = sd(shannon_length)
  )
View(shannon_length2)

barplot_shannon_length2 <- ggplot(shannon_length2, aes(x = length_range, y = shannon_mean, 
                                                       color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(0.7, 0.8, 1.1, 1.4, 1.6, 1.6, 1.6, 1.8, 1.4),
           label = c("A", "A", "A", "B", "BC", "C", "BC", "C", "ABC"),
           size = 3) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(from = 0, to = 2, by = 0.5)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Shannon's diversity (H')") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_shannon_length2

# BARPLOT - Mean Shannon x Infestation ####
shannon_inf2 <- abund_div_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_inf2)

barplot_shannon_inf2 <- ggplot(shannon_inf2, aes(x = inf.lvl, y = shannon_mean, 
                                                 color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(0.7, 0.9, 1.2, 1.3, 1.3),
           label = c("AB", "A", "B", "C", "ABC"),
           size = 3) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1.5, by = 0.5)) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean Shannon's diversity (H')") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = c(0.02, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(4, "mm"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_shannon_inf2

# BARPLOT - Mean Shannon x Degradation ####
shannon_deg2 <- abund_div_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    shannon_mean = mean(shannon_inf),
    shannon_sd = sd(shannon_inf)
  )
View(shannon_deg2)

barplot_shannon_deg2 <- ggplot(shannon_deg2, aes(x = inf.percent, y = shannon_mean, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(from = 0, to = 1.5, by = 0.5)) +
  labs(color = "SD Index", fill = "SD Index",
       x = "Shell degradation index", y = "\n Mean Shannon's diversity (H')") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = c(0.02, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(4, "mm"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_shannon_deg2

# BARPLOT - Mean Shannon x Position ####
shannon_pos2 <- abund_div_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    shannon_mean = mean(shannon_pos),
    shannon_sd = sd(shannon_pos)
  )
View(shannon_pos2)

barplot_shannon_pos2 <- ggplot(shannon_pos2, aes(x = position, y = shannon_mean, 
                                                 color = position, fill = position)) + 
  geom_errorbar(aes(ymin = shannon_mean, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3),
           y = c(0.75, 0.95, 0.75),
           label = c("AB", "A", "B"),
           size = 3) +
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean Shannon's diversity (H')") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_shannon_pos2

# BARPLOT - Mean Simpson x Shell length ####
simpson_length2 <- abund_div_length %>%
  dplyr::group_by(length_range) %>%
  summarise(
    simpson_mean = mean(simpson_length),
    simpson_sd = sd(simpson_length)
  )
View(simpson_length2)

barplot_simpson_length2 <- ggplot(simpson_length2, aes(x = length_range, y = simpson_mean, 
                                                       color = length_range, fill = length_range)) + 
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
           y = c(0.4, 0.5, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7),
           label = c("A", "A", "A", "B", "BC", "C", "BC", "BC", "ABC"),
           size = 3) +
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Simpson's diversity (λ)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_simpson_length2

# BARPLOT - Mean Simpson x Infestation ####
simpson_inf2 <- abund_div_inf %>%
  dplyr::group_by(inf.lvl) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_inf2)

barplot_simpson_inf2 <- ggplot(simpson_inf2, aes(x = inf.lvl, y = simpson_mean, 
                                                 color = inf.lvl, fill = inf.lvl)) + 
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  annotate("text",
           x = c(1, 2, 3, 4, 5),
           y = c(0.5, 0.6, 0.7, 0.8, 0.75),
           label = c("AB", "A", "B", "C", "ABC"),
           size = 3) +
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean Simpson's diversity (λ)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_simpson_inf2

# BARPLOT - Mean Simpson x Degradation ####
simpson_deg2 <- abund_div_inf %>%
  dplyr::group_by(inf.percent) %>%
  summarise(
    simpson_mean = mean(simpson_inf),
    simpson_sd = sd(simpson_inf)
  )
View(simpson_deg2)

barplot_simpson_deg2 <- ggplot(simpson_deg2, aes(x = inf.percent, y = simpson_mean, 
                                                 color = inf.percent, fill = inf.percent)) + 
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean Simpson's diversity (λ)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = "none")
barplot_simpson_deg2

# BARPLOT - Mean Simpson x Position ####
simpson_pos2 <- abund_div_pos %>%
  dplyr::group_by(position) %>%
  summarise(
    simpson_mean = mean(simpson_pos),
    simpson_sd = sd(simpson_pos)
  )
View(simpson_pos2)

barplot_simpson_pos2 <- ggplot(simpson_pos2, aes(x = position, y = simpson_mean, 
                                                 color = position, fill = position)) +
  geom_errorbar(aes(ymin = simpson_mean, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9), color = "black") +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean Simpson's diversity (λ)") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 9),
        legend.position = c(0.02, 1),
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(4, "mm"),
        legend.background = element_rect(fill='transparent')) +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2))
barplot_simpson_pos2

# BARPLOT - Mean Pielou x Shell length ####
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
  scale_colour_hue(h = c(180, 270)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(from = 0, to = 1, by = 0.25)) +
  scale_x_discrete(labels = c("40" = "]30-40]", "50" = "]40-50]",
                              "60" = "]50-60]","70" = "]60-70]",
                              "80" = "]70-80]", "90" = "]80-90]",
                              "100" = "]90-100]", "110" = "]100-110]",
                              "120" = "]110-120]")) +
  labs(color = "Shell length range", fill = "Shell length range",
       x = "Shell length range", y = "Mean Pielou's evenness (J')") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 9, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 7, hjust = 0.5),
        legend.position = "none")
barplot_pielou_length2

# BARPLOT - Mean Pielou x Infestation ####
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
  scale_color_manual(values = col_inf) +
  scale_fill_manual(values = col_inf) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Infestation levels", fill = "Infestation levels",
       x = "Infestation levels", y = "\n Mean Pielou's evenness (J')") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_pielou_inf2

# BARPLOT - Mean Pielou x Degradation ####
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
  scale_fill_manual(values = col_deg) +
  scale_color_manual(values = col_deg) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Shell degradation index", fill = "Shell degradation index",
       x = "Shell degradation index", y = "\n Mean Pielou's evenness (J')") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_pielou_deg2

# BARPLOT - Mean Pielou x Position ####
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
  scale_fill_manual(values = colors_position2) +
  scale_color_manual(values = colors_position2) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(from = 0, to = 1, by = 0.25)) +
  labs(color = "Position on the shell", fill = "Position on the shell",
       x = "Position on the shell", y = "\n Mean Pielou's evenness (J')") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
barplot_pielou_pos2