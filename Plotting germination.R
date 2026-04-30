############################################################
## PACKAGES
############################################################

packages <- c(
  "ggplot2", "dplyr", "tidyverse", "ggpattern", "colorspace",
  "glmmTMB", "DHARMa", "emmeans"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

dir.create("figures", showWarnings = FALSE)


############################################################
## DATA
############################################################

df <- read.csv("fullWAdata.csv")

# Convert germination Y/N to 1/0
df <- df %>%
  mutate(
    FOCAL_GERM = case_when(
      FOCAL_GERM == "Y" ~ 1,
      FOCAL_GERM == "N" ~ 0,
      TRUE ~ NA_real_
    )
  )

# Factor ordering
site_order <- c("BEN", "PJ", "NAM", "GH", "CS")

dens_trt_order <- c(
  "lambda",
  "High density monoculture",
  "Invasion when rare (IWR)",
  "IWR + Facilitator"
)

dens_trt_labels <- c(
  "lambda" = "λ",
  "High density monoculture" = "High\ndensity",
  "Invasion when rare (IWR)" = "IWR",
  "IWR + Facilitator" = "IWR +\nfacilitator"
)

focsp_labels <- c(
  "TRCY" = "italic('T. cyanopetala')",
  "TROR" = "italic('T. ornata')"
)

# Site colours
site_colors <- c(
  "BEN" = "#2B83BA",
  "PJ"  = "#ABDDA4",
  "NAM" = "#FFFFBF",
  "GH"  = "#FDAE61",
  "CS"  = "#D7191C"
)

# Shared dodge position so bars and error bars align
pd <- position_dodge2(width = 0.8, preserve = "single")


############################################################
## SHARED PLOT THEME
############################################################

plot_theme <- theme_bw(base_size = 14) +
  theme(
    text              = element_text(color = "black"),
    axis.title        = element_text(size = 13, color = "black"),
    axis.text         = element_text(size = 13, color = "black"),
    axis.text.x       = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y       = element_text(color = "black"),
    axis.ticks        = element_line(color = "black"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text.x = element_text(size = 13, color = "black"),
    strip.text.y      = element_text(size = 13, color = "black"),
    strip.text.y.left = element_text(size = 13, color = "black"),
    strip.background  = element_blank(),
    strip.placement   = "outside",
    legend.title      = element_text(size = 13, color = "black"),
    legend.text       = element_text(size = 13, color = "black")
  )


############################################################
## 1. MEAN GERMINATION BY COVER
############################################################
## Calculation check:
## - germ_mean is the proportion germinated.
## - germ_pct converts this to percent germination.
## - n_germ counts only non-missing germination observations.
## - SE uses binomial SE: sqrt(p * (1 - p) / n).

germ_cov_summary <- df %>%
  group_by(dens_trt, focsp, Site, Cov) %>%
  summarise(
    n_total   = n(),
    n_germ    = sum(!is.na(FOCAL_GERM)),
    germ_prop = mean(FOCAL_GERM, na.rm = TRUE),
    germ_pct  = 100 * germ_prop,
    germ_se   = sqrt(germ_prop * (1 - germ_prop) / n_germ),
    germ_se_pct = 100 * germ_se,
    .groups = "drop"
  ) %>%
  mutate(
    germ_prop   = ifelse(n_germ == 0, NA_real_, germ_prop),
    germ_pct    = ifelse(n_germ == 0, NA_real_, germ_pct),
    germ_se_pct = ifelse(n_germ == 0, NA_real_, germ_se_pct),
    ymin = pmax(germ_pct - germ_se_pct, 0),
    ymax = pmin(germ_pct + germ_se_pct, 100),
    Site     = factor(Site, levels = site_order),
    dens_trt = factor(dens_trt, levels = dens_trt_order),
    focsp    = factor(focsp),
    Cov      = factor(Cov, levels = c("SUN", "SH")),
    alpha_val = ifelse(Cov == "SUN", 0.5, 1)
  )

p_germ_cover <- ggplot(
  germ_cov_summary,
  aes(
    x     = Site,
    y     = germ_pct,
    fill  = Site,
    alpha = alpha_val,
    group = Cov
  )
) +
  geom_col(
    position  = pd,
    color     = "black",
    linewidth = 0.25,
    width     = 0.75
  ) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width     = 0.0,
    position = position_dodge(0.8), 
    linewidth = 0.35,
    color     = "black"
  ) +
  facet_grid(
    dens_trt ~ focsp,
    labeller = labeller(
      dens_trt = as_labeller(dens_trt_labels),
      focsp = as_labeller(focsp_labels, default = label_parsed)
    )
  )+
  scale_fill_manual(values = site_colors) +
  scale_alpha_identity(
    name   = "Cover",
    breaks = c(0.5, 1),
    labels = c("Sun", "Shade"),
    guide  = "legend"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    y = "Germination (%)",
    x = "Site"
  ) +
  plot_theme

p_germ_cover

ggsave(
  "figures/Mean_Germination_Cover_TRT.png",
  p_germ_cover,
  width = 7,
  height = 9,
  dpi = 300
)


############################################################
## 2. MEAN SEED PRODUCTION BY COVER
############################################################
## Calculation check:
## This calculates mean seed production only among germinated individuals.
## This is appropriate if the plot is intended to show fecundity conditional
## on germination. Non-germinated individuals are not included here.

df_germinated <- df %>%
  filter(FOCAL_GERM == 1)

fitness_summary <- df_germinated %>%
  group_by(dens_trt, Site, Cov, focsp) %>%
  summarise(
    n_total      = n(),
    n_fitness    = sum(!is.na(fitness)),
    mean_fitness = mean(fitness, na.rm = TRUE),
    se_fitness   = sd(fitness, na.rm = TRUE) / sqrt(n_fitness),
    .groups      = "drop"
  ) %>%
  mutate(
    se_fitness = ifelse(n_fitness <= 1, NA_real_, se_fitness),
    ymin = pmax(mean_fitness - se_fitness, 0),
    ymax = mean_fitness + se_fitness,
    Site      = factor(Site, levels = site_order),
    dens_trt  = factor(dens_trt, levels = dens_trt_order),
    focsp     = factor(focsp),
    Cov       = factor(Cov, levels = c("SUN", "SH")),
    alpha_val = ifelse(Cov == "SUN", 0.5, 1)
  )

p_seed <- ggplot(
  fitness_summary,
  aes(
    x     = Site,
    y     = mean_fitness,
    fill  = Site,
    alpha = alpha_val,
    group = Cov
  )
) +
  geom_col(
    position  = pd,
    color     = "black",
    linewidth = 0.25,
    width     = 0.75
  ) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width     = 0.0,
    position = position_dodge(0.8), 
    linewidth = 0.35,
    color     = "black"
  ) +
  facet_grid(
    dens_trt ~ focsp,
    labeller = labeller(
      dens_trt = as_labeller(dens_trt_labels),
      focsp = as_labeller(focsp_labels, default = label_parsed)
    )
  )  +
  scale_fill_manual(values = site_colors) +
  scale_alpha_identity(
    name   = "Cover",
    breaks = c(0.5, 1),
    labels = c("Sun", "Shade"),
    guide  = "legend"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    y = "Mean seed production per germinated plant (± SE)",
    x = "Site"
  ) +
  plot_theme

p_seed

ggsave(
  "figures/MeanSeedProduction_TRT.png",
  p_seed,
  width = 7,
  height = 9,
  dpi = 300
)


############################################################
## GLMMTMB — unchanged
############################################################

df$Site     <- factor(df$Site)
df$Cov      <- factor(df$Cov)
df$dens_trt <- factor(df$dens_trt, levels = dens_trt_order)
df$BLK      <- factor(df$BLK)

mod_nb <- glmmTMB(
  fitness ~ Site * Cov * dens_trt + (1 | BLK),
  family = nbinom2,
  data   = df
)

emm <- emmeans(mod_nb, ~ Site * Cov * dens_trt)
pairs(emm, by = c("Site", "Cov"))