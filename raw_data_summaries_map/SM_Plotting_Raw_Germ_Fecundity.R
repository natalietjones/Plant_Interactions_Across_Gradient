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

############################################################
## 0. PATHS
############################################################

normalise_existing <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

find_project_root <- function(data_file_name = "fullWAdata.csv") {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  candidates <- unique(normalise_existing(c(
    wd,
    dirname(wd),
    dirname(dirname(wd)),
    dirname(dirname(dirname(wd)))
  )))
  
  candidates <- candidates[dir.exists(candidates)]
  
  has_data_file <- vapply(
    candidates,
    function(x) file.exists(file.path(x, "data", "field", data_file_name)),
    logical(1)
  )
  
  hits <- candidates[has_data_file]
  
  if (length(hits) == 0) {
    stop(
      "Could not find data/field/", data_file_name, ".\n",
      "Current working directory: ", wd, "\n",
      "Checked:\n",
      paste(file.path(candidates, "data", "field", data_file_name), collapse = "\n"),
      "\n\nSet the working directory to the project root or to raw_data_summaries_map/."
    )
  }
  
  hits[[1]]
}

project_dir <- find_project_root("fullWAdata.csv")

DATA_FILE <- file.path(project_dir, "data", "field", "fullWAdata.csv")
FIGDIR    <- file.path(project_dir, "figures")

dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DATA_FILE)) {
  stop("Cannot find DATA_FILE: ", DATA_FILE)
}

cat("\nPath setup:\n")
cat("  project_dir: ", normalizePath(project_dir, winslash = "/"), "\n", sep = "")
cat("  DATA_FILE:   ", normalizePath(DATA_FILE, winslash = "/"), "\n", sep = "")
cat("  FIGDIR:      ", normalizePath(FIGDIR, winslash = "/"), "\n\n", sep = "")

############################################################
## DATA
############################################################

df <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

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
site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

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
  "BEN" = "#89A6FB",
  "GH"  = "#3AB795",
  "NAM" = "#FFCF56",
  "PJ"  = "#FC7A57",
  "CS"  = "#FE4A49"
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


FIGDIR <- file.path("..", "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(FIGDIR, "Mean_Germination_Cover_TRT.png"),
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

FIGDIR <- file.path("..", "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(FIGDIR, "MeanSeedProduction_TRT.png"),
  p_seed,
  width = 7,
  height = 9,
  dpi = 300
)

############################################################
## 3. FREQUENCY DISTRIBUTION OF GRASS ABUNDANCE
############################################################
## Plots the distribution of NO_GRASS separately for each
## site and cover treatment.
##
## Figure layout:
## - 5 sites as columns
## - Sun and Shade as rows
## - compact spacing for A4 supplementary figure

grass_freq_dat <- df %>%
  filter(
    !is.na(NO_GRASS),
    (NO_GRASS>0),
    !is.na(Site),
    !is.na(Cov)
  ) %>%
  mutate(
    NO_GRASS = as.numeric(NO_GRASS),
    Site = factor(Site, levels = site_order),
    Cov = factor(
      Cov,
      levels = c("SUN", "SH"),
      labels = c("Sun", "Shade")
    )
  )

grass_rep_counts <- grass_freq_dat %>%
  group_by(Site, Cov) %>%
  summarise(
    n_total = n(),
    n_with_grass = sum(NO_GRASS >= 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_label = paste0("n ≥ 1 grass = ", n_with_grass)
  )


p_grass_freq <- ggplot(
  grass_freq_dat,
  aes(x = NO_GRASS)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    fill = "grey75",
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = grass_rep_counts,
    aes(
      x = Inf,
      y = Inf,
      label = n_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.3,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  facet_grid(
    Cov ~ Site,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = function(x) pretty(x, n = 4),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 3),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = "Number of grasses",
    y = "Frequency"
  ) +
  theme_classic(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    
    axis.title = element_text(size = 15, color = "black"),
    axis.text  = element_text(size = 11, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.line  = element_line(color = "black", linewidth = 0.4),
    
    strip.text = element_text(size = 13, color = "black"),
    strip.background = element_blank(),
    
    panel.spacing = unit(0.15, "lines"),
    plot.margin = margin(4, 4, 4, 4, unit = "pt"),
    
    legend.position = "none"
  )

p_grass_freq


ggsave(
  file.path(FIGDIR, "NO_GRASS_frequency_by_site_cover.png"),
  p_grass_freq,
   width = 11.69,
  height = 5.2,
  units = "in",
  dpi = 600,
  bg = "white"
)