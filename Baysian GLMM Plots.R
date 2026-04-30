############################################################
## PLOT PREDICTED FECUNDITY UNDER LAMBDA / HIGH COMPETITION
## final clean version
############################################################

rm(list = ls())

############################
## 0. Packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "patchwork",
  "colorspace"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(brms)
library(patchwork)
library(colorspace)

############################
## 1. Paths and files
############################
analysis_id <- "env_fac_model_comparison_v1"
run_dir     <- file.path(getwd(), analysis_id)
model_dir   <- file.path(run_dir, "fac_context_test_with_CS")

best_name <- scan(
  file.path(model_dir, "best_fac_context_model_with_CS.txt"),
  what = "character",
  quiet = TRUE
)

fit_file_map <- c(
  nofac     = "fit_site_cov_nofac_with_CS.rds",
  linear    = "fit_site_cov_fac_linear_with_CS.rds",
  quad      = "fit_site_cov_fac_quad_with_CS.rds",
  site_quad = "fit_site_quadfac_with_CS.rds",
  cov_quad  = "fit_cov_quadfac_with_CS.rds",
  both_quad = "fit_both_quadfac_with_CS.rds"
)

best_fit <- readRDS(file.path(model_dir, fit_file_map[[best_name]]))

dat <- read.csv(file.path(run_dir, "fecundity_analysis_data_all_sites.csv"))

############################
## 2. Settings
############################
site_levels <- c("Ben", "GH", "Nam", "PJ", "CS")

site_colors <- c(
  "Ben" = "#6F2DBD",
  "GH"  = "#6184D8",
  "Nam" = "#48BF84",
  "PJ"  = "#FFA400",
  "CS"  = "#D7191C"
)

lambda_colors <- darken(site_colors, amount = 0.20)

dat <- dat %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    Cov   = factor(Cov, levels = c("Shade", "Sun")),
    focsp = factor(focsp, levels = c("TRCY", "TROR")),
    BLK   = factor(BLK)
  )

fac_center <- mean(dat$fac_sc, na.rm = TRUE)
fac_scale  <- sd(dat$fac_sc, na.rm = TRUE)

############################
## 3. Site x focal-species 90th percentiles
############################
q90_tbl <- dat %>%
  group_by(focsp, Site) %>%
  summarise(
    fac_90   = quantile(fac_sc,   probs = 0.9, na.rm = TRUE),
    intra_90 = quantile(intra_sc, probs = 0.9, na.rm = TRUE),
    inter_90 = quantile(inter_sc, probs = 0.9, na.rm = TRUE),
    .groups = "drop"
  )

############################
## 4. Prediction grid
##
## Conditions within each cover:
## 1) Lambda (no neighbours, no facilitator)
## 2) Intra_noFac
## 3) Intra_withFac
## 4) Inter_noFac
## 5) Inter_withFac
############################
cond_levels <- c(
  "Lambda",
  "Intra_noFac",
  "Intra_withFac",
  "Inter_noFac",
  "Inter_withFac"
)

pred_grid <- expand_grid(
  focsp     = factor(c("TRCY", "TROR"), levels = c("TRCY", "TROR")),
  Site      = factor(site_levels, levels = site_levels),
  Cov       = factor(c("Shade", "Sun"), levels = c("Shade", "Sun")),
  condition = factor(cond_levels, levels = cond_levels)
) %>%
  left_join(q90_tbl, by = c("focsp", "Site")) %>%
  mutate(
    fac_sc = case_when(
      condition %in% c("Lambda", "Intra_noFac", "Inter_noFac") ~ 0,
      condition %in% c("Intra_withFac", "Inter_withFac") ~ fac_90
    ),
    intra_sc = case_when(
      condition %in% c("Intra_noFac", "Intra_withFac") ~ intra_90,
      TRUE ~ 0
    ),
    inter_sc = case_when(
      condition %in% c("Inter_noFac", "Inter_withFac") ~ inter_90,
      TRUE ~ 0
    ),
    fac_z = (fac_sc - fac_center) / fac_scale,
    fac_z2 = fac_z^2,
    BLK = levels(dat$BLK)[1],
    comp_type = case_when(
      condition == "Lambda" ~ "Lambda",
      str_detect(condition, "^Intra") ~ "Intra",
      str_detect(condition, "^Inter") ~ "Inter"
    ),
    fac_state = case_when(
      condition == "Lambda" ~ "Lambda",
      str_detect(condition, "withFac") ~ "With facilitator",
      TRUE ~ "No facilitator"
    )
  )

############################
## 5. Predict mean fecundity on response scale
############################
pred_sum <- fitted(
  best_fit,
  newdata = pred_grid,
  re_formula = NA,
  summary = TRUE
)

plot_dat <- bind_cols(pred_grid, as.data.frame(pred_sum))

############################
## 6. Symmetric positions: 5 in Shade, 5 in Sun
############################
offset_key <- tribble(
  ~Cov,     ~condition,         ~offset,
  "Shade",  "Lambda",           -0.36,
  "Shade",  "Intra_noFac",      -0.28,
  "Shade",  "Intra_withFac",    -0.20,
  "Shade",  "Inter_noFac",      -0.12,
  "Shade",  "Inter_withFac",    -0.04,
  "Sun",    "Lambda",            0.04,
  "Sun",    "Intra_noFac",       0.12,
  "Sun",    "Intra_withFac",     0.20,
  "Sun",    "Inter_noFac",       0.28,
  "Sun",    "Inter_withFac",     0.36
) %>%
  mutate(
    Cov = as.character(Cov),
    condition = as.character(condition),
    offset = as.numeric(offset)
  )

plot_dat <- plot_dat %>%
  mutate(
    Cov = as.character(Cov),
    condition = as.character(condition)
  ) %>%
  left_join(offset_key, by = c("Cov", "condition")) %>%
  mutate(
    site_id = match(as.character(Site), site_levels),
    xpos = as.numeric(site_id) + as.numeric(offset),
    fill_key = case_when(
      comp_type == "Lambda" ~ paste0("lambda_", Site),
      fac_state == "With facilitator" ~ as.character(Site),
      TRUE ~ "open"
    )
  )

############################
## 7. Shade background
############################
shade_rects <- tibble(
  site_id = seq_along(site_levels),
  xmin = site_id - 0.40,
  xmax = site_id,
  ymin = -Inf,
  ymax = Inf
)

############################
## 8. Two-panel plotting helper
############################
make_species_plot <- function(species_code,
                              panel_title,
                              show_x_axis = FALSE,
                              show_y_title = TRUE) {
  
  dat_sp <- plot_dat %>%
    filter(focsp == species_code)
  
  ggplot() +
    geom_rect(
      data = shade_rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey93",
      colour = NA
    ) +
    geom_linerange(
      data = dat_sp,
      aes(x = xpos, ymin = Q2.5, ymax = Q97.5, colour = Site),
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    ## non-Lambda points
    geom_point(
      data = dat_sp %>% filter(comp_type != "Lambda"),
      aes(
        x = xpos,
        y = Estimate,
        shape = comp_type,
        colour = Site,
        fill = fill_key
      ),
      size = 2.5,   # ~20% smaller than before
      stroke = 1,
      alpha = 1
    ) +
    ## Lambda points: darker site colour, semi-transparent, no black outline
    geom_point(
      data = dat_sp %>% filter(comp_type == "Lambda"),
      aes(
        x = xpos,
        y = Estimate,
        shape = comp_type,
        colour = Site,
        fill = fill_key
      ),
      size = 2.5,
      stroke = 1,
      alpha = 1,
      color="black"
    ) +
    scale_x_continuous(
      breaks = seq_along(site_levels),
      labels = site_levels,
      limits = c(0.58, length(site_levels) + 0.42),
      expand = c(0, 0)
    ) +
    scale_colour_manual(values = site_colors, guide = "none") +
    scale_fill_manual(
      values = c(
        "open" = "white",
        site_colors,
        setNames(lambda_colors, paste0("lambda_", names(lambda_colors)))
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "Lambda" = 22,
        "Intra"  = 21,
        "Inter"  = 24
      ),
      guide = guide_legend(
        title = "Competition type",
        override.aes = list(
          fill = "white",
          colour = "black",
          alpha = 1,
          size = 3
        )
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.x = if (show_x_axis) {
        element_text(size = 16, colour = "black")
      } else {
        element_blank()
      },
      axis.text.x  = if (show_x_axis) {
        element_text(size = 13, margin = margin(t = 10), colour = "black")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) element_line(colour = "black") else element_blank(),
      axis.title.y = if (show_y_title) {
        element_text(size = 16, colour = "black")
      } else {
        element_blank()
      },
      axis.text.y  = element_text(size = 12, colour = "black"),
      axis.ticks.y = element_line(colour = "black"),
      plot.title   = element_text(size = 16, hjust = 0, face = "plain", colour = "black"),
      plot.margin  = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = panel_title,
      x = if (show_x_axis) "Site" else NULL,
      y = if (show_y_title) "Mean Seed Production" else NULL
    )
}

############################
## 9. Make panels
############################
p_trcy <- make_species_plot(
  species_code = "TRCY",
  panel_title = expression("(a) " * italic("T. cyanopetala")),
  show_x_axis = FALSE,
  show_y_title = TRUE
)

p_tror <- make_species_plot(
  species_code = "TROR",
  panel_title = expression("(b) " * italic("T. ornata")),
  show_x_axis = TRUE,
  show_y_title = TRUE
)

p_final <- (p_trcy / p_tror) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right")

print(p_final)

ggsave(
  file.path(run_dir, "Figure_predicted_fecundity_lambda_competition.png"),
  p_final,
  width = 12,
  height = 7,
  dpi = 300
)

############################################################
## COEFFICIENT PLOT
## species-specific effects with Ben added explicitly
############################################################

library(tidyverse)
library(brms)
library(posterior)

############################
## 1. Load final saved model
############################
analysis_id <- "env_fac_model_comparison_v1"
run_dir     <- file.path(getwd(), analysis_id)

load_final_glmm <- function(run_dir,
                            label = "with_CS",
                            model_set = c("fac_context", "broad")) {
  
  model_set <- match.arg(model_set)
  
  if (model_set == "fac_context") {
    out_dir <- file.path(run_dir, paste0("fac_context_test_", label))
    best_file <- file.path(out_dir, paste0("best_fac_context_model_", label, ".txt"))
    
    fit_file_map <- c(
      nofac     = paste0("fit_site_cov_nofac_", label, ".rds"),
      linear    = paste0("fit_site_cov_fac_linear_", label, ".rds"),
      quad      = paste0("fit_site_cov_fac_quad_", label, ".rds"),
      site_quad = paste0("fit_site_quadfac_", label, ".rds"),
      cov_quad  = paste0("fit_cov_quadfac_", label, ".rds"),
      both_quad = paste0("fit_both_quadfac_", label, ".rds")
    )
    
  } else {
    out_dir <- file.path(run_dir, paste0("models_", label))
    best_file <- file.path(out_dir, paste0("best_model_", label, ".txt"))
    
    fit_file_map <- c(
      base     = paste0("fit_base_", label, ".rds"),
      cov      = paste0("fit_cov_", label, ".rds"),
      site     = paste0("fit_site_", label, ".rds"),
      site_cov = paste0("fit_site_cov_", label, ".rds"),
      covfac   = paste0("fit_covfac_", label, ".rds"),
      sitefac  = paste0("fit_sitefac_", label, ".rds"),
      bothfac  = paste0("fit_bothfac_", label, ".rds")
    )
  }
  
  best_name <- scan(best_file, what = "character", quiet = TRUE)
  fit_path  <- file.path(out_dir, fit_file_map[[best_name]])
  fit       <- readRDS(fit_path)
  
  list(best_name = best_name, fit = fit, out_dir = out_dir)
}

final_glmm <- load_final_glmm(
  run_dir = run_dir,
  label = "with_CS",
  model_set = "fac_context"
)

best_fit_final <- final_glmm$fit

############################
## 2. Posterior draws helper
############################
draws <- posterior::as_draws_df(best_fit_final)

b <- function(term) {
  nm <- paste0("b_", term)
  if (!nm %in% names(draws)) stop("Missing coefficient: ", nm)
  draws[[nm]]
}

summ_term <- function(x, species, variable) {
  tibble(
    species  = species,
    variable = variable,
    Estimate = mean(x),
    lower    = quantile(x, 0.025),
    upper    = quantile(x, 0.975)
  )
}

############################
## 3. Build table
##
## Ben is added explicitly from the species baselines
############################
coef_dat <- bind_rows(
  
  ## ---------------- Site rows ----------------
  summ_term(b("focspTRCY"), "TRCY", "SiteBen"),
  summ_term(b("focspTRCY") + b("focspTRCY:SiteGH"), "TRCY", "SiteGH"),
  summ_term(b("focspTRCY") + b("focspTRCY:SiteNam"), "TRCY", "SiteNam"),
  summ_term(b("focspTRCY") + b("focspTRCY:SitePJ"), "TRCY", "SitePJ"),
  summ_term(b("focspTRCY") + b("focspTRCY:SiteCS"), "TRCY", "SiteCS"),
  
  summ_term(b("focspTROR"), "TROR", "SiteBen"),
  summ_term(b("focspTROR") + b("focspTROR:SiteGH"), "TROR", "SiteGH"),
  summ_term(b("focspTROR") + b("focspTROR:SiteNam"), "TROR", "SiteNam"),
  summ_term(b("focspTROR") + b("focspTROR:SitePJ"), "TROR", "SitePJ"),
  summ_term(b("focspTROR") + b("focspTROR:SiteCS"), "TROR", "SiteCS"),
  
  ## ---------------- Other rows ----------------
  summ_term(b("focspTRCY:CovSun"),   "TRCY", "CovSun"),
  summ_term(b("focspTROR:CovSun"),   "TROR", "CovSun"),
  
  summ_term(b("focspTRCY:fac_z"),    "TRCY", "fac_z"),
  summ_term(b("focspTROR:fac_z"),    "TROR", "fac_z"),
  
  summ_term(b("focspTRCY:fac_z2"),   "TRCY", "fac_z2"),
  summ_term(b("focspTROR:fac_z2"),   "TROR", "fac_z2"),
  
  summ_term(b("focspTRCY:intra_sc"), "TRCY", "intra_sc"),
  summ_term(b("focspTROR:intra_sc"), "TROR", "intra_sc"),
  
  summ_term(b("focspTRCY:inter_sc"), "TRCY", "inter_sc"),
  summ_term(b("focspTROR:inter_sc"), "TROR", "inter_sc")
)

############################
## 4. Nice labels
############################
label_map <- c(
  SiteBen  = "Site = Ben",
  SiteGH   = "Site = GH",
  SiteNam  = "Site = Nam",
  SitePJ   = "Site = PJ",
  SiteCS   = "Site = CS",
  CovSun   = "Cover = Sun",
  fac_z    = "Facilitator (linear)",
  fac_z2   = "Facilitator²",
  intra_sc = "Intraspecific density",
  inter_sc = "Interspecific density"
)

var_order <- c(
  "SiteBen",
  "SiteGH",
  "SiteNam",
  "SitePJ",
  "SiteCS",
  "CovSun",
  "fac_z",
  "fac_z2",
  "intra_sc",
  "inter_sc"
)

coef_dat <- coef_dat %>%
  mutate(
    variable = factor(variable, levels = rev(var_order)),
    variable_label = factor(label_map[as.character(variable)],
                            levels = rev(label_map[var_order])),
    species = factor(species, levels = c("TRCY", "TROR"))
  )

############################
## 5. Colours
############################
species_cols <- c(
  "TRCY" = "#065143",
  "TROR" = "#F2C14E"
)

############################
## 6. Plot
############################
pd <- position_dodge(width = 0.55)

p_coef_james <- ggplot(coef_dat, aes(x = Estimate, y = variable_label, colour = species)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0,
    linewidth = 0.9,
    position = pd
  ) +
  geom_point(
    size = 3.2,
    position = pd
  ) +
  scale_colour_manual(values = species_cols) +
  coord_cartesian(xlim = c(-3, 3), clip = "off") +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, colour = "black"),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    axis.text.x  = element_text(size = 14, colour = "black"),
    axis.text.y  = element_text(size = 13, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13, colour = "black")
  ) +
  labs(
    title = paste("Final fac_context GLMM:", final_glmm$best_name),
    x = "Posterior estimate",
    y = "Variable"
  )

print(p_coef_james)

ggsave(
  file.path(final_glmm$out_dir, "final_glmm_coefficient_plot.png"),
  p_coef_james,
  width = 8,
  height = 6.5,
  dpi = 300
)

write.csv(
  coef_dat,
  file.path(final_glmm$out_dir, "final_glmm_coefficient_plot.csv"),
  row.names = FALSE
)


############################################################
## RESPONSE SURFACE FIGURES FROM FINAL GLMM
## 5 rows (sites) x 4 cols (Shade/Sun x con/hetero)
## Separate figure for TRCY and TROR
## Uses observed site-specific ranges and overlays raw data
############################################################

rm(list = ls())

############################
## 0. Packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "viridis"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(brms)
library(viridis)

############################
## 1. Load final saved GLMM
############################
analysis_id <- "env_fac_model_comparison_v1"
run_dir     <- file.path(getwd(), analysis_id)

load_final_glmm <- function(run_dir,
                            label = "with_CS",
                            model_set = c("fac_context", "broad")) {
  
  model_set <- match.arg(model_set)
  
  if (model_set == "fac_context") {
    out_dir <- file.path(run_dir, paste0("fac_context_test_", label))
    best_file <- file.path(out_dir, paste0("best_fac_context_model_", label, ".txt"))
    
    fit_file_map <- c(
      nofac     = paste0("fit_site_cov_nofac_", label, ".rds"),
      linear    = paste0("fit_site_cov_fac_linear_", label, ".rds"),
      quad      = paste0("fit_site_cov_fac_quad_", label, ".rds"),
      site_quad = paste0("fit_site_quadfac_", label, ".rds"),
      cov_quad  = paste0("fit_cov_quadfac_", label, ".rds"),
      both_quad = paste0("fit_both_quadfac_", label, ".rds")
    )
    
  } else {
    out_dir <- file.path(run_dir, paste0("models_", label))
    best_file <- file.path(out_dir, paste0("best_model_", label, ".txt"))
    
    fit_file_map <- c(
      base     = paste0("fit_base_", label, ".rds"),
      cov      = paste0("fit_cov_", label, ".rds"),
      site     = paste0("fit_site_", label, ".rds"),
      site_cov = paste0("fit_site_cov_", label, ".rds"),
      covfac   = paste0("fit_covfac_", label, ".rds"),
      sitefac  = paste0("fit_sitefac_", label, ".rds"),
      bothfac  = paste0("fit_bothfac_", label, ".rds")
    )
  }
  
  if (!file.exists(best_file)) {
    stop("Cannot find best-model file: ", best_file)
  }
  
  best_name <- scan(best_file, what = "character", quiet = TRUE)
  
  if (!best_name %in% names(fit_file_map)) {
    stop("Best model name not recognised: ", best_name)
  }
  
  fit_path <- file.path(out_dir, fit_file_map[[best_name]])
  
  if (!file.exists(fit_path)) {
    stop("Cannot find saved model file: ", fit_path)
  }
  
  fit <- readRDS(fit_path)
  
  list(
    best_name = best_name,
    fit = fit,
    out_dir = out_dir
  )
}

final_glmm <- load_final_glmm(
  run_dir = run_dir,
  label = "with_CS",
  model_set = "fac_context"
)

best_fit_final <- final_glmm$fit
cat("Loaded final GLMM:", final_glmm$best_name, "\n")

############################
## 2. Load fecundity data
############################
dat_fec <- read.csv(file.path(run_dir, "fecundity_analysis_data_all_sites.csv"))

site_levels <- c("Ben", "GH", "Nam", "PJ", "CS")
cov_levels  <- c("Shade", "Sun")
sp_levels   <- c("TRCY", "TROR")

dat_fec <- dat_fec %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    Cov   = factor(Cov, levels = cov_levels),
    focsp = factor(focsp, levels = sp_levels),
    BLK   = factor(BLK)
  )

fac_center <- mean(dat_fec$fac_sc, na.rm = TRUE)
fac_scale  <- sd(dat_fec$fac_sc, na.rm = TRUE)

if (!is.finite(fac_scale) || fac_scale == 0) {
  stop("fac_sc has zero or undefined SD; cannot compute fac_z.")
}

############################
## 3. Helpers
############################
safe_seq <- function(x, n = 70) {
  rng <- range(x, na.rm = TRUE)
  
  if (!all(is.finite(rng))) {
    return(rep(0, n))
  }
  
  if (rng[1] == rng[2]) {
    eps <- ifelse(rng[1] == 0, 0.1, abs(rng[1]) * 0.05)
    rng <- c(rng[1] - eps, rng[2] + eps)
  }
  
  seq(rng[1], rng[2], length.out = n)
}

############################
## 4. Build response-surface data
##
## For each panel:
## - x-axis = observed facilitator range
## - y-axis = observed density range for the focal panel
## - unused density held at observed median for that site/cov/species
############################
if (!requireNamespace("ggh4x", quietly = TRUE)) {
  install.packages("ggh4x")
}
library(ggh4x)

make_surface_data <- function(species_code,
                              n_fac = 70,
                              n_dens = 70) {
  
  raw_sp <- dat_fec %>%
    filter(focsp == species_code)
  
  panel_defs <- tribble(
    ~Cov,     ~density_type,      ~yvar,       ~other_yvar,   ~panel,
    "Shade",  "Conspecific",      "intra_sc",  "inter_sc",    "Shade\nconspecific",
    "Shade",  "Heterospecific",   "inter_sc",  "intra_sc",    "Shade\nheterospecific",
    "Sun",    "Conspecific",      "intra_sc",  "inter_sc",    "Sun\nconspecific",
    "Sun",    "Heterospecific",   "inter_sc",  "intra_sc",    "Sun\nheterospecific"
  ) %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          "Shade\nconspecific",
          "Shade\nheterospecific",
          "Sun\nconspecific",
          "Sun\nheterospecific"
        )
      )
    )
  
  grid_list <- list()
  
  for (i in seq_len(nrow(panel_defs))) {
    this_cov    <- panel_defs$Cov[i]
    this_type   <- panel_defs$density_type[i]
    this_yvar   <- panel_defs$yvar[i]
    this_other  <- panel_defs$other_yvar[i]
    this_panel  <- panel_defs$panel[i]
    
    for (this_site in site_levels) {
      
      dat_sub <- raw_sp %>%
        filter(Site == this_site, Cov == this_cov)
      
      if (nrow(dat_sub) == 0) next
      
      fac_seq  <- safe_seq(dat_sub$fac_sc, n = n_fac)
      dens_seq <- safe_seq(dat_sub[[this_yvar]], n = n_dens)
      
      other_fixed <- median(dat_sub[[this_other]], na.rm = TRUE)
      
      g <- expand_grid(
        focsp    = factor(species_code, levels = sp_levels),
        Site     = factor(this_site, levels = site_levels),
        Cov      = factor(this_cov, levels = cov_levels),
        fac_sc   = fac_seq,
        dens_val = dens_seq
      )
      
      if (this_yvar == "intra_sc") {
        g <- g %>%
          mutate(
            intra_sc = dens_val,
            inter_sc = other_fixed
          )
      } else {
        g <- g %>%
          mutate(
            intra_sc = other_fixed,
            inter_sc = dens_val
          )
      }
      
      g <- g %>%
        mutate(
          fac_z  = (fac_sc - fac_center) / fac_scale,
          fac_z2 = fac_z^2,
          BLK = levels(dat_fec$BLK)[1],
          density_type = this_type,
          other_fixed = other_fixed,
          panel = this_panel
        )
      
      grid_list[[length(grid_list) + 1]] <- g
    }
  }
  
  pred_grid <- bind_rows(grid_list)
  
  pred <- fitted(
    best_fit_final,
    newdata = pred_grid,
    re_formula = NA,
    summary = TRUE
  )
  
  bind_cols(pred_grid, as.data.frame(pred))
}

############################
## 5. Build raw-data overlay
##
## Each observation appears in both its conspecific panel
## and its heterospecific panel
############################
make_surface_obs <- function(species_code) {
  
  raw_sp <- dat_fec %>%
    filter(focsp == species_code) %>%
    mutate(
      Site = factor(Site, levels = site_levels),
      Cov  = factor(Cov, levels = cov_levels)
    )
  
  obs_cons <- raw_sp %>%
    transmute(
      focsp,
      Site,
      Cov,
      fac_sc,
      dens_val = intra_sc,
      fitness,
      panel = case_when(
        Cov == "Shade" ~ "Shade\nconspecific",
        Cov == "Sun"   ~ "Sun\nconspecific"
      )
    )
  
  obs_hetero <- raw_sp %>%
    transmute(
      focsp,
      Site,
      Cov,
      fac_sc,
      dens_val = inter_sc,
      fitness,
      panel = case_when(
        Cov == "Shade" ~ "Shade\nheterospecific",
        Cov == "Sun"   ~ "Sun\nheterospecific"
      )
    )
  
  bind_rows(obs_cons, obs_hetero) %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          "Shade\nconspecific",
          "Shade\nheterospecific",
          "Sun\nconspecific",
          "Sun\nheterospecific"
        )
      )
    )
}

############################
## 6. Plot helper
############################
plot_surface_figure <- function(surface_dat,
                                obs_dat,
                                species_code,
                                fill_limits = NULL) {
  
  title_lab <- if (species_code == "TRCY") {
    expression(italic("T. cyanopetala"))
  } else {
    expression(italic("T. ornata"))
  }
  
  y_axis_lab <- if (species_code == "TRCY") {
    "Density (conspecific = TRCY; heterospecific = TROR)"
  } else {
    "Density (conspecific = TROR; heterospecific = TRCY)"
  }
  
  p <- ggplot(surface_dat, aes(x = fac_sc, y = dens_val, fill = Estimate)) +
    geom_tile() +
    geom_contour(
      aes(z = Estimate),
      colour = "white",
      alpha = 0.45,
      linewidth = 0.2
    ) +
    geom_point(
      data = obs_dat,
      aes(x = fac_sc, y = dens_val, size = fitness),
      inherit.aes = FALSE,
      shape = 1,
      colour = "grey20",
      stroke = 0.3,
      alpha = 0.8
    ) +
    ggh4x::facet_grid2(
      Site ~ panel,
      scales = "free",
      independent = "all"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_size_continuous(
      name = "Observed\nseed number",
      range = c(0.6, 3.2)
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.text.x = element_text(size = 13),
      strip.text.y = element_text(size = 13),
      axis.title = element_text(size = 16, colour = "black"),
      axis.text = element_text(size = 11, colour = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 18, colour = "black")
    ) +
    labs(
      title = title_lab,
      x = "Facilitator density",
      y = y_axis_lab,
      fill = "Predicted\nfecundity"
    )
  
  if (is.null(fill_limits)) {
    p <- p + scale_fill_viridis_c(option = "C")
  } else {
    p <- p + scale_fill_viridis_c(option = "C", limits = fill_limits)
  }
  
  p
}
############################
## 7. Build surfaces and overlays
############################
surf_trcy <- make_surface_data("TRCY", n_fac = 70, n_dens = 70)
surf_tror <- make_surface_data("TROR", n_fac = 70, n_dens = 70)

obs_trcy <- make_surface_obs("TRCY")
obs_tror <- make_surface_obs("TROR")

############################
## 8. Common fill scale across both species
############################
common_fill_limits <- range(
  c(surf_trcy$Estimate, surf_tror$Estimate),
  na.rm = TRUE
)

############################
## 9. Plot
############################
p_surf_trcy <- plot_surface_figure(
  surface_dat = surf_trcy,
  obs_dat = obs_trcy,
  species_code = "TRCY",
  fill_limits = common_fill_limits
)

p_surf_tror <- plot_surface_figure(
  surface_dat = surf_tror,
  obs_dat = obs_tror,
  species_code = "TROR",
  fill_limits = common_fill_limits
)

print(p_surf_trcy)
print(p_surf_tror)

############################
## 10. Save outputs
############################
ggsave(
  file.path(final_glmm$out_dir, "response_surface_TRCY_obs_overlay.png"),
  p_surf_trcy,
  width = 12,
  height = 13,
  dpi = 300
)

ggsave(
  file.path(final_glmm$out_dir, "response_surface_TROR_obs_overlay.png"),
  p_surf_tror,
  width = 12,
  height = 13,
  dpi = 300
)

write.csv(
  surf_trcy,
  file.path(final_glmm$out_dir, "response_surface_TRCY_obs_overlay_data.csv"),
  row.names = FALSE
)

write.csv(
  surf_tror,
  file.path(final_glmm$out_dir, "response_surface_TROR_obs_overlay_data.csv"),
  row.names = FALSE
)

cat("Saved figures to:\n")
cat(file.path(final_glmm$out_dir, "response_surface_TRCY_obs_overlay.png"), "\n")
cat(file.path(final_glmm$out_dir, "response_surface_TROR_obs_overlay.png"), "\n")