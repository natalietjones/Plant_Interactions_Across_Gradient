############################################################
## 03_invasion_growth_by_site_cover_fac_q50_sitefill.R
##
## Updates:
## - evaluates "without facilitator" at fac_sc = 0
## - evaluates "with facilitator" as the mean fac_sc among
##   observations at or above the 50th percentile (within Site x Cov),
##   restricting to positive facilitator densities when available
## - panels split by Cover (Shade / Sun)
## - x-order within each Site is:
##     without facilitator: TRCY, TROR
##     with facilitator:    TRCY, TROR
## - shape encodes focal species (TRCY circle, TROR triangle)
## - filled vs open encodes facilitator treatment
## - point fill for "with facilitator" uses the site colour
## - removes grey background bands
############################################################

rm(list = ls())

required_packages <- c("tidyverse", "brms", "posterior")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(brms)
library(posterior)

options(mc.cores = max(1, parallel::detectCores() - 1))
options(brms.backend = "rstan")

############################
## 1. User options
############################
use_germination <- FALSE
structure_dir   <- "01_alpha_structure_selection"
out_dir         <- "03_invasion_growth"

site_cols <- c(
  "BEN" = "#6F2DBD",
  "GH"  = "#5E81D8",
  "NAM" = "#44C07A",
  "PJ"  = "#F59E00",
  "CS"  = "#E31A1C"
)

site_labels <- c(
  "BEN" = "Ben",
  "GH"  = "GH",
  "NAM" = "Nam",
  "PJ"  = "PJ",
  "CS"  = "CS"
)

############################
## 2. Output folders
############################
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "draws"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)

############################
## 3. Load selected model and data
############################
alpha_info <- readRDS(file.path(structure_dir, "selection", "alpha_structure_selection_info.rds"))

if (!identical(alpha_info$selected_growth_family, "BH")) {
  stop("This script currently supports BH only. Selected family: ", alpha_info$selected_growth_family)
}

dat_stage1 <- readRDS(alpha_info$stage1_data_rds)
fit_final  <- readRDS(alpha_info$selected_fit_file)
site_order <- alpha_info$site_order

family_tag <- alpha_info$selected_growth_family
error_tag  <- ifelse(alpha_info$selected_error_model == "gaussian_logfit",
                     "loggaussian",
                     alpha_info$selected_error_model)
prefix_tag <- ifelse(use_germination, "germ", "nogerm")
file_prefix <- paste0(family_tag, "_", error_tag, "_", prefix_tag, "_igr_q50_sitefill")

############################
## 4. Helper: matrix -> long data frame
############################
mat_to_long <- function(mat, meta, value_name) {
  meta <- meta %>% mutate(row_id = row_number())
  colnames(mat) <- paste0("V", seq_len(ncol(mat)))

  as_tibble(mat, .name_repair = "minimal") %>%
    mutate(draw = row_number()) %>%
    pivot_longer(
      cols = -draw,
      names_to = "col",
      values_to = value_name
    ) %>%
    mutate(row_id = as.integer(sub("^V", "", col))) %>%
    select(-col) %>%
    left_join(meta, by = "row_id") %>%
    select(draw, everything(), -row_id)
}

############################
## 5. Select environments to evaluate
############################
selected_env_zero <- dat_stage1 %>%
  distinct(Site, Cov) %>%
  mutate(
    fac_sc = 0,
    fac_label = "Without facilitator",
    fac_source = "fixed_zero"
  )

selected_env_with <- dat_stage1 %>%
  group_by(Site, Cov) %>%
  summarise(
    fac_q50 = quantile(fac_sc, 0.50, na.rm = TRUE),
    fac_sc = {
      x <- fac_sc
      x_hi_pos <- x[x >= fac_q50 & x > 0]
      if (length(x_hi_pos) > 0) {
        mean(x_hi_pos, na.rm = TRUE)
      } else {
        x_hi <- x[x >= fac_q50]
        if (length(x_hi) > 0) mean(x_hi, na.rm = TRUE) else 0
      }
    },
    n_used = {
      x <- fac_sc
      x_hi_pos <- x[x >= fac_q50 & x > 0]
      if (length(x_hi_pos) > 0) length(x_hi_pos) else sum(x >= fac_q50, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  mutate(
    fac_label = "With facilitator",
    fac_source = "mean_ge_q50"
  )

selected_env <- bind_rows(selected_env_zero, selected_env_with) %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov = factor(as.character(Cov), levels = c("SH", "SUN")),
    fac_label = factor(fac_label, levels = c("Without facilitator", "With facilitator"))
  )

############################
## 6. Prediction grids
############################
lambda_grid <- expand_grid(
  focsp = levels(dat_stage1$focsp),
  selected_env %>% select(Site, Cov, fac_sc, fac_label)
) %>%
  mutate(
    BLK      = levels(dat_stage1$BLK)[1],
    intra_sc = 0,
    inter_sc = 0
  )

alpha_grid <- expand_grid(
  focsp = levels(dat_stage1$focsp),
  Site  = levels(dat_stage1$Site),
  Cov   = levels(dat_stage1$Cov)
) %>%
  mutate(
    fac_sc   = 0,
    BLK      = levels(dat_stage1$BLK)[1],
    intra_sc = 0,
    inter_sc = 0
  )

############################
## 7. Posterior draws of nonlinear predictors
############################
etalam_draws <- posterior_linpred(
  fit_final,
  newdata    = lambda_grid,
  nlpar      = "etalam",
  transform  = FALSE,
  re_formula = NA
)

etaai_draws <- posterior_linpred(
  fit_final,
  newdata    = alpha_grid,
  nlpar      = "etaai",
  transform  = FALSE,
  re_formula = NA
)

etaaj_draws <- posterior_linpred(
  fit_final,
  newdata    = alpha_grid,
  nlpar      = "etaaj",
  transform  = FALSE,
  re_formula = NA
)

############################
## 8. Natural-scale lambda and alpha
############################
lambda_draws <- mat_to_long(etalam_draws, lambda_grid, "etalam") %>%
  mutate(lambda = exp(etalam)) %>%
  select(draw, focsp, Site, Cov, fac_sc, fac_label, lambda)

alpha_ii_draws <- mat_to_long(etaai_draws, alpha_grid, "etaai") %>%
  mutate(alpha_ii = exp(etaai)) %>%
  select(draw, focsp, Site, Cov, alpha_ii)

alpha_ij_draws <- mat_to_long(etaaj_draws, alpha_grid, "etaaj") %>%
  mutate(alpha_ij = exp(etaaj)) %>%
  select(draw, focsp, Site, Cov, alpha_ij)

alpha_draws <- alpha_ii_draws %>%
  left_join(alpha_ij_draws, by = c("draw", "focsp", "Site", "Cov"))

############################
## 9. Germination draws or constants
############################
n_draws <- n_distinct(lambda_draws$draw)

if (use_germination) {
  if (file.exists(alpha_info$germination_site_priors_rds)) {
    germ_priors <- readRDS(alpha_info$germination_site_priors_rds) %>% as_tibble()
  } else {
    germ_info <- readRDS(alpha_info$germination_priors_rds)
    germ_priors <- germ_info$site_priors %>% as_tibble()
  }

  germ_priors <- germ_priors %>%
    mutate(
      Site = as.character(Site),
      focsp = as.character(focsp)
    )

  set.seed(123)
  germ_draws <- expand_grid(
    draw = seq_len(n_draws),
    germ_priors %>% select(Site, focsp, g_alpha_final, g_beta_final, g_mean_final)
  ) %>%
    mutate(g_draw = rbeta(n = n(), shape1 = g_alpha_final, shape2 = g_beta_final))

  germ_wide <- germ_draws %>%
    select(draw, Site, focsp, g_draw) %>%
    pivot_wider(
      names_from = focsp,
      values_from = g_draw,
      names_glue = "g_{focsp}"
    )
} else {
  germ_wide <- expand_grid(
    draw = seq_len(n_draws),
    Site = levels(dat_stage1$Site)
  ) %>%
    mutate(
      g_TRCY = 1,
      g_TROR = 1
    )
}

############################
## 10. Wide draws object
############################
lambda_wide <- lambda_draws %>%
  pivot_wider(
    names_from  = focsp,
    values_from = lambda,
    names_glue  = "lambda_{focsp}"
  )

alpha_wide <- alpha_draws %>%
  pivot_wider(
    names_from  = focsp,
    values_from = c(alpha_ii, alpha_ij),
    names_glue  = "{.value}_{focsp}"
  )

igr_draws <- lambda_wide %>%
  left_join(alpha_wide, by = c("draw", "Site", "Cov")) %>%
  left_join(germ_wide, by = c("draw", "Site"))

############################
## 11. Eta terms and resident equilibrium densities
############################
s_vals <- c(TRCY = 0, TROR = 0)

igr_draws <- igr_draws %>%
  mutate(
    eta_TRCY = lambda_TRCY * g_TRCY / (1 - s_vals["TRCY"] * (1 - g_TRCY)),
    eta_TROR = lambda_TROR * g_TROR / (1 - s_vals["TROR"] * (1 - g_TROR)),
    Nstar_TRCY = ifelse(eta_TRCY > 1, (eta_TRCY - 1) / (alpha_ii_TRCY * g_TRCY), 0),
    Nstar_TROR = ifelse(eta_TROR > 1, (eta_TROR - 1) / (alpha_ii_TROR * g_TROR), 0)
  )

############################
## 12. Invasion growth rates when rare
############################
igr_draws <- igr_draws %>%
  mutate(
    growth_factor_TRCY_invades =
      (1 - g_TRCY) * s_vals["TRCY"] +
      (g_TRCY * lambda_TRCY) / (1 + alpha_ij_TRCY * g_TROR * Nstar_TROR),

    growth_factor_TROR_invades =
      (1 - g_TROR) * s_vals["TROR"] +
      (g_TROR * lambda_TROR) / (1 + alpha_ij_TROR * g_TRCY * Nstar_TRCY),

    igr_TRCY = log(growth_factor_TRCY_invades),
    igr_TROR = log(growth_factor_TROR_invades)
  )

igr_long <- igr_draws %>%
  select(draw, Site, Cov, fac_sc, fac_label, igr_TRCY, igr_TROR) %>%
  pivot_longer(
    cols = c(igr_TRCY, igr_TROR),
    names_to = "species_code",
    values_to = "igr"
  ) %>%
  mutate(
    species_code = recode(species_code, igr_TRCY = "TRCY", igr_TROR = "TROR"),
    species = recode(species_code, TRCY = "T. cyanopetala", TROR = "T. ornata"),
    species = factor(species, levels = c("T. cyanopetala", "T. ornata")),
    species_code = factor(species_code, levels = c("TRCY", "TROR")),
    Site = factor(as.character(Site), levels = site_order),
    Cov = factor(as.character(Cov), levels = c("SH", "SUN"))
  ) %>%
  filter(is.finite(igr))

############################
## 13. Summaries
############################
x_levels <- unlist(lapply(site_order, function(st) {
  c(
    paste(st, "Without facilitator", "TRCY", sep = "__"),
    paste(st, "Without facilitator", "TROR", sep = "__"),
    paste(st, "With facilitator", "TRCY", sep = "__"),
    paste(st, "With facilitator", "TROR", sep = "__")
  )
}))

igr_summary <- igr_long %>%
  group_by(Site, Cov, fac_label, fac_sc, species, species_code) %>%
  summarise(
    igr_med = median(igr),
    igr_lwr = quantile(igr, 0.025),
    igr_upr = quantile(igr, 0.975),
    p_positive = mean(igr > 0),
    .groups = "drop"
  ) %>%
  mutate(
    x_group = paste(as.character(Site), as.character(fac_label), as.character(species_code), sep = "__"),
    x_group = factor(x_group, levels = x_levels),
    site_colour = site_cols[as.character(Site)],
    site_label = site_labels[as.character(Site)],
    fill_value = ifelse(fac_label == "With facilitator", site_colour, "white")
  )

write.csv(
  igr_summary,
  file.path(out_dir, "tables", paste0(file_prefix, "_summary.csv")),
  row.names = FALSE
)

saveRDS(
  igr_long,
  file.path(out_dir, "draws", paste0(file_prefix, "_draws.rds"))
)

############################
## 14. Plot
############################

cover_labels <- c("SH" = "Shade", "SUN" = "Sun")
shape_vals <- c("TRCY" = 21, "TROR" = 24)

site_centers <- tibble(
  Site = factor(site_order, levels = site_order),
  site_label = site_labels[site_order],
  x = seq(2.5, by = 4, length.out = length(site_order))
)

site_boundaries <- tibble(xint = seq(4.5, by = 4, length.out = length(site_order) - 1))

shape_legend_df <- tibble(
  species_code = factor(c("TRCY", "TROR"), levels = c("TRCY", "TROR")),
  x_group = factor(x_levels[1], levels = x_levels),
  igr_med = NA_real_
)

fill_legend_df <- tibble(
  fac_label = factor(c("Without facilitator", "With facilitator"),
                     levels = c("Without facilitator", "With facilitator")),
  x_group = factor(x_levels[1], levels = x_levels),
  igr_med = NA_real_
)

p_igr <- ggplot(igr_summary, aes(x = x_group, y = igr_med)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey35") +
  geom_vline(
    data = site_boundaries,
    aes(xintercept = xint),
    inherit.aes = FALSE,
    linewidth = 0.3,
    colour = "grey85"
  ) +
  geom_errorbar(
    aes(ymin = igr_lwr, ymax = igr_upr, colour = Site),
    width = 0.0,
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = igr_summary %>% filter(fac_label == "Without facilitator"),
    aes(shape = species_code, colour = Site),
    fill = "white",
    size = 3.1,
    stroke = 1,
    show.legend = FALSE
  ) +
  geom_point(
    data = igr_summary %>% filter(fac_label == "With facilitator"),
    aes(shape = species_code, colour = Site, fill = Site),
    size = 3.1,
    stroke = 1,
    show.legend = FALSE
  ) +
  ## dummy legends
  geom_point(
    data = shape_legend_df,
    aes(shape = species_code),
    x = NA, y = NA,
    inherit.aes = FALSE,
    size = 3.1,
    stroke = 1,
    fill = "white",
    colour = "black",
    show.legend = TRUE
  ) +
  geom_point(
    data = fill_legend_df,
    aes(fill = fac_label),
    x = NA, y = NA,
    inherit.aes = FALSE,
    shape = 21,
    size = 3.1,
    stroke = 1,
    colour = "black",
    show.legend = TRUE
  ) +
  facet_wrap(~ Cov, ncol = 1, scales = "free_y", labeller = labeller(Cov = cover_labels)) +
  scale_x_discrete(labels = rep("", length(x_levels))) +
  theme(
    strip.placement = "outside"
  )+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+
  scale_shape_manual(
    values = shape_vals,
    name = "Species",
    labels = c(TRCY = "TRCY", TROR = "TROR")
  ) +
  scale_colour_manual(values = site_cols, drop = FALSE, guide = "none") +
  scale_fill_manual(
    values = c(
      site_cols,
      "Without facilitator" = "white",
      "With facilitator" = "black"
    ),
    breaks = c("Without facilitator", "With facilitator"),
    name = "Facilitator",
    labels = c("Without facilitator" = "No facilitator", "With facilitator" = "Facilitator")
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(colour = "black", fill = "white")),
    fill = guide_legend(order = 2, override.aes = list(shape = 21, colour = "black", size = 3.1))
  ) +
  labs(
    x = NULL,
    y = "Invasion growth rate when rare (log)"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
    strip.text = element_text(face = "plain"),
    legend.position = "right",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    theme_bw(base_size = 13) +
      theme(
        panel.spacing = unit(0.2, "lines")
      )+
    plot.margin = margin(10, 10, 26, 10),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)
  )

p_igr <- p_igr +
  geom_text(
    data = site_centers,
    inherit.aes = FALSE,
    aes(x = x, y = -Inf, label = site_label, colour = Site),
    vjust = 1.7,
    fontface = "bold",
    size = 4,
    show.legend = FALSE
  )

ggsave(
  file.path(out_dir, "figures", paste0(file_prefix, "_plot.png")),
  p_igr,
  width = 8,
  height = 6,
  dpi = 300
)

capture.output(
  {
    cat("\nUsing germination:\n")
    print(use_germination)

    cat("\nSelected environments (fac_sc = 0 vs mean fac_sc among values >= 50th percentile in each Site x Cov):\n")
    print(selected_env %>% arrange(Site, Cov, fac_label))

    cat("\nSite and cover levels used:\n")
    print(levels(dat_stage1$Site))
    print(levels(dat_stage1$Cov))

    cat("\nSummary of posterior support for positive invasion growth:\n")
    print(
      igr_summary %>%
        select(Site, Cov, fac_label, species_code, fac_sc, igr_med, igr_lwr, igr_upr, p_positive) %>%
        arrange(Cov, Site, fac_label, species_code)
    )
  },
  file = file.path(out_dir, "summaries", paste0(file_prefix, "_checks.txt"))
)

print(p_igr)
cat("\nDone. Outputs saved in:\n", normalizePath(out_dir), "\n")
