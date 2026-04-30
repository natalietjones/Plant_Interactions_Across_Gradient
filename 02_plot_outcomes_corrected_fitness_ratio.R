############################################################
## 02_plot_outcomes.R
##
## Purpose:
## 1. Load the selected structure model from 01_alpha_structure_selection/
## 2. Extract posterior draws for lambda and alpha
## 3. Optionally load site-specific germination priors and draw g
## 4. Propagate uncertainty into niche and fitness differences
## 5. Save summary tables for outcomes
## 6. Make outcome figures
##
## Main setting:
##   use_germination <- FALSE
##
## This makes the primary analysis use the competition experiment
## alone (post-germination fitness), while allowing a germination-
## weighted sensitivity analysis by setting use_germination <- TRUE.
##
## Note:
## - no seed bank in either version (s = 0)
## - this script is currently written for the BH branch only
############################################################

rm(list = ls())

############################
## 0. Install/load packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "posterior",
  "ggrepel",
  "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(brms)
library(posterior)
library(ggrepel)
library(patchwork)

options(mc.cores = max(1, parallel::detectCores() - 1))
options(brms.backend = "rstan")
# options(brms.backend = "cmdstanr")

############################
## 1. User options
############################
## Primary analysis = FALSE
## Sensitivity analysis = TRUE
use_germination <- TRUE

############################
## 2. Set input/output folders
############################
structure_dir <- "01_alpha_structure_selection"
outcomes_dir  <- "02_outcomes"

dir.create(outcomes_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outcomes_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outcomes_dir, "draws"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outcomes_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outcomes_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)

############################
## 3. Load selected model and prepared data
############################
alpha_info <- readRDS(
  file.path(structure_dir, "selection", "alpha_structure_selection_info.rds")
)

if (!identical(alpha_info$selected_growth_family, "BH")) {
  stop(
    paste0(
      "This script currently implements BH coexistence calculations only.\n",
      "Selected family is: ", alpha_info$selected_growth_family
    )
  )
}

dat_stage1 <- readRDS(alpha_info$stage1_data_rds)
fit_final  <- readRDS(alpha_info$selected_fit_file)
site_order <- alpha_info$site_order

family_tag <- alpha_info$selected_growth_family
error_tag  <- ifelse(alpha_info$selected_error_model == "gaussian_logfit",
                     "loggaussian",
                     alpha_info$selected_error_model)

prefix_tag <- ifelse(use_germination, "germ", "nogerm")
file_prefix <- paste0(family_tag, "_", error_tag, "_", prefix_tag)

cat("\nLoaded structure-selection info:\n")
print(alpha_info)

cat("\nUsing germination in outcome calculations:\n")
print(use_germination)

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
## 5. Build prediction grids
############################
lambda_grid <- expand_grid(
  focsp  = levels(dat_stage1$focsp),
  Site   = levels(dat_stage1$Site),
  Cov    = levels(dat_stage1$Cov),
  fac_sc = sort(unique(dat_stage1$fac_sc))
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
    fac_sc   = sort(unique(dat_stage1$fac_sc))[1],
    BLK      = levels(dat_stage1$BLK)[1],
    intra_sc = 0,
    inter_sc = 0
  )

############################
## 6. Extract posterior draws of nonlinear predictors
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
## 7. Convert to lambda and alpha on natural scale
############################
lambda_draws <- mat_to_long(etalam_draws, lambda_grid, "etalam") %>%
  mutate(lambda = exp(etalam)) %>%
  select(draw, focsp, Site, Cov, fac_sc, lambda)

alpha_ii_draws <- mat_to_long(etaai_draws, alpha_grid, "etaai") %>%
  mutate(alpha_ii = exp(etaai)) %>%
  select(draw, focsp, Site, Cov, alpha_ii)

alpha_ij_draws <- mat_to_long(etaaj_draws, alpha_grid, "etaaj") %>%
  mutate(alpha_ij = exp(etaaj)) %>%
  select(draw, focsp, Site, Cov, alpha_ij)

alpha_draws <- alpha_ii_draws %>%
  left_join(alpha_ij_draws, by = c("draw", "focsp", "Site", "Cov"))

############################
## 8. Save summaries of lambda and alpha
############################
lambda_summary <- lambda_draws %>%
  group_by(focsp, Site, Cov, fac_sc) %>%
  summarise(
    lambda_med = median(lambda),
    lambda_lwr = quantile(lambda, 0.025),
    lambda_upr = quantile(lambda, 0.975),
    .groups = "drop"
  )

alpha_summary <- alpha_draws %>%
  group_by(focsp, Site, Cov) %>%
  summarise(
    alpha_ii_med = median(alpha_ii),
    alpha_ii_lwr = quantile(alpha_ii, 0.025),
    alpha_ii_upr = quantile(alpha_ii, 0.975),
    alpha_ij_med = median(alpha_ij),
    alpha_ij_lwr = quantile(alpha_ij, 0.025),
    alpha_ij_upr = quantile(alpha_ij, 0.975),
    .groups = "drop"
  )

write.csv(
  lambda_summary,
  file.path(outcomes_dir, "tables", paste0(file_prefix, "_lambda_summary.csv")),
  row.names = FALSE
)

write.csv(
  alpha_summary,
  file.path(outcomes_dir, "tables", paste0(file_prefix, "_alpha_summary.csv")),
  row.names = FALSE
)

############################
## 9. Germination draws OR no-germination default
############################
n_draws <- n_distinct(lambda_draws$draw)

if (use_germination) {
  
  if (file.exists(alpha_info$germination_site_priors_rds)) {
    germ_priors <- readRDS(alpha_info$germination_site_priors_rds) %>%
      as_tibble()
  } else {
    germ_info <- readRDS(alpha_info$germination_priors_rds)
    germ_priors <- germ_info$site_priors %>%
      as_tibble()
  }
  
  germ_priors <- germ_priors %>%
    mutate(
      Site  = as.character(Site),
      focsp = as.character(focsp)
    )
  
  required_g_cols <- c("Site", "focsp", "g_alpha_final", "g_beta_final", "g_mean_final")
  missing_g_cols <- setdiff(required_g_cols, names(germ_priors))
  
  if (length(missing_g_cols) > 0) {
    stop(
      "Missing required germination columns: ",
      paste(missing_g_cols, collapse = ", ")
    )
  }
  
  set.seed(123)
  
  germ_draws <- expand_grid(
    draw = seq_len(n_draws),
    germ_priors %>% select(Site, focsp, g_alpha_final, g_beta_final, g_mean_final)
  ) %>%
    mutate(
      g_draw = rbeta(
        n = n(),
        shape1 = g_alpha_final,
        shape2 = g_beta_final
      )
    )
  
  germ_summary <- germ_draws %>%
    group_by(Site, focsp) %>%
    summarise(
      g_input_mean = first(g_mean_final),
      g_draw_med   = median(g_draw),
      g_draw_lwr   = quantile(g_draw, 0.025),
      g_draw_upr   = quantile(g_draw, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Site = factor(Site, levels = site_order)) %>%
    arrange(focsp, Site)
  
  germ_wide <- germ_draws %>%
    select(draw, Site, focsp, g_draw) %>%
    pivot_wider(
      names_from  = focsp,
      values_from = g_draw,
      names_glue  = "g_{focsp}"
    )
  
} else {
  
  germ_summary <- expand_grid(
    Site = levels(dat_stage1$Site),
    focsp = levels(dat_stage1$focsp)
  ) %>%
    mutate(
      g_input_mean = 1,
      g_draw_med   = 1,
      g_draw_lwr   = 1,
      g_draw_upr   = 1
    ) %>%
    mutate(Site = factor(Site, levels = site_order)) %>%
    arrange(focsp, Site)
  
  germ_wide <- expand_grid(
    draw = seq_len(n_draws),
    Site = levels(dat_stage1$Site)
  ) %>%
    mutate(
      g_TRCY = 1,
      g_TROR = 1
    )
}

cat("\nGermination summary used in outcome calculations:\n")
print(germ_summary)

write.csv(
  germ_summary,
  file.path(outcomes_dir, "tables", paste0(file_prefix, "_germination_summary.csv")),
  row.names = FALSE
)

saveRDS(
  germ_wide,
  file.path(outcomes_dir, "draws", paste0(file_prefix, "_germination_draws_or_constants.rds"))
)

############################
## 10. Wide format for coexistence calculations
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

coex_draws <- lambda_wide %>%
  left_join(alpha_wide, by = c("draw", "Site", "Cov")) %>%
  left_join(germ_wide,  by = c("draw", "Site"))

############################
## 11. Germination / seed survival
############################
## No seed bank in either version
s_vals <- c(TROR = 0, TRCY = 0)

coex_draws <- coex_draws %>%
  mutate(
    eta_TROR = lambda_TROR * g_TROR /
      (1 - s_vals["TROR"] * (1 - g_TROR)),
    eta_TRCY = lambda_TRCY * g_TRCY /
      (1 - s_vals["TRCY"] * (1 - g_TRCY))
  )

############################
## 12. Niche and fitness differences
############################
coex_draws <- coex_draws %>%
  mutate(
    rho = sqrt(
      (alpha_ij_TROR * alpha_ij_TRCY) /
        (alpha_ii_TROR * alpha_ii_TRCY)
    ),
    niche_diff = 1 - rho,
    
    ## Following Godoy & Levine (2014) / Van Dyke et al. (2022):
    ## kappa_i / kappa_j = ((eta_i - 1)/(eta_j - 1)) * sqrt((alpha_ij/alpha_jj) * (alpha_ii/alpha_ji))
    ## Here i = TRCY and j = TROR.
    fitness_ratio_TRCY_over_TROR =
      ((eta_TRCY - 1) / (eta_TROR - 1)) *
      sqrt(
        (alpha_ij_TRCY / alpha_ii_TROR) *
          (alpha_ii_TRCY / alpha_ij_TROR)
      ),
    
    inv_rho = 1 / rho,
    
    coexist = rho < fitness_ratio_TRCY_over_TROR &
      fitness_ratio_TRCY_over_TROR < inv_rho,
    
    TRCY_wins = fitness_ratio_TRCY_over_TROR > inv_rho,
    TROR_wins = fitness_ratio_TRCY_over_TROR < rho,
    priority_effect = !coexist & !TRCY_wins & !TROR_wins
  )

############################
## 13. Restrict to observed environments only
############################
obs_env <- dat_stage1 %>%
  distinct(Site, Cov, fac_sc)

coex_draws_obs <- coex_draws %>%
  semi_join(obs_env, by = c("Site", "Cov", "fac_sc"))

env_counts <- dat_stage1 %>%
  count(Site, Cov, fac_sc, name = "n_obs")

############################
## 14. Summarise propagated uncertainty by environment
############################
coex_summary <- coex_draws_obs %>%
  group_by(Site, Cov, fac_sc) %>%
  summarise(
    niche_med = median(niche_diff, na.rm = TRUE),
    niche_lwr = quantile(niche_diff, 0.025, na.rm = TRUE),
    niche_upr = quantile(niche_diff, 0.975, na.rm = TRUE),
    
    fit_med = median(fitness_ratio_TRCY_over_TROR, na.rm = TRUE),
    fit_lwr = quantile(fitness_ratio_TRCY_over_TROR, 0.025, na.rm = TRUE),
    fit_upr = quantile(fitness_ratio_TRCY_over_TROR, 0.975, na.rm = TRUE),
    
    p_coexist = mean(coexist, na.rm = TRUE),
    p_TRCY_wins = mean(TRCY_wins, na.rm = TRUE),
    p_TROR_wins = mean(TROR_wins, na.rm = TRUE),
    p_priority = mean(priority_effect, na.rm = TRUE),
    
    mean_g_TROR = mean(g_TROR, na.rm = TRUE),
    mean_g_TRCY = mean(g_TRCY, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  left_join(env_counts, by = c("Site", "Cov", "fac_sc")) %>%
  rowwise() %>%
  mutate(
    max_support = max(c(p_coexist, p_TRCY_wins, p_TROR_wins, p_priority)),
    most_supported_outcome = c("coexist", "TRCY_wins", "TROR_wins", "priority_effect")[
      which.max(c(p_coexist, p_TRCY_wins, p_TROR_wins, p_priority))
    ],
    outcome_90 = ifelse(max_support >= 0.90, most_supported_outcome, "ambiguous")
  ) %>%
  ungroup() %>%
  mutate(
    Site = factor(Site, levels = c("BEN", "PJ", "NAM", "GH", "CS")),
    Cov  = factor(Cov, levels = c("SH", "SUN"))
  )

write.csv(
  coex_summary,
  file.path(outcomes_dir, "tables", paste0(file_prefix, "_coexistence_summary.csv")),
  row.names = FALSE
)

saveRDS(
  lambda_draws,
  file.path(outcomes_dir, "draws", paste0(file_prefix, "_lambda_draws.rds"))
)
saveRDS(
  alpha_draws,
  file.path(outcomes_dir, "draws", paste0(file_prefix, "_alpha_draws.rds"))
)
saveRDS(
  coex_draws_obs,
  file.path(outcomes_dir, "draws", paste0(file_prefix, "_coexistence_draws_observed_env.rds"))
)

############################
## 15. Save quick checks and notes
############################
capture.output(
  {
    cat("\nUsing germination:\n")
    print(use_germination)
    
    cat("\nOutcome-90 counts:\n")
    print(coex_summary %>% count(outcome_90))
    
    cat("\nMost-supported-outcome counts:\n")
    print(coex_summary %>% count(most_supported_outcome))
    
    cat("\nSummary of max support:\n")
    print(summary(coex_summary$max_support))
    
    cat("\nTop environments by max support:\n")
    print(
      coex_summary %>%
        arrange(desc(max_support)) %>%
        select(Site, Cov, fac_sc, most_supported_outcome, max_support,
               p_coexist, p_TRCY_wins, p_TROR_wins, p_priority,
               mean_g_TROR, mean_g_TRCY, n_obs)
    )
  },
  file = file.path(outcomes_dir, "summaries", paste0(file_prefix, "_outcome_checks.txt"))
)


site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

coex_summary <- coex_summary %>%
  mutate(
    Site = factor(as.character(Site), levels = c("BEN", "GH", "NAM", "PJ", "CS")),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN"))
  )

coex_draws_obs <- coex_draws_obs %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN"))
  )


############################
## 16. Posterior support figure
############################
plot_dat <- coex_summary %>%
  select(Site, Cov, fac_sc, p_coexist, p_TRCY_wins, p_TROR_wins, p_priority) %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "outcome",
    values_to = "prob"
  ) %>%
  mutate(
    outcome = recode(
      outcome,
      p_coexist   = "Coexist",
      p_TRCY_wins = "TRCY wins",
      p_TROR_wins = "TROR wins",
      p_priority  = "Priority effect"
    ),
    outcome = factor(
      outcome,
      levels = c("Coexist", "TRCY wins", "TROR wins", "Priority effect")
    )
  )

p_support <- ggplot(plot_dat, aes(x = fac_sc, y = prob, colour = outcome)) +
  geom_line(linewidth = 0.8, alpha = 0.6) +
  geom_point(size = 1.8) +
  facet_grid(Cov ~ Site) +
  scale_colour_manual(
    values = c(
      "Coexist" = "#09BC8A",
      "TRCY wins" = "#FA7921",
      "TROR wins" = "#446DF6",
      "Priority effect" = "#9F9F92"
    )
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Facilitator density",
    y = "Posterior probability of outcome",
    colour = "Outcome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text = element_text(colour = "black"),
    legend.position = "bottom"
  )

ggsave(
  file.path(outcomes_dir, "figures", paste0(file_prefix, "_coexistence_support_by_environment.png")),
  p_support,
  width = 12,
  height = 7,
  dpi = 300
)

############################
## 17. Coexistence plane with posterior draws
############################
selected_env <- bind_rows(
  coex_summary %>%
    group_by(Site, Cov) %>%
    slice_min(fac_sc, n = 1, with_ties = FALSE) %>%
    mutate(fac_level = "Low"),
  coex_summary %>%
    group_by(Site, Cov) %>%
    slice_max(fac_sc, n = 1, with_ties = FALSE) %>%
    mutate(fac_level = "High")
) %>%
  ungroup() %>%
  distinct(Site, Cov, fac_sc, .keep_all = TRUE) %>%
  mutate(
    condition = paste(Cov, fac_level, sep = "-")
  )

draw_plot <- coex_draws_obs %>%
  inner_join(
    selected_env %>%
      select(Site, Cov, fac_sc, fac_level, condition),
    by = c("Site", "Cov", "fac_sc")
  ) %>%
  filter(
    is.finite(niche_diff),
    is.finite(fitness_ratio_TRCY_over_TROR),
    fitness_ratio_TRCY_over_TROR > 0,
    niche_diff >= 0,
    niche_diff <= 0.25
  )

set.seed(123)
draw_plot_thin <- draw_plot %>%
  split(list(.$Site, .$condition), drop = TRUE) %>%
  purrr::map_dfr(~ dplyr::slice_sample(.x, n = min(1500, nrow(.x))))

point_plot <- selected_env %>%
  filter(
    is.finite(niche_med),
    is.finite(fit_med),
    fit_med > 0,
    niche_med >= 0,
    niche_med <= 0.25
  )

boundary_dat <- tibble(
  niche_diff = seq(0, 0.25, length.out = 1000)
) %>%
  mutate(
    lower = 1 - niche_diff,
    upper = 1 / (1 - niche_diff)
  ) %>%
  filter(lower > 0, upper > 0)

region_labels <- tibble(
  niche_diff = c(0.06, 0.18, 0.18, 0.12),
  fit = c(1.05, 1.8, 0.7, 1.0),
  lab = c("Priority\neffect", "TRCY wins", "TROR wins", "Coexistence")
)

cond_cols <- c(
  "SH-Low"   = "#1F78B4",
  "SH-High"  = "#08306B",
  "SUN-Low"  = "#E6AB02",
  "SUN-High" = "#D95F02"
)

cond_shapes <- c(
  "SH-Low"   = 16,
  "SH-High"  = 17,
  "SUN-Low"  = 15,
  "SUN-High" = 18
)

p_plane_draws <- ggplot() +
  geom_ribbon(
    data = boundary_dat,
    aes(x = niche_diff, ymin = lower, ymax = upper),
    fill = "grey90"
  ) +
  geom_line(
    data = boundary_dat,
    aes(x = niche_diff, y = lower),
    linewidth = 0.8
  ) +
  geom_line(
    data = boundary_dat,
    aes(x = niche_diff, y = upper),
    linewidth = 0.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(
    data = draw_plot_thin,
    aes(
      x = niche_diff,
      y = fitness_ratio_TRCY_over_TROR,
      colour = condition
    ),
    alpha = 0.08,
    size = 0.5
  ) +
  geom_point(
    data = point_plot,
    aes(
      x = niche_med,
      y = fit_med,
      colour = condition,
      shape = condition
    ),
    size = 2.8,
    stroke = 0.9
  ) +
  geom_text_repel(
    data = point_plot,
    aes(
      x = niche_med,
      y = fit_med,
      label = condition,
      colour = condition
    ),
    size = 2.8,
    seed = 123,
    box.padding = 0.2,
    point.padding = 0.12,
    min.segment.length = 0,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  geom_text(
    data = region_labels,
    aes(x = niche_diff, y = fit, label = lab),
    size = 3.1
  ) +
  facet_wrap(~ Site, nrow = 1) +
  scale_y_log10() +
  scale_colour_manual(values = cond_cols, drop = FALSE) +
  scale_shape_manual(values = cond_shapes, drop = FALSE) +
  labs(
    x = "Niche difference (1 - rho)",
    y = "Fitness difference (TRCY / TROR)",
    colour = "Condition",
    shape = "Condition"
  ) +
  coord_cartesian(xlim = c(0, 0.25)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey95")
  )

ggsave(
  file.path(outcomes_dir, "figures", paste0(file_prefix, "_coexistence_plane_draws_by_condition.png")),
  p_plane_draws,
  width = 14,
  height = 4.8,
  dpi = 300
)

############################
## 18. LDGR distributions across sites
############################
site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

coex_summary <- coex_summary %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN")),
    fac_sc = as.numeric(fac_sc)
  )

coex_draws_obs <- coex_draws_obs %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN")),
    fac_sc = as.numeric(fac_sc)
  )

selected_env_ldgr <- bind_rows(
  coex_summary %>%
    group_by(Site, Cov) %>%
    slice_min(fac_sc, n = 1, with_ties = FALSE) %>%
    mutate(fac_level = "Low"),
  coex_summary %>%
    group_by(Site, Cov) %>%
    slice_max(fac_sc, n = 1, with_ties = FALSE) %>%
    mutate(fac_level = "High")
) %>%
  ungroup() %>%
  distinct(Site, Cov, fac_sc, .keep_all = TRUE) %>%
  mutate(
    panel_row = paste(Cov, fac_level, sep = "-"),
    panel_row = factor(panel_row, levels = c("SH-Low", "SH-High", "SUN-Low", "SUN-High"))
  )

ldgr_base <- coex_draws_obs %>%
  inner_join(
    selected_env_ldgr %>%
      select(Site, Cov, fac_sc, fac_level, panel_row),
    by = c("Site", "Cov", "fac_sc")
  )

ldgr_draws <- ldgr_base %>%
  mutate(
    Nstar_TROR = ifelse(
      eta_TROR > 1,
      (eta_TROR - 1) / (alpha_ii_TROR * g_TROR),
      0
    ),
    Nstar_TRCY = ifelse(
      eta_TRCY > 1,
      (eta_TRCY - 1) / (alpha_ii_TRCY * g_TRCY),
      0
    ),
    
    ldgr_TROR = log(
      (1 - g_TROR) * s_vals["TROR"] +
        (g_TROR * lambda_TROR) /
        (1 + alpha_ii_TROR * g_TROR * 1 + alpha_ij_TROR * g_TRCY * Nstar_TRCY)
    ),
    
    ldgr_TRCY = log(
      (1 - g_TRCY) * s_vals["TRCY"] +
        (g_TRCY * lambda_TRCY) /
        (1 + alpha_ii_TRCY * g_TRCY * 1 + alpha_ij_TRCY * g_TROR * Nstar_TROR)
    )
  )

ldgr_long <- ldgr_draws %>%
  select(draw, Site, Cov, fac_sc, fac_level, panel_row, ldgr_TROR, ldgr_TRCY) %>%
  pivot_longer(
    cols = c(ldgr_TROR, ldgr_TRCY),
    names_to = "species",
    values_to = "ldgr"
  ) %>%
  mutate(
    species = recode(
      species,
      ldgr_TROR = "TROR",
      ldgr_TRCY = "TRCY"
    )
  ) %>%
  filter(is.finite(ldgr))

ldgr_points <- ldgr_long %>%
  group_by(Site, panel_row, species) %>%
  summarise(
    ldgr_point = median(ldgr),
    .groups = "drop"
  )

p_coexist_ldgr <- ldgr_draws %>%
  group_by(Site, panel_row) %>%
  summarise(
    p_coexist = mean(ldgr_TROR > 0 & ldgr_TRCY > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prob_lab = paste0(sprintf("%.1f", 100 * p_coexist), "%")
  )

ldgr_xpos <- ldgr_long %>%
  group_by(Site, panel_row) %>%
  summarise(
    x_lab = quantile(ldgr, 0.98, na.rm = TRUE),
    .groups = "drop"
  )

prob_plot <- p_coexist_ldgr %>%
  left_join(ldgr_xpos, by = c("Site", "panel_row"))

species_cols <- c(
  "TROR" = "#004346",
  "TRCY" = "#F2C14E"
)


p_ldgr <- ggplot(ldgr_long, aes(x = ldgr, colour = species, alpha = 0.6)) +
  geom_density(linewidth = 0.8, adjust = 1.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(
    data = ldgr_points,
    aes(x = ldgr_point, y = 0, colour = species),
    inherit.aes = FALSE,
    size = 1.8
  ) +
  geom_text(
    data = prob_plot,
    aes(x = x_lab, y = Inf, label = prob_lab),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.2,
    size = 3.1
  ) +
  facet_grid(panel_row ~ Site, scales = "free_y") +
  scale_colour_manual(values = species_cols, drop = FALSE) +
  labs(
    x = "Low-density growth rate (LDGR)",
    y = "Density",
    colour = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "top"
  )

ggsave(
  file.path(outcomes_dir, "figures", paste0(file_prefix, "_LDGR_distributions_sites_lowhighfac.png")),
  p_ldgr,
  width = 12,
  height = 8,
  dpi = 300
)

p_ldgr

############################
## 19. Finish
############################
cat("\nDone.\n")
cat("Using germination:\n", use_germination, "\n")
cat("Selected family:\n", alpha_info$selected_growth_family, "\n")
cat("Selected structure model:\n", alpha_info$selected_structure_model, "\n")
cat("Outputs saved in:\n", normalizePath(outcomes_dir), "\n")
