############################################################
## 02_BH_Bayesian_plot_outcomes.R
##
## Main checks:
## 1. Coexistence formulas for rho and kappa ratio checked against
##    Godoy & Levine (2014) / Bowler et al. (2022) annual-plant BH equations.
## 2. LDGR plotting section corrected so invasion growth when rare does NOT
##    include invader self-competition (alpha_ii term).
## 
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

############################
## 1. User options
############################
use_germination <- FALSE

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

############################
## 3b. Define facilitator levels for none vs high plots
############################
selected_fac_levels <- bind_rows(
  dat_stage1 %>%
    distinct(Site, Cov) %>%
    mutate(
      fac_sc = 0,
      fac_level = "None"
    ),
  dat_stage1 %>%
    group_by(Site, Cov) %>%
    summarise(
      fac_sc = {
        x <- fac_sc[is.finite(fac_sc) & fac_sc > 0]
        if (length(x) == 0) {
          0
        } else {
          q <- quantile(x, 0.90, na.rm = TRUE)
          mean(x[x >= q], na.rm = TRUE)
        }
      },
      n_positive = sum(fac_sc > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      fac_level = "High"
    )
) %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN")),
    fac_level = factor(fac_level, levels = c("None", "High"))
  )

selected_fac_levels %>%
  pivot_wider(names_from = fac_level, values_from = fac_sc) %>%
  mutate(delta = High - None) %>%
  arrange(Site, Cov)

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
fac_values_all <- sort(unique(c(
  dat_stage1$fac_sc,
  selected_fac_levels$fac_sc
)))

lambda_grid <- expand_grid(
  focsp  = levels(dat_stage1$focsp),
  Site   = levels(dat_stage1$Site),
  Cov    = levels(dat_stage1$Cov),
  fac_sc = fac_values_all
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

## alpha_ij_focspTRCY = effect of TROR on TRCY
## alpha_ij_focspTROR = effect of TRCY on TROR
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
## Godoy & Levine / Bowler annual-plant equations:
## rho = sqrt( (alpha_ij * alpha_ji) / (alpha_ii * alpha_jj) )
## kappa_TRCY / kappa_TROR = ((eta_TRCY - 1)/(eta_TROR - 1)) *
##   sqrt( (alpha_TROR,TRCY / alpha_TRCY,TRCY) * (alpha_TROR,TROR / alpha_TRCY,TROR) )
coex_draws <- coex_draws %>%
  mutate(
    rho = sqrt(
      (alpha_ij_TROR * alpha_ij_TRCY) /
        (alpha_ii_TROR * alpha_ii_TRCY)
    ),
    niche_diff = 1 - rho,
    
    fitness_ratio_TRCY_over_TROR =
      ((eta_TRCY - 1) / (eta_TROR - 1)) *
      sqrt(
        (alpha_ij_TROR * alpha_ii_TROR) /
          (alpha_ii_TRCY * alpha_ij_TRCY)
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
    Site = factor(as.character(Site), levels = c("BEN", "GH", "NAM", "PJ", "CS")),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN"))
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
## 15. Quick checks
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
  },
  file = file.path(outcomes_dir, "summaries", paste0(file_prefix, "_outcome_checks_checked.txt"))
)

############################
## 16. LDGR distributions across sites
## Matched to q50 none-vs-high facilitator definition
############################

site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

selected_env_ldgr <- selected_fac_levels %>%
  mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = c("SH", "SUN")),
    panel_row = paste(Cov, fac_level, sep = "-"),
    panel_row = factor(panel_row, levels = c("SH-None", "SH-High", "SUN-None", "SUN-High"))
  )

## IMPORTANT:
## Use coex_draws (not coex_draws_obs), because the q50 "High"
## fac_sc value is often not an observed value.
ldgr_base <- coex_draws %>%
  inner_join(
    selected_env_ldgr %>%
      select(Site, Cov, fac_sc, fac_level, panel_row),
    by = c("Site", "Cov", "fac_sc")
  )

ldgr_draws <- ldgr_base %>%
  mutate(
    ## Resident monoculture equilibria
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
    
    ## Invasion growth when rare: omit invader self-competition
    ldgr_TROR = log(
      (1 - g_TROR) * s_vals["TROR"] +
        (g_TROR * lambda_TROR) /
        (1 + alpha_ij_TROR * g_TRCY * Nstar_TRCY)
    ),
    
    ldgr_TRCY = log(
      (1 - g_TRCY) * s_vals["TRCY"] +
        (g_TRCY * lambda_TRCY) /
        (1 + alpha_ij_TRCY * g_TROR * Nstar_TROR)
    )
  )

site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

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
    ),
    Site = factor(as.character(Site), levels = site_order),
    panel_row = factor(panel_row, levels = c("SH-None", "SH-High", "SUN-None", "SUN-High"))
  ) %>%
  filter(is.finite(ldgr))
ldgr_points <- ldgr_long %>%
  group_by(Site, panel_row, species) %>%
  summarise(
    ldgr_point = median(ldgr),
    .groups = "drop"
  ) %>%
  mutate(Site = factor(as.character(Site), levels = site_order))

p_coexist_ldgr <- ldgr_draws %>%
  group_by(Site, panel_row) %>%
  summarise(
    p_coexist = mean(ldgr_TROR > 0 & ldgr_TRCY > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prob_lab = paste0(round(100 * p_coexist), "%"),
    Site = factor(as.character(Site), levels = site_order)
  )

ldgr_xpos <- ldgr_long %>%
  group_by(Site, panel_row) %>%
  summarise(
    x_lab = quantile(ldgr, 0.005, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Site = factor(as.character(Site), levels = site_order))

prob_plot <- p_coexist_ldgr %>%
  left_join(ldgr_xpos, by = c("Site", "panel_row")) %>%
  mutate(Site = factor(as.character(Site), levels = site_order))

species_cols <- c(
  "TRCY" = "#065143",
  "TROR" = "#F2C14E"
)

p_ldgr <- ggplot(ldgr_long, aes(x = ldgr, colour = species)) +
  geom_density(linewidth = 1, adjust = 1.2, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(
    data = ldgr_points,
    aes(x = ldgr_point, y = 0, colour = species),
    inherit.aes = FALSE,
    size = 2,
    alpha = 0.6
  ) +
  geom_text(
    data = prob_plot,
    aes(x = x_lab, y = Inf, label = prob_lab),
    inherit.aes = FALSE,
    hjust = 1.3,
    vjust = 1.4,
    size = 4,
    colour = "black"
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
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "top",
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x  = element_text(size = 11, colour = "black"),
    axis.text.y  = element_text(size = 11, colour = "black"),
    strip.text   = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    legend.text  = element_text(size = 11, colour = "black")
  ) +
  guides(alpha = "none")

ggsave(
  file.path(outcomes_dir, "figures", paste0(file_prefix, "_LDGR_distributions_sites_q90.png")),
  p_ldgr,
  width = 8,
  height = 7,
  dpi = 300
)

print(p_ldgr)

