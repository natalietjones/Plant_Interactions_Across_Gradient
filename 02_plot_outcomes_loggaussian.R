############################################################
## 02_plot_outcomes_loggaussian.R
##
## Purpose:
## 1. Load selected model and prepared data from Script 1
## 2. Extract posterior draws
## 3. Calculate lambda, alpha, niche differences, fitness differences
## 4. Summarise posterior outcome probabilities
## 5. Plot outcomes
############################################################

rm(list = ls())

############################
## 0. Install/load packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "posterior",
  "ggrepel"
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

############################
## 1. Load outputs from Script 1
############################
dat_bh <- readRDS("dat_bh_prepared.rds")
model_selection_info <- readRDS("model_selection_info.rds")
fit_final <- readRDS(model_selection_info$selected_fit_file)
site_order <- model_selection_info$site_order

############################
## 2. Helper: matrix -> long data frame
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
## 3. Build prediction grids
############################
lambda_grid <- expand_grid(
  focsp  = levels(dat_bh$focsp),
  Site   = levels(dat_bh$Site),
  Cov    = levels(dat_bh$Cov),
  fac_sc = sort(unique(dat_bh$fac_sc))
) %>%
  mutate(
    BLK      = levels(dat_bh$BLK)[1],
    intra_sc = 0,
    inter_sc = 0
  )

alpha_grid <- tibble(
  focsp = levels(dat_bh$focsp)
) %>%
  mutate(
    Site     = levels(dat_bh$Site)[1],
    Cov      = levels(dat_bh$Cov)[1],
    fac_sc   = sort(unique(dat_bh$fac_sc))[1],
    BLK      = levels(dat_bh$BLK)[1],
    intra_sc = 0,
    inter_sc = 0
  )

############################
## 4. Extract posterior draws of nonlinear predictors
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
## 5. Convert to lambda and alpha on natural scale
############################
lambda_draws <- mat_to_long(etalam_draws, lambda_grid, "etalam") %>%
  mutate(lambda = exp(etalam)) %>%
  select(draw, focsp, Site, Cov, fac_sc, lambda)

alpha_ii_draws <- mat_to_long(etaai_draws, alpha_grid, "etaai") %>%
  mutate(alpha_ii = exp(etaai)) %>%
  select(draw, focsp, alpha_ii)

alpha_ij_draws <- mat_to_long(etaaj_draws, alpha_grid, "etaaj") %>%
  mutate(alpha_ij = exp(etaaj)) %>%
  select(draw, focsp, alpha_ij)

alpha_draws <- alpha_ii_draws %>%
  left_join(alpha_ij_draws, by = c("draw", "focsp"))

############################
## 6. Save summaries of lambda and alpha
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
  group_by(focsp) %>%
  summarise(
    alpha_ii_med = median(alpha_ii),
    alpha_ii_lwr = quantile(alpha_ii, 0.025),
    alpha_ii_upr = quantile(alpha_ii, 0.975),
    alpha_ij_med = median(alpha_ij),
    alpha_ij_lwr = quantile(alpha_ij, 0.025),
    alpha_ij_upr = quantile(alpha_ij, 0.975),
    .groups = "drop"
  )

write.csv(lambda_summary, "BH_M1_lambda_summary.csv", row.names = FALSE)
write.csv(alpha_summary,  "BH_M1_alpha_summary.csv",  row.names = FALSE)

############################
## 7. Wide format for coexistence calculations
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
  left_join(alpha_wide, by = "draw")

############################
## 8. Germination / seed survival
############################
## Replace if values become available.
## Current version = no-seed-bank approximation
g_vals <- c(TROR = 1, TRCY = 1)
s_vals <- c(TROR = 0, TRCY = 0)

coex_draws <- coex_draws %>%
  mutate(
    eta_TROR = lambda_TROR * g_vals["TROR"] /
      (1 - s_vals["TROR"] * (1 - g_vals["TROR"])),
    eta_TRCY = lambda_TRCY * g_vals["TRCY"] /
      (1 - s_vals["TRCY"] * (1 - g_vals["TRCY"]))
  )

############################
## 9. Niche and fitness differences
############################
## i = TROR
## j = TRCY

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
## 10. Restrict to observed environments only
############################
obs_env <- dat_bh %>%
  distinct(Site, Cov, fac_sc)

coex_draws_obs <- coex_draws %>%
  semi_join(obs_env, by = c("Site", "Cov", "fac_sc"))

env_counts <- dat_bh %>%
  count(Site, Cov, fac_sc, name = "n_obs")

############################
## 11. Summarise propagated uncertainty by environment
############################
coex_summary <- coex_draws_obs %>%
  group_by(Site, Cov, fac_sc) %>%
  summarise(
    niche_med = median(niche_diff),
    niche_lwr = quantile(niche_diff, 0.025),
    niche_upr = quantile(niche_diff, 0.975),
    
    fit_med = median(fitness_ratio_TRCY_over_TROR),
    fit_lwr = quantile(fitness_ratio_TRCY_over_TROR, 0.025),
    fit_upr = quantile(fitness_ratio_TRCY_over_TROR, 0.975),
    
    p_coexist = mean(coexist),
    p_TRCY_wins = mean(TRCY_wins),
    p_TROR_wins = mean(TROR_wins),
    p_priority = mean(priority_effect),
    
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
    Site = factor(Site, levels = site_order)
  )

coex_draws_obs <- coex_draws_obs %>%
  mutate(
    Site = factor(Site, levels = site_order)
  )

write.csv(coex_summary, "BH_M1_coexistence_summary.csv", row.names = FALSE)
saveRDS(lambda_draws,   "BH_M1_lambda_draws.rds")
saveRDS(alpha_draws,    "BH_M1_alpha_draws.rds")
saveRDS(coex_draws_obs, "BH_M1_coexistence_draws_observed_env.rds")

############################
## 12. Quick checks
############################
print(coex_summary %>% count(outcome_90))
print(coex_summary %>% count(most_supported_outcome))
print(summary(coex_summary$max_support))

############################
## 13. Posterior support figure
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
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(Cov ~ Site) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Facilitator density",
    y = "Posterior probability of outcome",
    colour = "Outcome"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

print(p_support)

ggsave(
  "BH_M1_coexistence_support_by_environment.png",
  p_support,
  width = 12,
  height = 7,
  dpi = 300
)

############################
## 14. Coexistence plane summary plot
############################
plot_env <- coex_summary %>%
  mutate(
    Cov = factor(Cov),
    fac_sc = as.numeric(fac_sc)
  ) %>%
  filter(
    is.finite(niche_med),
    is.finite(fit_med),
    fit_med > 0,
    is.finite(max_support),
    niche_med >= 0,
    niche_med <= 0.25
  ) %>%
  arrange(Site, Cov, fac_sc)

boundary_dat <- tibble(
  niche_diff = seq(0, 0.25, length.out = 1000)
) %>%
  mutate(
    lower = 1 - niche_diff,
    upper = 1 / (1 - niche_diff)
  ) %>%
  filter(lower > 0, upper > 0)

region_labels <- tibble(
  niche_diff = c(0.08, 0.18, 0.18, 0.12),
  fit = c(1.05, 1.8, 0.7, 1.0),
  lab = c("Priority\neffect", "TRCY wins", "TROR wins", "Coexistence")
)

p_plane <- ggplot() +
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
  geom_path(
    data = plot_env,
    aes(
      x = niche_med,
      y = fit_med,
      group = interaction(Site, Cov),
      linetype = Cov
    ),
    linewidth = 0.7,
    alpha = 0.6,
    colour = "grey30"
  ) +
  geom_point(
    data = plot_env,
    aes(
      x = niche_med,
      y = fit_med,
      colour = most_supported_outcome,
      shape = Cov,
      size = max_support
    ),
    alpha = 0.9
  ) +
  geom_text_repel(
    data = plot_env,
    aes(
      x = niche_med,
      y = fit_med,
      label = fac_sc
    ),
    size = 2.8,
    seed = 123,
    box.padding = 0.2,
    point.padding = 0.1,
    min.segment.length = 0,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  geom_text(
    data = region_labels,
    aes(x = niche_diff, y = fit, label = lab),
    size = 3.5
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = c(
      "coexist" = "#F4A261",
      "TRCY_wins" = "#2A9D8F",
      "TROR_wins" = "#457B9D",
      "priority_effect" = "grey50"
    )
  ) +
  scale_size_continuous(range = c(2, 5), limits = c(0, 1)) +
  labs(
    x = "Niche difference (1 - rho)",
    y = "Fitness difference (TRCY / TROR)",
    colour = "Most supported outcome",
    shape = "Cov",
    linetype = "Cov",
    size = "Max support"
  ) +
  coord_cartesian(xlim = c(0, 0.25)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_plane)

ggsave(
  "BH_M1_coexistence_plane_summary.png",
  p_plane,
  width = 10,
  height = 7,
  dpi = 300
)