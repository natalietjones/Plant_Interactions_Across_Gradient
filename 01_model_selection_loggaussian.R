############################################################
## 01_model_selection_loggaussian.R
##
## Purpose:
## 1. Prepare data
## 2. Fit/load candidate log-Gaussian Beverton-Holt models
## 3. Compare models with LOO
## 4. Save model-selection outputs for downstream plotting
############################################################

rm(list = ls())

############################
## 0. Install/load packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "loo",
  "posterior"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(brms)
library(loo)
library(posterior)

options(mc.cores = max(1, parallel::detectCores() - 1))
options(brms.backend = "rstan")
# options(brms.backend = "cmdstanr")

############################
## 1. Read and prepare data
############################
dat0 <- read.csv("fullWAdata.csv")

dat_bh <- dat0 %>%
  mutate(
    Site     = factor(Site),
    Cov      = factor(Cov),
    BLK      = factor(BLK),
    focsp    = factor(focsp),
    dens_trt = factor(dens_trt),
    germ     = FOCAL_GERM %in% c("Y", "y", 1, "1", TRUE)
  ) %>%
  filter(germ) %>%
  mutate(
    intra = case_when(
      focsp == "TROR" ~ TROR,
      focsp == "TRCY" ~ TRCY,
      TRUE ~ NA_real_
    ),
    inter = case_when(
      focsp == "TROR" ~ TRCY,
      focsp == "TRCY" ~ TROR,
      TRUE ~ NA_real_
    ),
    fac = NO_GRASS
  ) %>%
  filter(
    !is.na(fitness),
    !is.na(intra),
    !is.na(inter),
    !is.na(fac),
    !is.na(Site),
    !is.na(Cov),
    !is.na(BLK),
    !is.na(focsp)
  ) %>%
  mutate(
    intra_sc    = intra,
    inter_sc    = inter,
    fac_sc      = fac,
    fitness_raw = fitness,
    fitness_log = log(fitness + 1)
  ) %>%
  droplevels()

saveRDS(dat_bh, "dat_bh_prepared.rds")

############################
## 2. Quick checks
############################
cat("\nSummary of raw fitness:\n")
print(summary(dat_bh$fitness_raw))

cat("\nSummary of log(fitness + 1):\n")
print(summary(dat_bh$fitness_log))

cat("\nLevels:\n")
print(levels(dat_bh$focsp))
print(levels(dat_bh$Site))
print(levels(dat_bh$Cov))

############################
## 3. Family, priors, settings
############################
bh_family <- gaussian()

priors_common <- c(
  prior(normal(0, 1.5), nlpar = "etalam", class = "b"),
  prior(normal(-3.5, 1.0), nlpar = "etaai", class = "b"),
  prior(normal(-3.5, 1.0), nlpar = "etaaj", class = "b"),
  prior(exponential(1), class = "sd", nlpar = "etalam"),
  prior(exponential(1), class = "sigma")
)

ctrl <- list(adapt_delta = 0.99, max_treedepth = 15)

############################
## 4. Model formulas
############################

## M1: lambda varies by environment; alpha constant by species
bf_bh_m1 <- bf(
  fitness_log ~ log(
    exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
  ),
  etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
    focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## M2: lambda varies by environment; alpha varies by species + Site + Cov
bf_bh_m2 <- bf(
  fitness_log ~ log(
    exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
  ),
  etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
    focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
  etaai ~ 0 + focsp + focsp:Site + focsp:Cov,
  etaaj ~ 0 + focsp + focsp:Site + focsp:Cov,
  nl = TRUE
)

## M3: lambda varies by environment; alpha varies by species + Cov only
bf_bh_m3 <- bf(
  fitness_log ~ log(
    exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
  ),
  etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
    focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
  etaai ~ 0 + focsp + focsp:Cov,
  etaaj ~ 0 + focsp + focsp:Cov,
  nl = TRUE
)

############################
## 5. Helper functions
############################

fit_or_load_brms <- function(formula_obj, data, family, prior, seed, file_base, control) {
  rds_file <- paste0(file_base, ".rds")
  
  if (file.exists(rds_file)) {
    cat("\nLoading existing model:", rds_file, "\n")
    fit <- readRDS(rds_file)
  } else {
    cat("\nFitting model and saving to:", rds_file, "\n")
    fit <- brm(
      formula = formula_obj,
      data = data,
      family = family,
      prior = prior,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      control = control,
      seed = seed,
      save_pars = save_pars(all = TRUE),
      file = file_base
    )
  }
  
  fit
}

loo_or_load <- function(fit, loo_file, use_reloo = FALSE) {
  if (file.exists(loo_file)) {
    cat("\nLoading existing LOO object:", loo_file, "\n")
    loo_obj <- readRDS(loo_file)
  } else {
    cat("\nComputing LOO and saving to:", loo_file, "\n")
    if (use_reloo) {
      loo_obj <- loo(fit, reloo = TRUE)
    } else {
      loo_obj <- loo(fit, moment_match = TRUE)
    }
    saveRDS(loo_obj, loo_file)
  }
  
  loo_obj
}

############################
## 6. Fit or load models
############################
fit_bh_m1 <- fit_or_load_brms(
  formula_obj = bf_bh_m1,
  data        = dat_bh,
  family      = bh_family,
  prior       = priors_common,
  seed        = 1234,
  file_base   = "fit_bh_m1_loggaussian",
  control     = ctrl
)

fit_bh_m2 <- fit_or_load_brms(
  formula_obj = bf_bh_m2,
  data        = dat_bh,
  family      = bh_family,
  prior       = priors_common,
  seed        = 1235,
  file_base   = "fit_bh_m2_loggaussian",
  control     = ctrl
)

fit_bh_m3 <- fit_or_load_brms(
  formula_obj = bf_bh_m3,
  data        = dat_bh,
  family      = bh_family,
  prior       = priors_common,
  seed        = 1236,
  file_base   = "fit_bh_m3_loggaussian",
  control     = ctrl
)

############################
## 7. LOO: compute or load
############################
loo_bh_m1 <- loo_or_load(
  fit      = fit_bh_m1,
  loo_file = "loo_bh_m1_loggaussian.rds",
  use_reloo = FALSE
)

loo_bh_m2 <- loo_or_load(
  fit      = fit_bh_m2,
  loo_file = "loo_bh_m2_loggaussian.rds",
  use_reloo = FALSE
)

loo_bh_m3 <- loo_or_load(
  fit      = fit_bh_m3,
  loo_file = "loo_bh_m3_loggaussian.rds",
  use_reloo = FALSE
)

cat("\n====================\nLOO summaries\n====================\n")
print(loo_bh_m1)
print(loo_bh_m2)
print(loo_bh_m3)

############################
## 8. LOO comparison and final model choice
############################
loo_list <- list(
  m1_species  = loo_bh_m1,
  m2_site_cov = loo_bh_m2,
  m3_cov      = loo_bh_m3
)

loo_tab <- loo_compare(loo_list)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

write.csv(
  as.data.frame(loo_tab) %>%
    tibble::rownames_to_column("model"),
  "BH_loggaussian_model_comparison_LOO.csv",
  row.names = FALSE
)

## Set final model manually here if needed
selected_model <- "m1_species"
selected_fit_file <- switch(
  selected_model,
  m1_species  = "fit_bh_m1_loggaussian.rds",
  m2_site_cov = "fit_bh_m2_loggaussian.rds",
  m3_cov      = "fit_bh_m3_loggaussian.rds"
)

model_selection_info <- list(
  selected_model    = selected_model,
  selected_fit_file = selected_fit_file,
  site_order        = c("BEN", "GH", "NAM", "PJ", "CS")
)

saveRDS(model_selection_info, "model_selection_info.rds")

############################
## 9. Save summaries
############################
sink("BH_loggaussian_model_summaries.txt")
cat("\n====================\nModel 1 summary\n====================\n")
print(summary(fit_bh_m1))

cat("\n====================\nModel 2 summary\n====================\n")
print(summary(fit_bh_m2))

cat("\n====================\nModel 3 summary\n====================\n")
print(summary(fit_bh_m3))

cat("\n====================\nLOO summaries\n====================\n")
print(loo_bh_m1)
print(loo_bh_m2)
print(loo_bh_m3)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

cat("\n====================\nSelected model\n====================\n")
print(model_selection_info)
sink()

############################
## 10. Posterior predictive checks
############################
pdf("BH_loggaussian_pp_checks.pdf", width = 8, height = 6)
print(pp_check(fit_bh_m1, ndraws = 100))
print(pp_check(fit_bh_m2, ndraws = 100))
print(pp_check(fit_bh_m3, ndraws = 100))
dev.off()

############################
## 11. Final object in console
############################
loo_tab
model_selection_info