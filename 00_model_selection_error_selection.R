############################################################
## 00_family_error_selection_with_germination.R
##
## Purpose:
## 1. Prepare competition analysis data
## 2. Load and validate germination priors
## 3. Compare Beverton-Holt vs Ricker crossed with
##    four observation/error models
## 4. Save fits, LOO objects, and a clear record of the
##    selected competition-family + error-model combination
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
## 1. Set output folders
############################
base_dir <- "00_family_error_selection"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "fits"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "loo"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "selection"), showWarnings = FALSE, recursive = TRUE)

############################
## 2. Read and prepare competition data
############################
dat0 <- read.csv("fullWAdata.csv")

dat_comp <- dat0 %>%
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
    intra_sc      = intra,
    inter_sc      = inter,
    fac_sc        = fac,
    fitness_raw   = fitness,
    fitness_log   = log(fitness + 1),
    fitness_count = pmax(0, round(fitness))
  ) %>%
  droplevels()

write.csv(
  dat_comp,
  file.path(base_dir, "data", "analysis_data_family_error_compare.csv"),
  row.names = FALSE
)
saveRDS(
  dat_comp,
  file.path(base_dir, "data", "analysis_data_family_error_compare.rds")
)

############################
## 3. Load and validate germination priors
############################
## Germination is NOT used in this stage-0 fitting likelihood.
## It is loaded here so the selected model object can carry
## the correct germination prior file forward into downstream
## coexistence calculations.

germination_priors_file <- file.path(
  "00_germination_priors", "selection", "germination_priors.rds"
)

if (!file.exists(germination_priors_file)) {
  stop(
    "Could not find germination prior file:\n  ",
    germination_priors_file,
    "\nRun 00_germination_priors.R first."
  )
}

germ_info <- readRDS(germination_priors_file)

if (!"site_priors" %in% names(germ_info)) {
  stop("germination_priors.rds does not contain an element named 'site_priors'.")
}

germ_priors <- germ_info$site_priors %>%
  as_tibble() %>%
  mutate(
    Site  = as.character(Site),
    focsp = as.character(focsp)
  )

required_germ_cols <- c("Site", "focsp", "g_alpha_final", "g_beta_final", "g_mean_final")
missing_germ_cols <- setdiff(required_germ_cols, names(germ_priors))

if (length(missing_germ_cols) > 0) {
  stop(
    "The germination priors object is missing required columns:\n  ",
    paste(missing_germ_cols, collapse = ", ")
  )
}

expected_site_sp <- expand.grid(
  Site  = levels(dat_comp$Site),
  focsp = levels(dat_comp$focsp),
  stringsAsFactors = FALSE
) %>%
  as_tibble()

germ_check <- expected_site_sp %>%
  left_join(germ_priors, by = c("Site", "focsp"))

if (any(is.na(germ_check$g_alpha_final) | is.na(germ_check$g_beta_final))) {
  stop(
    "Missing germination priors for at least one Site x focsp combination used in fullWAdata.csv.\n",
    "Check 00_germination_priors.R output."
  )
}

write.csv(
  germ_check,
  file.path(base_dir, "data", "germination_site_priors_used.csv"),
  row.names = FALSE
)
saveRDS(
  germ_check,
  file.path(base_dir, "data", "germination_site_priors_used.rds")
)

############################
## 4. Quick data checks
############################
cat("\nSummary of raw fitness:\n")
print(summary(dat_comp$fitness_raw))

cat("\nCheck whether fitness is integer-valued:\n")
print(all(abs(dat_comp$fitness_raw - round(dat_comp$fitness_raw)) < 1e-8))

cat("\nSummary of rounded count fitness:\n")
print(summary(dat_comp$fitness_count))

capture.output(
  {
    cat("\nSummary of raw fitness:\n")
    print(summary(dat_comp$fitness_raw))
    cat("\nCheck whether fitness is integer-valued:\n")
    print(all(abs(dat_comp$fitness_raw - round(dat_comp$fitness_raw)) < 1e-8))
    cat("\nSummary of rounded count fitness:\n")
    print(summary(dat_comp$fitness_count))
    cat("\nGermination priors loaded from:\n")
    print(normalizePath(germination_priors_file))
  },
  file = file.path(base_dir, "summaries", "data_checks.txt")
)

############################
## 5. Priors
############################
## Same weakly regularizing priors used for BH and Ricker.
## etalam is on log(lambda)
## etaai and etaaj are on log(alpha)

base_priors <- c(
  prior(normal(1.5, 1), class = "b", nlpar = "etalam"),
  prior(normal(-3, 1), class = "b", nlpar = "etaai"),
  prior(normal(-3, 1), class = "b", nlpar = "etaaj"),
  prior(exponential(1), class = "sd", nlpar = "etalam")
)

priors_gaussian <- c(
  base_priors,
  prior(exponential(1), class = "sigma")
)

priors_poisson <- base_priors

priors_negbin <- c(
  base_priors,
  prior(exponential(1), class = "shape")
)

priors_loggaussian <- c(
  base_priors,
  prior(exponential(1), class = "sigma")
)

############################
## 6. Model formulas
############################
## Beverton-Holt:
##   F = lambda / (1 + alpha_ii * N_i + alpha_ij * N_j)
##
## Ricker:
##   F = lambda * exp(- alpha_ii * N_i - alpha_ij * N_j)

## ---- Beverton-Holt + Gaussian raw
form_bh_gaussian <- bf(
  fitness_raw ~ exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + Poisson count
form_bh_poisson <- bf(
  fitness_count ~ exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + NegBin count
form_bh_negbin <- bf(
  fitness_count ~ exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + Gaussian logfit
form_bh_loggaussian <- bf(
  fitness_log ~ log(
    exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
  ),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Gaussian raw
form_ricker_gaussian <- bf(
  fitness_raw ~ exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Poisson count
form_ricker_poisson <- bf(
  fitness_count ~ exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + NegBin count
form_ricker_negbin <- bf(
  fitness_count ~ exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Gaussian logfit
form_ricker_loggaussian <- bf(
  fitness_log ~ log(
    exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
  ),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    (1 | BLK),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

############################
## 7. Fit settings
############################
fit_args <- list(
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = min(4, max(1, parallel::detectCores() - 1)),
  refresh = 50,
  control = list(adapt_delta = 0.995, max_treedepth = 14),
  save_pars = save_pars(all = TRUE)
)

############################
## 8. Helper: fit or load brms model
############################
fit_or_load <- function(formula, data, family, priors, model_name, seed, base_dir) {
  
  fit_rds  <- file.path(base_dir, "fits", paste0(model_name, ".rds"))
  fit_file <- file.path(base_dir, "fits", model_name)
  
  if (file.exists(fit_rds)) {
    cat("\nLoading existing fit:", model_name, "\n")
    fit <- readRDS(fit_rds)
  } else {
    cat("\nFitting model:", model_name, "\n")
    fit <- brm(
      formula = formula,
      data = data,
      family = family,
      prior = priors,
      chains = fit_args$chains,
      iter = fit_args$iter,
      warmup = fit_args$warmup,
      cores = fit_args$cores,
      refresh = fit_args$refresh,
      control = fit_args$control,
      save_pars = fit_args$save_pars,
      seed = seed,
      file = fit_file,
      file_refit = "on_change"
    )
    saveRDS(fit, fit_rds)
  }
  
  try(
    capture.output(
      print(summary(fit), digits = 2),
      file = file.path(base_dir, "summaries", paste0(model_name, "_summary.txt"))
    ),
    silent = TRUE
  )
  
  try({
    png(
      file.path(base_dir, "figures", paste0(model_name, "_ppcheck.png")),
      width = 1400, height = 1000, res = 160
    )
    print(pp_check(fit))
    dev.off()
  }, silent = TRUE)
  
  try(
    writeLines(
      variables(fit),
      con = file.path(base_dir, "summaries", paste0(model_name, "_variables.txt"))
    ),
    silent = TRUE
  )
  
  invisible(fit)
}

############################
## 9. Helper: compute or load LOO
############################
loo_or_load <- function(fit, loo_name, base_dir) {
  
  loo_rds <- file.path(base_dir, "loo", paste0(loo_name, ".rds"))
  
  if (file.exists(loo_rds)) {
    cat("\nLoading existing LOO:", loo_name, "\n")
    loo_obj <- readRDS(loo_rds)
  } else {
    cat("\nComputing LOO:", loo_name, "\n")
    loo_obj <- loo(fit, moment_match = TRUE)
    saveRDS(loo_obj, loo_rds)
  }
  
  invisible(loo_obj)
}

############################
## 10. Model specification table
############################
model_specs <- list(
  list(
    key = "BH_gaussian_raw",
    growth_family = "BH",
    error_model = "gaussian_raw",
    formula = form_bh_gaussian,
    family = gaussian(),
    priors = priors_gaussian,
    fit_name = "fit_BH_gaussian_raw",
    loo_name = "loo_BH_gaussian_raw",
    seed = 101
  ),
  list(
    key = "BH_poisson_count",
    growth_family = "BH",
    error_model = "poisson_count",
    formula = form_bh_poisson,
    family = poisson(link = "identity"),
    priors = priors_poisson,
    fit_name = "fit_BH_poisson_count",
    loo_name = "loo_BH_poisson_count",
    seed = 202
  ),
  list(
    key = "BH_negbin_count",
    growth_family = "BH",
    error_model = "negbin_count",
    formula = form_bh_negbin,
    family = negbinomial(link = "identity"),
    priors = priors_negbin,
    fit_name = "fit_BH_negbin_count",
    loo_name = "loo_BH_negbin_count",
    seed = 303
  ),
  list(
    key = "BH_gaussian_logfit",
    growth_family = "BH",
    error_model = "gaussian_logfit",
    formula = form_bh_loggaussian,
    family = gaussian(),
    priors = priors_loggaussian,
    fit_name = "fit_BH_gaussian_logfit",
    loo_name = "loo_BH_gaussian_logfit",
    seed = 404
  ),
  list(
    key = "Ricker_gaussian_raw",
    growth_family = "Ricker",
    error_model = "gaussian_raw",
    formula = form_ricker_gaussian,
    family = gaussian(),
    priors = priors_gaussian,
    fit_name = "fit_Ricker_gaussian_raw",
    loo_name = "loo_Ricker_gaussian_raw",
    seed = 505
  ),
  list(
    key = "Ricker_poisson_count",
    growth_family = "Ricker",
    error_model = "poisson_count",
    formula = form_ricker_poisson,
    family = poisson(link = "identity"),
    priors = priors_poisson,
    fit_name = "fit_Ricker_poisson_count",
    loo_name = "loo_Ricker_poisson_count",
    seed = 606
  ),
  list(
    key = "Ricker_negbin_count",
    growth_family = "Ricker",
    error_model = "negbin_count",
    formula = form_ricker_negbin,
    family = negbinomial(link = "identity"),
    priors = priors_negbin,
    fit_name = "fit_Ricker_negbin_count",
    loo_name = "loo_Ricker_negbin_count",
    seed = 707
  ),
  list(
    key = "Ricker_gaussian_logfit",
    growth_family = "Ricker",
    error_model = "gaussian_logfit",
    formula = form_ricker_loggaussian,
    family = gaussian(),
    priors = priors_loggaussian,
    fit_name = "fit_Ricker_gaussian_logfit",
    loo_name = "loo_Ricker_gaussian_logfit",
    seed = 808
  )
)

spec_lookup <- setNames(model_specs, vapply(model_specs, `[[`, character(1), "key"))

model_spec_df <- tibble(
  model = vapply(model_specs, `[[`, character(1), "key"),
  growth_family = vapply(model_specs, `[[`, character(1), "growth_family"),
  error_model   = vapply(model_specs, `[[`, character(1), "error_model"),
  fit_name      = vapply(model_specs, `[[`, character(1), "fit_name"),
  loo_name      = vapply(model_specs, `[[`, character(1), "loo_name"),
  seed          = vapply(model_specs, `[[`, numeric(1), "seed")
)

write.csv(
  model_spec_df,
  file.path(base_dir, "selection", "candidate_family_error_models.csv"),
  row.names = FALSE
)

############################
## 11. Fit/load all candidate models
############################
fit_list <- list()

for (spec in model_specs) {
  fit_list[[spec$key]] <- fit_or_load(
    formula = spec$formula,
    data = dat_comp,
    family = spec$family,
    priors = spec$priors,
    model_name = spec$fit_name,
    seed = spec$seed,
    base_dir = base_dir
  )
}

############################
## 12. Compute/load LOO for all models
############################
loo_list <- list()

for (spec in model_specs) {
  loo_list[[spec$key]] <- loo_or_load(
    fit = fit_list[[spec$key]],
    loo_name = spec$loo_name,
    base_dir = base_dir
  )
}

############################
## 13. LOO comparison
############################
loo_tab <- loo_compare(loo_list)

print(loo_tab)

capture.output(
  print(loo_tab),
  file = file.path(base_dir, "selection", "loo_compare_family_error_models.txt")
)

write.csv(
  as.data.frame(loo_tab) %>%
    tibble::rownames_to_column("model"),
  file.path(base_dir, "selection", "loo_compare_family_error_models.csv"),
  row.names = FALSE
)

weights <- loo_model_weights(
  loo_list,
  method = "stacking"
)

weights_df <- tibble(
  model = names(weights),
  stacking_weight = as.numeric(weights)
) %>%
  left_join(model_spec_df, by = "model") %>%
  arrange(desc(stacking_weight))

write.csv(
  weights_df,
  file.path(base_dir, "selection", "model_weights_family_error_models.csv"),
  row.names = FALSE
)

############################
## 14. Record selected family + error model
############################
best_model <- rownames(loo_tab)[1]
best_spec  <- spec_lookup[[best_model]]

selection_info <- list(
  best_model = best_model,
  best_growth_family = best_spec$growth_family,
  best_error_model = best_spec$error_model,
  best_fit_rds = file.path(base_dir, "fits", paste0(best_spec$fit_name, ".rds")),
  best_fit_name = best_spec$fit_name,
  comparison_csv = file.path(base_dir, "selection", "loo_compare_family_error_models.csv"),
  weights_csv = file.path(base_dir, "selection", "model_weights_family_error_models.csv"),
  candidate_models_csv = file.path(base_dir, "selection", "candidate_family_error_models.csv"),
  data_rds = file.path(base_dir, "data", "analysis_data_family_error_compare.rds"),
  germination_priors_rds = normalizePath(germination_priors_file),
  germination_site_priors_rds = file.path(base_dir, "data", "germination_site_priors_used.rds"),
  note = "Germination priors are loaded and validated here but are NOT used in stage-0 fitting; use them downstream for lambda_eff = g * lambda."
)

saveRDS(
  selection_info,
  file.path(base_dir, "selection", "family_error_model_selection_info.rds")
)

writeLines(
  best_model,
  con = file.path(base_dir, "selection", "best_family_error_model.txt")
)

note_lines <- c(
  paste("Selected model:", best_model),
  paste("Selected competition family:", best_spec$growth_family),
  paste("Selected error model:", best_spec$error_model),
  "",
  "Interpretation notes:",
  "1. Use the LOO comparison table for formal family + error model selection.",
  "2. Use pp_check plots and summaries as complementary diagnostics.",
  "3. If models are very close, prefer the biologically appropriate and stable model.",
  "4. Count-based models here were fit to rounded fitness_count; note this in methods if discussed.",
  "5. Germination priors were loaded and validated here, but not used in the stage-0 likelihood.",
  "6. Use germination downstream in coexistence calculations via lambda_eff = g * lambda."
)

writeLines(
  note_lines,
  con = file.path(base_dir, "selection", "model_selection_notes.txt")
)

############################
## 15. Save session info
############################
capture.output(
  sessionInfo(),
  file = file.path(base_dir, "selection", "sessionInfo.txt")
)

############################
## 16. Finish
############################
cat("\nDone.\n")
cat("Outputs saved in:\n", normalizePath(base_dir), "\n")
cat("Selected model:\n", best_model, "\n")
cat("Selected competition family:\n", best_spec$growth_family, "\n")
cat("Selected error model:\n", best_spec$error_model, "\n")
