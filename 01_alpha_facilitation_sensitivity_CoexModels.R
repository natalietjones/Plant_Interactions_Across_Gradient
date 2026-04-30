############################################################
## alpha_facilitation_sensitivity.R
##
## Purpose:
## Compare the selected baseline alpha-structure model (m1_species)
## against targeted sensitivity models that allow facilitator density
## (fac_sc) to modify alpha_ij only, or both alpha_ii and alpha_ij.
##
## Biological motivation:
## - Baseline m1 assumes facilitators act through lambda only.
## - m4 tests whether facilitators modify interspecific competition.
## - m5 tests whether facilitators modify both intra- and interspecific
##   competition.
##
## Recommended use:
## Run this AFTER 01_alpha_structure_selection_with_germination.R has
## completed and selected m1_species as the primary model.
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
## 1. User options
############################
## Leave NULL to use the previously selected family
family_override <- NULL

############################
## 2. Paths
############################
base_stage1_dir <- "01_alpha_structure_selection"
base_dir        <- "01b_alpha_facilitation_sensitivity"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "fits"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "loo"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "selection"), showWarnings = FALSE, recursive = TRUE)

############################
## 3. Load previous selection info and data
############################
selection_info <- readRDS(
  file.path(base_stage1_dir, "selection", "alpha_structure_selection_info.rds")
)

dat_stage1 <- readRDS(selection_info$stage1_data_rds)

family_to_use <- if (is.null(family_override)) {
  selection_info$selected_growth_family
} else {
  family_override
}

if (!family_to_use %in% c("BH", "Ricker")) {
  stop("family_to_use must be either 'BH' or 'Ricker'.")
}

if (!identical(selection_info$selected_structure_model, "m1_species")) {
  warning(
    "The previously selected model was not m1_species. \n",
    "This script will still run, but its main purpose is to compare m1_species \n",
    "to facilitator-on-alpha sensitivity models."
  )
}

cat("\nLoaded selection info:\n")
print(selection_info)

cat("\nFamily used in sensitivity analysis:\n")
print(family_to_use)

############################
## 4. Family, priors, settings
############################
model_family <- gaussian()

priors_common <- c(
  prior(normal(0, 1.5), nlpar = "etalam", class = "b"),
  prior(normal(-3.5, 1.0), nlpar = "etaai", class = "b"),
  prior(normal(-3.5, 1.0), nlpar = "etaaj", class = "b"),
  prior(exponential(1), class = "sd", nlpar = "etalam"),
  prior(exponential(1), class = "sigma")
)

ctrl <- list(adapt_delta = 0.99, max_treedepth = 15)

############################
## 5. Model formulas
############################
## Baseline m1: facilitators act through lambda only
## m4_fac_inter: facilitators act through lambda + alpha_ij
## m5_fac_both : facilitators act through lambda + alpha_ii + alpha_ij

if (family_to_use == "BH") {

  bf_m1 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m4 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )

  bf_m5 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:fac_sc,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )

} else {

  bf_m1 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m4 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )

  bf_m5 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:fac_sc,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )
}

############################
## 6. Helper functions
############################
fit_or_load_brms <- function(formula_obj, data, family, prior, seed, file_base, control, base_dir) {
  rds_file <- file.path(base_dir, "fits", paste0(file_base, ".rds"))
  fit_file <- file.path(base_dir, "fits", file_base)

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
      file = fit_file,
      file_refit = "on_change"
    )
    saveRDS(fit, rds_file)
  }

  capture.output(
    print(summary(fit), digits = 2),
    file = file.path(base_dir, "summaries", paste0(file_base, "_summary.txt"))
  )

  try({
    png(
      file.path(base_dir, "figures", paste0(file_base, "_ppcheck.png")),
      width = 1400, height = 1000, res = 160
    )
    print(pp_check(fit))
    dev.off()
  }, silent = TRUE)

  invisible(fit)
}

loo_or_load <- function(fit, loo_file, base_dir) {
  loo_path <- file.path(base_dir, "loo", loo_file)

  if (file.exists(loo_path)) {
    cat("\nLoading existing LOO object:", loo_path, "\n")
    loo_obj <- readRDS(loo_path)
  } else {
    cat("\nComputing LOO and saving to:", loo_path, "\n")
    loo_obj <- loo(fit, moment_match = TRUE)
    saveRDS(loo_obj, loo_path)
  }

  loo_obj
}

############################
## 7. Fit/load baseline m1 from previous stage or refit locally
############################
fit_m1_path <- selection_info$selected_fit_file

if (file.exists(fit_m1_path)) {
  cat("\nLoading baseline selected m1 model from:\n", fit_m1_path, "\n")
  fit_m1 <- readRDS(fit_m1_path)
} else {
  warning(
    "Could not find the saved selected m1 fit at:\n", fit_m1_path,
    "\nRefitting m1 in this sensitivity-analysis directory."
  )
  fit_m1 <- fit_or_load_brms(
    formula_obj = bf_m1,
    data        = dat_stage1,
    family      = model_family,
    prior       = priors_common,
    seed        = 2001,
    file_base   = paste0("fit_", family_to_use, "_m1_baseline_loggaussian"),
    control     = ctrl,
    base_dir    = base_dir
  )
}

## Save a local summary of baseline m1 regardless
capture.output(
  print(summary(fit_m1), digits = 2),
  file = file.path(base_dir, "summaries", paste0("fit_", family_to_use, "_m1_baseline_loggaussian_summary.txt"))
)

############################
## 8. Fit/load new facilitator-on-alpha models
############################
fit_base_m4 <- paste0("fit_", family_to_use, "_m4_fac_inter_loggaussian")
fit_base_m5 <- paste0("fit_", family_to_use, "_m5_fac_both_loggaussian")

fit_m4 <- fit_or_load_brms(
  formula_obj = bf_m4,
  data        = dat_stage1,
  family      = model_family,
  prior       = priors_common,
  seed        = 2004,
  file_base   = fit_base_m4,
  control     = ctrl,
  base_dir    = base_dir
)

fit_m5 <- fit_or_load_brms(
  formula_obj = bf_m5,
  data        = dat_stage1,
  family      = model_family,
  prior       = priors_common,
  seed        = 2005,
  file_base   = fit_base_m5,
  control     = ctrl,
  base_dir    = base_dir
)

############################
## 9. LOO comparison
############################
loo_m1 <- loo_or_load(
  fit      = fit_m1,
  loo_file = paste0("loo_", family_to_use, "_m1_baseline_loggaussian.rds"),
  base_dir = base_dir
)

loo_m4 <- loo_or_load(
  fit      = fit_m4,
  loo_file = paste0("loo_", family_to_use, "_m4_fac_inter_loggaussian.rds"),
  base_dir = base_dir
)

loo_m5 <- loo_or_load(
  fit      = fit_m5,
  loo_file = paste0("loo_", family_to_use, "_m5_fac_both_loggaussian.rds"),
  base_dir = base_dir
)

loo_list <- list(
  m1_lambda_only = loo_m1,
  m4_fac_on_alphaij = loo_m4,
  m5_fac_on_alphaii_alphaij = loo_m5
)

loo_tab <- loo_compare(loo_list)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

write.csv(
  as.data.frame(loo_tab) %>% tibble::rownames_to_column("model"),
  file.path(base_dir, "selection", paste0(family_to_use, "_facilitation_alpha_sensitivity_LOO.csv")),
  row.names = FALSE
)

selected_model <- rownames(loo_tab)[1]

selection_out <- list(
  family_used = family_to_use,
  baseline_model = "m1_lambda_only",
  candidate_models = c(
    "m1_lambda_only",
    "m4_fac_on_alphaij",
    "m5_fac_on_alphaii_alphaij"
  ),
  selected_model = selected_model,
  baseline_fit_file = fit_m1_path,
  m4_fit_file = file.path(base_dir, "fits", paste0(fit_base_m4, ".rds")),
  m5_fit_file = file.path(base_dir, "fits", paste0(fit_base_m5, ".rds")),
  note = paste(
    "m1 tests facilitator effects through lambda only;",
    "m4 tests facilitator effects through lambda plus alpha_ij;",
    "m5 tests facilitator effects through lambda plus alpha_ii and alpha_ij."
  )
)

saveRDS(
  selection_out,
  file.path(base_dir, "selection", paste0("facilitation_alpha_sensitivity_selection_", family_to_use, ".rds"))
)

############################
## 10. Save text summary
############################
summary_txt <- file.path(
  base_dir,
  "summaries",
  paste0(family_to_use, "_facilitation_alpha_sensitivity_summary.txt")
)

sink(summary_txt)
cat("\n====================\nFamily used\n====================\n")
print(family_to_use)

cat("\n====================\nBaseline m1 summary\n====================\n")
print(summary(fit_m1))

cat("\n====================\nm4 summary: fac_sc on alpha_ij\n====================\n")
print(summary(fit_m4))

cat("\n====================\nm5 summary: fac_sc on alpha_ii and alpha_ij\n====================\n")
print(summary(fit_m5))

cat("\n====================\nLOO summaries\n====================\n")
print(loo_m1)
print(loo_m4)
print(loo_m5)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

cat("\n====================\nSelected model\n====================\n")
print(selection_out)
sink()

############################
## 11. Console output
############################
cat("\nDone.\n")
cat("Family used:\n", family_to_use, "\n")
cat("Selected model:\n", selected_model, "\n")
cat("Outputs saved in:\n", normalizePath(base_dir), "\n")

loo_tab
selection_out
