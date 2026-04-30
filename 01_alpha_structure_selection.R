############################################################
## 01_alpha_structure_selection_with_germination.R
##
## Purpose:
## 1. Load selected family+error outputs from 00_family_error_selection/
## 2. Load and validate germination priors
## 3. Fit/load candidate alpha-structure models under the selected
##    log-Gaussian competition family branch (BH primary by default)
## 4. Compare models with LOO
## 5. Save structure-selection outputs for downstream coexistence
##    calculations that will incorporate germination
##
## Notes:
## - Germination priors are loaded and saved forward here, but NOT
##   inserted into the stage-1 likelihood.
## - This script defaults to the best family from stage 0.
## - To run the Ricker sensitivity version, set:
##     family_override <- "Ricker"
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
## Leave NULL to use the best family from stage 0
## Set to "Ricker" to run the sensitivity branch
family_override <- NULL

############################
## 2. Set input/output folders
############################
error_dir <- "00_family_error_selection"
base_dir  <- "01_alpha_structure_selection"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "fits"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "loo"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "selection"), showWarnings = FALSE, recursive = TRUE)

############################
## 3. Load selected family+error info and prepared data
############################
selection_info <- readRDS(
  file.path(error_dir, "selection", "family_error_model_selection_info.rds")
)

dat_comp <- readRDS(
  file.path(error_dir, "data", "analysis_data_family_error_compare.rds")
)

cat("\nLoaded family+error selection info:\n")
print(selection_info)

family_to_use <- if (is.null(family_override)) {
  selection_info$best_growth_family
} else {
  family_override
}

error_to_use <- selection_info$best_error_model

cat("\nFamily branch to use in stage 1:\n")
print(family_to_use)

cat("\nError model from stage 0:\n")
print(error_to_use)

if (!identical(error_to_use, "gaussian_logfit")) {
  stop(
    paste0(
      "The selected stage-0 error model is '", error_to_use, "', not 'gaussian_logfit'.\n",
      "This stage-1 script currently supports the log-Gaussian branch only."
    )
  )
}

if (!family_to_use %in% c("BH", "Ricker")) {
  stop("family_to_use must be either 'BH' or 'Ricker'.")
}

############################
## 4. Load and validate germination priors
############################
germination_priors_file <- selection_info$germination_priors_rds

if (!file.exists(germination_priors_file)) {
  stop(
    "Could not find germination priors file:\n  ",
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

## Save a copy used in this stage
write.csv(
  germ_check,
  file.path(base_dir, "data", "germination_site_priors_used.csv"),
  row.names = FALSE
)
saveRDS(
  germ_check,
  file.path(base_dir, "data", "germination_site_priors_used.rds")
)

## Join germination means to the data for bookkeeping only
## (not used in model formulas at this stage)
site_levels <- levels(dat_comp$Site)
focsp_levels <- levels(dat_comp$focsp)

dat_stage1 <- dat_comp %>%
  mutate(
    Site  = as.character(Site),
    focsp = as.character(focsp)
  ) %>%
  left_join(
    germ_check %>%
      select(Site, focsp, g_alpha_final, g_beta_final, g_mean_final),
    by = c("Site", "focsp")
  ) %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    focsp = factor(focsp, levels = focsp_levels)
  )

if (any(is.na(dat_stage1$g_mean_final))) {
  stop("Germination join failed for some rows in dat_stage1.")
}

write.csv(
  dat_stage1,
  file.path(base_dir, "data", "analysis_data_stage1_with_germination.csv"),
  row.names = FALSE
)
saveRDS(
  dat_stage1,
  file.path(base_dir, "data", "analysis_data_stage1_with_germination.rds")
)

############################
## 5. Quick checks
############################
cat("\nSummary of log(fitness + 1):\n")
print(summary(dat_stage1$fitness_log))

cat("\nLevels:\n")
print(levels(dat_stage1$focsp))
print(levels(dat_stage1$Site))
print(levels(dat_stage1$Cov))

cat("\nGermination priors loaded:\n")
print(
  germ_check %>%
    select(Site, focsp, g_mean_final) %>%
    arrange(focsp, Site)
)

capture.output(
  {
    cat("\nFamily branch used:\n")
    print(family_to_use)
    cat("\nError model used:\n")
    print(error_to_use)
    cat("\nSummary of log(fitness + 1):\n")
    print(summary(dat_stage1$fitness_log))
    cat("\nLevels:\n")
    print(levels(dat_stage1$focsp))
    print(levels(dat_stage1$Site))
    print(levels(dat_stage1$Cov))
    cat("\nGermination priors loaded:\n")
    print(
      germ_check %>%
        select(Site, focsp, g_mean_final) %>%
        arrange(focsp, Site)
    )
  },
  file = file.path(base_dir, "summaries", paste0("data_checks_", family_to_use, "_loggaussian.txt"))
)

############################
## 6. Family, priors, settings
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
## 7. Model formulas
############################
## M1: lambda varies by environment; alpha constant by species
## M2: lambda varies by environment; alpha varies by species + Site + Cov
## M3: lambda varies by environment; alpha varies by species + Cov only

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
  
  bf_m2 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:Site + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Site + focsp:Cov,
    nl = TRUE
  )
  
  bf_m3 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Cov,
    nl = TRUE
  )
  
} else if (family_to_use == "Ricker") {
  
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
  
  bf_m2 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:Site + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Site + focsp:Cov,
    nl = TRUE
  )
  
  bf_m3 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | BLK),
    etaai ~ 0 + focsp + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Cov,
    nl = TRUE
  )
}

############################
## 8. Helper functions
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

loo_or_load <- function(fit, loo_file, use_reloo = FALSE, base_dir) {
  loo_path <- file.path(base_dir, "loo", loo_file)
  
  if (file.exists(loo_path)) {
    cat("\nLoading existing LOO object:", loo_path, "\n")
    loo_obj <- readRDS(loo_path)
  } else {
    cat("\nComputing LOO and saving to:", loo_path, "\n")
    if (use_reloo) {
      loo_obj <- loo(fit, reloo = TRUE)
    } else {
      loo_obj <- loo(fit, moment_match = TRUE)
    }
    saveRDS(loo_obj, loo_path)
  }
  
  loo_obj
}

############################
## 9. File tags
############################
family_tag <- family_to_use
error_tag  <- "loggaussian"

fit_base_m1 <- paste0("fit_", family_tag, "_m1_", error_tag)
fit_base_m2 <- paste0("fit_", family_tag, "_m2_", error_tag)
fit_base_m3 <- paste0("fit_", family_tag, "_m3_", error_tag)

loo_file_m1 <- paste0("loo_", family_tag, "_m1_", error_tag, ".rds")
loo_file_m2 <- paste0("loo_", family_tag, "_m2_", error_tag, ".rds")
loo_file_m3 <- paste0("loo_", family_tag, "_m3_", error_tag, ".rds")

############################
## 10. Fit or load models
############################
fit_m1 <- fit_or_load_brms(
  formula_obj = bf_m1,
  data        = dat_stage1,
  family      = model_family,
  prior       = priors_common,
  seed        = 1234,
  file_base   = fit_base_m1,
  control     = ctrl,
  base_dir    = base_dir
)

fit_m2 <- fit_or_load_brms(
  formula_obj = bf_m2,
  data        = dat_stage1,
  family      = model_family,
  prior       = priors_common,
  seed        = 1235,
  file_base   = fit_base_m2,
  control     = ctrl,
  base_dir    = base_dir
)

fit_m3 <- fit_or_load_brms(
  formula_obj = bf_m3,
  data        = dat_stage1,
  family      = model_family,
  prior       = priors_common,
  seed        = 1236,
  file_base   = fit_base_m3,
  control     = ctrl,
  base_dir    = base_dir
)

############################
## 11. LOO: compute or load
############################
loo_m1 <- loo_or_load(
  fit       = fit_m1,
  loo_file  = loo_file_m1,
  use_reloo = FALSE,
  base_dir  = base_dir
)

loo_m2 <- loo_or_load(
  fit       = fit_m2,
  loo_file  = loo_file_m2,
  use_reloo = FALSE,
  base_dir  = base_dir
)

loo_m3 <- loo_or_load(
  fit       = fit_m3,
  loo_file  = loo_file_m3,
  use_reloo = FALSE,
  base_dir  = base_dir
)

cat("\n====================\nLOO summaries\n====================\n")
print(loo_m1)
print(loo_m2)
print(loo_m3)

############################
## 12. LOO comparison and final model choice
############################
loo_list <- list(
  m1_species  = loo_m1,
  m2_site_cov = loo_m2,
  m3_cov      = loo_m3
)

loo_tab <- loo_compare(loo_list)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

loo_csv_file <- file.path(
  base_dir,
  "selection",
  paste0(family_tag, "_", error_tag, "_structure_comparison_LOO.csv")
)

write.csv(
  as.data.frame(loo_tab) %>%
    tibble::rownames_to_column("model"),
  loo_csv_file,
  row.names = FALSE
)

## choose best by LOO automatically
selected_model <- rownames(loo_tab)[1]

selected_fit_file <- switch(
  selected_model,
  m1_species  = file.path(base_dir, "fits", paste0(fit_base_m1, ".rds")),
  m2_site_cov = file.path(base_dir, "fits", paste0(fit_base_m2, ".rds")),
  m3_cov      = file.path(base_dir, "fits", paste0(fit_base_m3, ".rds"))
)

alpha_structure_selection_info <- list(
  selected_structure_model = selected_model,
  selected_fit_file        = selected_fit_file,
  selected_growth_family   = family_tag,
  selected_error_model     = error_to_use,
  site_order               = c("BEN", "GH", "NAM", "PJ", "CS"),
  source_stage0_dir        = error_dir,
  source_stage0_best_model = selection_info$best_model,
  source_stage0_data_rds   = file.path(error_dir, "data", "analysis_data_family_error_compare.rds"),
  stage1_data_rds          = file.path(base_dir, "data", "analysis_data_stage1_with_germination.rds"),
  germination_priors_rds   = germination_priors_file,
  germination_site_priors_rds = file.path(base_dir, "data", "germination_site_priors_used.rds"),
  note = paste(
    "Germination priors are loaded and validated in stage 1 but are not used in the alpha-structure likelihood.",
    "Use them downstream in coexistence calculations via lambda_eff = g * lambda."
  )
)

## family-specific save
selection_rds_file <- file.path(
  base_dir,
  "selection",
  paste0("alpha_structure_selection_info_", family_tag, "_", error_tag, ".rds")
)

saveRDS(
  alpha_structure_selection_info,
  selection_rds_file
)

## also save a generic alias if this matches the primary stage-0 family
if (identical(family_tag, selection_info$best_growth_family)) {
  saveRDS(
    alpha_structure_selection_info,
    file.path(base_dir, "selection", "alpha_structure_selection_info.rds")
  )
}

############################
## 13. Save summaries
############################
summary_txt <- file.path(
  base_dir,
  "summaries",
  paste0(family_tag, "_", error_tag, "_structure_model_summaries.txt")
)

sink(summary_txt)
cat("\n====================\nFamily branch used\n====================\n")
print(family_tag)

cat("\n====================\nModel 1 summary\n====================\n")
print(summary(fit_m1))

cat("\n====================\nModel 2 summary\n====================\n")
print(summary(fit_m2))

cat("\n====================\nModel 3 summary\n====================\n")
print(summary(fit_m3))

cat("\n====================\nLOO summaries\n====================\n")
print(loo_m1)
print(loo_m2)
print(loo_m3)

cat("\n====================\nLOO comparison table\n====================\n")
print(loo_tab)

cat("\n====================\nSelected structure model\n====================\n")
print(alpha_structure_selection_info)
sink()

############################
## 14. Save session info
############################
capture.output(
  sessionInfo(),
  file = file.path(base_dir, "selection", paste0("sessionInfo_", family_tag, "_", error_tag, ".txt"))
)

############################
## 15. Final objects in console
############################
cat("\nDone.\n")
cat("Family branch used:\n", family_tag, "\n")
cat("Selected structure model:\n", selected_model, "\n")
cat("Selection object saved to:\n", normalizePath(selection_rds_file), "\n")

loo_tab
alpha_structure_selection_info
