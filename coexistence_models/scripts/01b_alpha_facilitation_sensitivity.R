############################################################
## 01b_alpha_facilitation_sensitivity_checked.R
##
## Purpose:
## Compare the selected baseline alpha-structure model (m1_species)
## against targeted sensitivity models that allow facilitator density
## (fac_sc) to modify alpha_ij only, or both alpha_ii and alpha_ij.
##
## Key corrections in this version:
## 1. Uses robust project paths.
## 2. Resolves stale saved file paths from previous runs.
## 3. Standardises Site, Cov, focsp, BLK and SiteBlock labels.
## 4. Uses SiteBlock as the random intercept, not BLK.
## 5. Uses priors consistent with the corrected model-selection scripts.
## 6. Forces refits/recomputed LOO by default so cached BLK models are not reused.
############################################################

rm(list = ls())

############################
## 0. Install/load packages
############################
required_packages <- c(
  "tidyverse",
  "brms",
  "loo",
  "posterior",
  "here"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

options(mc.cores = max(1, parallel::detectCores() - 1))
options(brms.backend = "rstan")
# options(brms.backend = "cmdstanr")

############################
## 1. User options
############################

## Leave NULL to use the previously selected family from alpha_structure_selection_info.rds.
family_override <- NULL

## Important: keep these TRUE for the first corrected rerun.
## Otherwise the script may silently reuse old cached fits/LOO objects
force_refit_baseline_m1 <- FALSE
force_refit_models      <- FALSE
force_recompute_loo     <- FALSE

############################
## 2. Paths
############################

## Robustly handle cases where the R project root is:
##   - the folder containing coexistence_models/
##   - coexistence_models/ itself
##   - the folder containing 01_alpha_structure_selection/ directly
if (dir.exists(here::here("coexistence_models", "01_alpha_structure_selection"))) {
  project_dir <- here::here("coexistence_models")
} else if (dir.exists(here::here("01_alpha_structure_selection"))) {
  project_dir <- here::here()
} else {
  stop(
    "Could not find 01_alpha_structure_selection. Check your R project root.\n",
    "Current here::here() is: ", here::here()
  )
}

base_stage1_dir <- file.path(project_dir, "01_alpha_structure_selection")
base_dir        <- file.path(project_dir, "01b_alpha_facilitation_sensitivity")

fits_dir      <- file.path(base_dir, "fits")
loo_dir       <- file.path(base_dir, "loo")
summaries_dir <- file.path(base_dir, "summaries")
figures_dir   <- file.path(base_dir, "figures")
selection_dir <- file.path(base_dir, "selection")

dir.create(base_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(fits_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(loo_dir,       showWarnings = FALSE, recursive = TRUE)
dir.create(summaries_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(selection_dir, showWarnings = FALSE, recursive = TRUE)

resolve_existing_file <- function(path, fallback_dirs = character(), label = "file") {
  valid_path <- !is.null(path) && length(path) == 1 && !is.na(path) && nzchar(as.character(path))

  candidates <- character()
  if (valid_path) {
    path <- as.character(path)
    candidates <- c(candidates, path)

    if (length(fallback_dirs) > 0) {
      candidates <- c(candidates, file.path(fallback_dirs, basename(path)))
    }
  }

  candidates <- unique(candidates)
  existing <- candidates[file.exists(candidates)]

  if (length(existing) > 0) {
    return(normalizePath(existing[1], winslash = "/", mustWork = TRUE))
  }

  stop(
    "Could not find ", label, ".\n",
    "Original path: ", ifelse(valid_path, path, "<NULL/empty>"), "\n",
    "Tried:\n", paste(candidates, collapse = "\n")
  )
}

cat("\nProject dir:      ", normalizePath(project_dir), "\n", sep = "")
cat("Stage-1 dir:      ", normalizePath(base_stage1_dir), "\n", sep = "")
cat("Sensitivity dir:  ", normalizePath(base_dir), "\n\n", sep = "")

############################
## 3. Load previous selection info and data
############################

selection_info_file <- file.path(base_stage1_dir, "selection", "alpha_structure_selection_info.rds")
if (!file.exists(selection_info_file)) {
  stop("Cannot find alpha_structure_selection_info.rds at: ", selection_info_file)
}

selection_info <- readRDS(selection_info_file)

stage1_data_file <- resolve_existing_file(
  selection_info$stage1_data_rds,
  fallback_dirs = c(
    file.path(base_stage1_dir, "data"),
    file.path(base_stage1_dir, "draws"),
    base_stage1_dir,
    project_dir
  ),
  label = "stage-1 data RDS"
)

dat_stage1 <- readRDS(stage1_data_file)

required_cols <- c("fitness_log", "focsp", "Site", "Cov", "BLK", "fac_sc", "intra_sc", "inter_sc")
missing_cols <- setdiff(required_cols, names(dat_stage1))
if (length(missing_cols) > 0) {
  stop("dat_stage1 is missing required columns: ", paste(missing_cols, collapse = ", "))
}

## Standardise labels to match the corrected structure-selection scripts.
dat_stage1 <- dat_stage1 %>%
  dplyr::mutate(
    Site  = toupper(as.character(Site)),
    focsp = toupper(as.character(focsp)),
    Cov = dplyr::case_when(
      as.character(Cov) %in% c("SH", "SHADE", "Shade", "shade") ~ "Shade",
      as.character(Cov) %in% c("SUN", "Sun", "sun") ~ "Sun",
      TRUE ~ as.character(Cov)
    ),
    BLK = factor(BLK),
    SiteBlock = interaction(Site, as.character(BLK), drop = TRUE, sep = "_")
  )

allowed_sites <- c("BEN", "GH", "NAM", "PJ", "CS")
allowed_cov   <- c("Shade", "Sun")
allowed_focsp <- c("TRCY", "TROR")

bad_sites <- setdiff(unique(dat_stage1$Site), allowed_sites)
bad_cov   <- setdiff(unique(dat_stage1$Cov), allowed_cov)
bad_focsp <- setdiff(unique(dat_stage1$focsp), allowed_focsp)

if (length(bad_sites) > 0) stop("Unexpected Site values after standardisation: ", paste(bad_sites, collapse = ", "))
if (length(bad_cov) > 0)   stop("Unexpected Cov values after standardisation: ", paste(bad_cov, collapse = ", "))
if (length(bad_focsp) > 0) stop("Unexpected focsp values after standardisation: ", paste(bad_focsp, collapse = ", "))

site_order <- allowed_sites[allowed_sites %in% unique(dat_stage1$Site)]
cov_order  <- allowed_cov[allowed_cov %in% unique(dat_stage1$Cov)]

dat_stage1 <- dat_stage1 %>%
  dplyr::mutate(
    Site      = factor(Site, levels = site_order),
    Cov       = factor(Cov, levels = cov_order),
    focsp     = factor(focsp, levels = allowed_focsp),
    SiteBlock = factor(SiteBlock)
  )

if (!any(abs(dat_stage1$fac_sc) < .Machine$double.eps^0.5, na.rm = TRUE)) {
  warning(
    "No fac_sc values equal to 0 were found. This script assumes fac_sc = 0 means no facilitator. ",
    "If fac_sc was centred or z-scaled, the no-facilitator interpretation is wrong."
  )
}

family_to_use <- if (is.null(family_override)) {
  as.character(selection_info$selected_growth_family)
} else {
  as.character(family_override)
}

if (!family_to_use %in% c("BH", "Ricker")) {
  stop("family_to_use must be either 'BH' or 'Ricker'. Current value: ", family_to_use)
}

if (!identical(selection_info$selected_structure_model, "m1_species")) {
  warning(
    "The previously selected model was not m1_species.\n",
    "This script will still run, but its main purpose is to compare m1_species ",
    "to facilitator-on-alpha sensitivity models."
  )
}

cat("\nLoaded selection info:\n")
print(selection_info)

cat("\nFamily used in sensitivity analysis:\n")
print(family_to_use)

cat("\nStage-1 data file:\n")
print(stage1_data_file)

cat("\nRandom effect used in this script: (1 | SiteBlock)\n")

############################
## 4. Family, priors, settings
############################

model_family <- gaussian()

## Match corrected model-selection scripts.
priors_common <- c(
  prior(normal(1.5, 1), class = "b", nlpar = "etalam"),
  prior(normal(-3, 1),  class = "b", nlpar = "etaai"),
  prior(normal(-3, 1),  class = "b", nlpar = "etaaj"),
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
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m4 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )

  bf_m5 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
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
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m4 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )

  bf_m5 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp + focsp:fac_sc,
    etaaj ~ 0 + focsp + focsp:fac_sc,
    nl = TRUE
  )
}

############################
## 6. Helper functions
############################

fit_or_load_brms <- function(formula_obj, data, family, prior, seed, file_base,
                             control, fits_dir, summaries_dir, figures_dir,
                             force_refit = FALSE) {
  rds_file <- file.path(fits_dir, paste0(file_base, ".rds"))
  fit_file <- file.path(fits_dir, file_base)

  if (file.exists(rds_file) && !force_refit) {
    cat("\nLoading existing model:", rds_file, "\n")
    fit <- readRDS(rds_file)
  } else {
    if (file.exists(rds_file) && force_refit) {
      cat("\nForce refit is TRUE; ignoring existing model:", rds_file, "\n")
    }

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
      file_refit = "always"
    )
    saveRDS(fit, rds_file)
  }

  capture.output(
    print(summary(fit), digits = 2),
    file = file.path(summaries_dir, paste0(file_base, "_summary.txt"))
  )

  try({
    png(
      file.path(figures_dir, paste0(file_base, "_ppcheck.png")),
      width = 1400, height = 1000, res = 160
    )
    print(pp_check(fit))
    dev.off()
  }, silent = TRUE)

  invisible(fit)
}

loo_or_load <- function(fit, loo_file, loo_dir, force_recompute = FALSE) {
  loo_path <- file.path(loo_dir, loo_file)

  if (file.exists(loo_path) && !force_recompute) {
    cat("\nLoading existing LOO object:", loo_path, "\n")
    loo_obj <- readRDS(loo_path)
  } else {
    if (file.exists(loo_path) && force_recompute) {
      cat("\nForce recompute is TRUE; ignoring existing LOO object:", loo_path, "\n")
    }

    cat("\nComputing LOO and saving to:", loo_path, "\n")
    loo_obj <- loo(fit, moment_match = TRUE)
    saveRDS(loo_obj, loo_path)
  }

  loo_obj
}

############################
## 7. Fit/load baseline m1
############################
## Do not compare m4/m5 against an old m1 cached with (1 | BLK).
## For the first corrected rerun, force_refit_baseline_m1 should be TRUE.

fit_base_m1 <- paste0("fit_", family_to_use, "_m1_baseline_loggaussian")

if (!force_refit_baseline_m1) {
  fit_m1_path <- resolve_existing_file(
    selection_info$selected_fit_file,
    fallback_dirs = c(file.path(base_stage1_dir, "fits"), fits_dir, base_stage1_dir, project_dir),
    label = "selected baseline m1 fit"
  )

  cat("\nLoading baseline selected m1 model from:\n", fit_m1_path, "\n")
  fit_m1 <- readRDS(fit_m1_path)
  baseline_fit_file_out <- fit_m1_path
} else {
  cat("\nForce refit baseline m1 is TRUE; refitting m1 with the corrected SiteBlock structure.\n")
  fit_m1 <- fit_or_load_brms(
    formula_obj   = bf_m1,
    data          = dat_stage1,
    family        = model_family,
    prior         = priors_common,
    seed          = 2001,
    file_base     = fit_base_m1,
    control       = ctrl,
    fits_dir      = fits_dir,
    summaries_dir = summaries_dir,
    figures_dir   = figures_dir,
    force_refit   = force_refit_models
  )
  baseline_fit_file_out <- file.path(fits_dir, paste0(fit_base_m1, ".rds"))
}

capture.output(
  print(summary(fit_m1), digits = 2),
  file = file.path(summaries_dir, paste0(fit_base_m1, "_summary.txt"))
)

############################
## 8. Fit/load new facilitator-on-alpha models
############################

fit_base_m4 <- paste0("fit_", family_to_use, "_m4_fac_inter_loggaussian")
fit_base_m5 <- paste0("fit_", family_to_use, "_m5_fac_both_loggaussian")

fit_m4 <- fit_or_load_brms(
  formula_obj   = bf_m4,
  data          = dat_stage1,
  family        = model_family,
  prior         = priors_common,
  seed          = 2004,
  file_base     = fit_base_m4,
  control       = ctrl,
  fits_dir      = fits_dir,
  summaries_dir = summaries_dir,
  figures_dir   = figures_dir,
  force_refit   = force_refit_models
)

fit_m5 <- fit_or_load_brms(
  formula_obj   = bf_m5,
  data          = dat_stage1,
  family        = model_family,
  prior         = priors_common,
  seed          = 2005,
  file_base     = fit_base_m5,
  control       = ctrl,
  fits_dir      = fits_dir,
  summaries_dir = summaries_dir,
  figures_dir   = figures_dir,
  force_refit   = force_refit_models
)

############################
## 9. LOO comparison
############################

loo_m1 <- loo_or_load(
  fit = fit_m1,
  loo_file = paste0("loo_", family_to_use, "_m1_baseline_loggaussian.rds"),
  loo_dir = loo_dir,
  force_recompute = force_recompute_loo
)

loo_m4 <- loo_or_load(
  fit = fit_m4,
  loo_file = paste0("loo_", family_to_use, "_m4_fac_inter_loggaussian.rds"),
  loo_dir = loo_dir,
  force_recompute = force_recompute_loo
)

loo_m5 <- loo_or_load(
  fit = fit_m5,
  loo_file = paste0("loo_", family_to_use, "_m5_fac_both_loggaussian.rds"),
  loo_dir = loo_dir,
  force_recompute = force_recompute_loo
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
  file.path(selection_dir, paste0(family_to_use, "_facilitation_alpha_sensitivity_LOO.csv")),
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
  stage1_data_file = stage1_data_file,
  baseline_fit_file = baseline_fit_file_out,
  m4_fit_file = file.path(fits_dir, paste0(fit_base_m4, ".rds")),
  m5_fit_file = file.path(fits_dir, paste0(fit_base_m5, ".rds")),
  random_effect = "(1 | SiteBlock)",
  priors = "etalam: normal(1.5, 1); etaai/etaaj: normal(-3, 1); sd/sigma: exponential(1)",
  note = paste(
    "m1 tests facilitator effects through lambda only;",
    "m4 tests facilitator effects through lambda plus alpha_ij;",
    "m5 tests facilitator effects through lambda plus alpha_ii and alpha_ij."
  )
)

saveRDS(
  selection_out,
  file.path(selection_dir, paste0("facilitation_alpha_sensitivity_selection_", family_to_use, ".rds"))
)

############################
## 10. Save text summary
############################

summary_txt <- file.path(
  summaries_dir,
  paste0(family_to_use, "_facilitation_alpha_sensitivity_summary.txt")
)

sink(summary_txt)
cat("\n====================\nFamily used\n====================\n")
print(family_to_use)

cat("\n====================\nStage-1 data file\n====================\n")
print(stage1_data_file)

cat("\n====================\nRandom effect\n====================\n")
print("(1 | SiteBlock)")

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
