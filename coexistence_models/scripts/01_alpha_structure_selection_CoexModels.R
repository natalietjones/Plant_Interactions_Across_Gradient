############################################################
## 01_alpha_structure_selection_with_germination.R
##
## Purpose:
## 1. Load selected family+error outputs from 00_family_error_selection/
## 2. Load and validate germination priors
## 3. Fit/load candidate alpha-structure models under the selected
##    log-Gaussian competition family branch
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
  "posterior",
  "here",
  "readr"
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
library(here)
library(readr)

options(mc.cores = max(1, parallel::detectCores() - 1))
options(brms.backend = "rstan")
# options(brms.backend = "cmdstanr")

############################
## 1. User options
############################
## Leave NULL to use the best family from stage 0.
## Set to "BH" or "Ricker" to force one branch.
family_override <- NULL

## Set TRUE only if you intentionally want to refit or recompute.
force_refit_models <- TRUE
force_recompute_loo <- TRUE

## After this first corrected rerun, you can set both switches back to FALSE.
## They are TRUE here because older cached fits/LOO may use (1 | BLK)
## rather than the corrected (1 | SiteBlock) structure.

############################
## 2. Set input/output folders
############################
## Robust to either:
##   A) R project root = parent folder containing coexistence_models/
##   B) R project root = coexistence_models/
##   C) R project root = folder containing 00_family_error_selection/ directly

project_root <- here::here()

if (basename(normalizePath(project_root, winslash = "/", mustWork = FALSE)) == "coexistence_models") {
  coex_dir <- project_root
} else if (dir.exists(file.path(project_root, "coexistence_models"))) {
  coex_dir <- file.path(project_root, "coexistence_models")
} else {
  coex_dir <- project_root
}

error_dir <- file.path(coex_dir, "00_family_error_selection")
base_dir  <- file.path(coex_dir, "01_alpha_structure_selection")

message("\nProject root: ", normalizePath(project_root, winslash = "/", mustWork = FALSE))
message("Coexistence model directory: ", normalizePath(coex_dir, winslash = "/", mustWork = FALSE))
message("Stage-0 directory: ", normalizePath(error_dir, winslash = "/", mustWork = FALSE))
message("Stage-1 directory: ", normalizePath(base_dir, winslash = "/", mustWork = FALSE))

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "fits"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "loo"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "selection"), showWarnings = FALSE, recursive = TRUE)

############################
## 3. Path helper functions
############################
resolve_existing_file <- function(candidates, label = "file") {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  found <- candidates[file.exists(candidates)]

  if (length(found) == 0) {
    stop(
      "Could not find ", label, ". Tried:\n  ",
      paste(normalizePath(candidates, winslash = "/", mustWork = FALSE), collapse = "\n  "),
      call. = FALSE
    )
  }

  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

basename_or_empty <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) {
    character(0)
  } else {
    basename(x)
  }
}

############################
## 4. Load selected family+error info and prepared data
############################
selection_info_file <- resolve_existing_file(
  c(
    file.path(error_dir, "selection", "family_error_model_selection_info.rds"),
    file.path(project_root, "00_family_error_selection", "selection", "family_error_model_selection_info.rds"),
    file.path(project_root, "coexistence_models", "00_family_error_selection", "selection", "family_error_model_selection_info.rds")
  ),
  label = "stage-0 family/error selection info"
)

selection_info <- readRDS(selection_info_file)

required_selection_fields <- c("best_growth_family", "best_error_model")
missing_selection_fields <- setdiff(required_selection_fields, names(selection_info))
if (length(missing_selection_fields) > 0) {
  stop(
    "The stage-0 selection object is missing required fields:\n  ",
    paste(missing_selection_fields, collapse = ", "),
    call. = FALSE
  )
}

stage0_data_from_selection <- if ("data_rds" %in% names(selection_info)) {
  selection_info$data_rds
} else if ("stage0_data_rds" %in% names(selection_info)) {
  selection_info$stage0_data_rds
} else {
  NA_character_
}

stage0_data_basename <- basename_or_empty(stage0_data_from_selection)

stage0_data_file <- resolve_existing_file(
  c(
    stage0_data_from_selection,
    if (length(stage0_data_basename) > 0) file.path(error_dir, "data", stage0_data_basename),
    file.path(error_dir, "data", "analysis_data_family_error_compare.rds"),
    file.path(project_root, "00_family_error_selection", "data", "analysis_data_family_error_compare.rds"),
    file.path(project_root, "coexistence_models", "00_family_error_selection", "data", "analysis_data_family_error_compare.rds")
  ),
  label = "stage-0 analysis data"
)

dat_comp <- readRDS(stage0_data_file)

cat("\nLoaded family+error selection info from:\n")
cat(selection_info_file, "\n")
print(selection_info)

cat("\nLoaded stage-0 analysis data from:\n")
cat(stage0_data_file, "\n")

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
    ),
    call. = FALSE
  )
}

if (!family_to_use %in% c("BH", "Ricker")) {
  stop("family_to_use must be either 'BH' or 'Ricker'.", call. = FALSE)
}

############################
## 5. Validate and standardise stage-0 data
############################
if (!"fitness_log" %in% names(dat_comp) && "fitness" %in% names(dat_comp)) {
  message("fitness_log was not present; creating it as log(fitness + 1).")
  dat_comp <- dat_comp %>%
    mutate(fitness_log = log(fitness + 1))
}

required_dat_cols <- c(
  "fitness_log", "focsp", "Site", "Cov", "BLK",
  "intra_sc", "inter_sc", "fac_sc"
)

missing_dat_cols <- setdiff(required_dat_cols, names(dat_comp))
if (length(missing_dat_cols) > 0) {
  stop(
    "The stage-0 analysis data is missing required columns:\n  ",
    paste(missing_dat_cols, collapse = ", "),
    call. = FALSE
  )
}

factor_with_preferred_levels <- function(x, preferred_levels) {
  x_chr <- as.character(x)
  obs <- unique(x_chr[!is.na(x_chr)])

  if (all(obs %in% preferred_levels)) {
    factor(x_chr, levels = preferred_levels[preferred_levels %in% obs])
  } else {
    factor(x_chr)
  }
}

## Match the stage-0 coding exactly:
##   Site: BEN, GH, NAM, PJ, CS
##   Cov:  Shade, Sun
##   focsp: TRCY, TROR
## Stage-0 data should already be standardised, but these recodes make the
## script safer if an older/intermediate object is loaded.
standardise_site <- function(x) {
  x_chr <- stringr::str_trim(as.character(x))
  x_up  <- stringr::str_to_upper(x_chr)

  dplyr::case_when(
    x_up %in% c("BEN", "BENDERING") ~ "BEN",
    x_up %in% c("GH")               ~ "GH",
    x_up %in% c("NAM")              ~ "NAM",
    x_up %in% c("PJ")               ~ "PJ",
    x_up %in% c("CS")               ~ "CS",
    TRUE                            ~ x_chr
  )
}

standardise_cov <- function(x) {
  x_chr <- stringr::str_trim(as.character(x))
  x_up  <- stringr::str_to_upper(x_chr)

  dplyr::case_when(
    x_up %in% c("SH", "SHADE") ~ "Shade",
    x_up %in% c("SUN")         ~ "Sun",
    TRUE                       ~ x_chr
  )
}

standardise_focsp <- function(x) {
  stringr::str_to_upper(stringr::str_trim(as.character(x)))
}

site_order_default  <- c("BEN", "GH", "NAM", "PJ", "CS")
focsp_order_default <- c("TRCY", "TROR")
cov_order_default   <- c("Shade", "Sun")

dat_comp <- dat_comp %>%
  mutate(
    Site  = factor_with_preferred_levels(standardise_site(Site), site_order_default),
    focsp = factor_with_preferred_levels(standardise_focsp(focsp), focsp_order_default),
    Cov   = factor_with_preferred_levels(standardise_cov(Cov), cov_order_default),
    BLK   = factor(BLK),
    SiteBlock = interaction(Site, BLK, drop = TRUE, sep = "_")
  )

if (!"SiteBlock" %in% names(dat_comp)) {
  stop("SiteBlock was not created correctly.", call. = FALSE)
}

cat("\nRandom-effect grouping check: SiteBlock levels by Site\n")
print(dat_comp %>% distinct(Site, BLK, SiteBlock) %>% count(Site, name = "n_site_blocks"))

if (any(is.na(dat_comp$fitness_log))) {
  stop("fitness_log contains NA values. Check fitness and the log transformation.", call. = FALSE)
}

if (any(dat_comp$intra_sc < 0, na.rm = TRUE) || any(dat_comp$inter_sc < 0, na.rm = TRUE)) {
  warning("Some scaled competitor densities are negative. This is only appropriate if the model was intentionally centred/scaled.")
}

if (any(dat_comp$fac_sc < 0, na.rm = TRUE)) {
  warning(
    "fac_sc contains negative values. If fac_sc is facilitator density, this means it has probably been centred.\n",
    "That is fine for fitted effects, but fac_sc = 0 will not mean 'no facilitator'."
  )
}

site_levels  <- levels(dat_comp$Site)
focsp_levels <- levels(dat_comp$focsp)

if (length(site_levels) == 0 || length(focsp_levels) == 0) {
  stop("Site or focsp levels were not set correctly after factor conversion.", call. = FALSE)
}

############################
## 6. Load and validate germination priors
############################
germination_priors_from_selection <- if ("germination_priors_rds" %in% names(selection_info)) {
  selection_info$germination_priors_rds
} else {
  NA_character_
}

germination_priors_basename <- basename_or_empty(germination_priors_from_selection)

germination_priors_file <- resolve_existing_file(
  c(
    germination_priors_from_selection,
    if (length(germination_priors_basename) > 0) file.path(coex_dir, "00_germination_priors", "data", germination_priors_basename),
    if (length(germination_priors_basename) > 0) file.path(coex_dir, "00_germination_priors", germination_priors_basename),
    if (length(germination_priors_basename) > 0) file.path(error_dir, "data", germination_priors_basename),
    file.path(coex_dir, "00_germination_priors", "data", "germination_priors.rds"),
    file.path(coex_dir, "00_germination_priors", "germination_priors.rds"),
    file.path(coex_dir, "00_germination_priors", "data", "germination_site_priors.rds"),
    file.path(error_dir, "data", "germination_priors.rds"),
    file.path(error_dir, "data", "germination_site_priors.rds")
  ),
  label = "germination priors RDS"
)

germ_info <- readRDS(germination_priors_file)

if (!"site_priors" %in% names(germ_info)) {
  stop("The germination priors RDS does not contain an element named 'site_priors'.", call. = FALSE)
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
    paste(missing_germ_cols, collapse = ", "),
    call. = FALSE
  )
}

expected_site_sp <- tidyr::expand_grid(
  Site  = site_levels,
  focsp = focsp_levels
)

germ_check <- expected_site_sp %>%
  left_join(germ_priors, by = c("Site", "focsp"))

if (any(is.na(germ_check$g_alpha_final) | is.na(germ_check$g_beta_final))) {
  missing_germ <- germ_check %>%
    filter(is.na(g_alpha_final) | is.na(g_beta_final)) %>%
    select(Site, focsp)

  print(missing_germ)

  stop(
    "Missing germination priors for at least one Site x focsp combination used in the stage-0 data.\n",
    "Check whether Site/focsp labels match exactly between the model data and 00_germination_priors output.",
    call. = FALSE
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

## Join germination means to the data for bookkeeping only.
## These columns are not used in the model formulas at this stage.
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
    focsp = factor(focsp, levels = focsp_levels),
    Cov   = factor(standardise_cov(Cov), levels = levels(dat_comp$Cov)),
    BLK   = factor(BLK),
    SiteBlock = interaction(Site, BLK, drop = TRUE, sep = "_")
  )

if (any(is.na(dat_stage1$g_mean_final))) {
  stop("Germination join failed for some rows in dat_stage1.", call. = FALSE)
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
## 7. Quick checks
############################
cat("\nSummary of log(fitness + 1):\n")
print(summary(dat_stage1$fitness_log))

cat("\nLevels:\n")
print(levels(dat_stage1$focsp))
print(levels(dat_stage1$Site))
print(levels(dat_stage1$Cov))

cat("\nGermination priors loaded from:\n")
cat(germination_priors_file, "\n")

cat("\nGermination priors used:\n")
print(
  germ_check %>%
    select(Site, focsp, g_mean_final) %>%
    arrange(focsp, Site)
)

capture.output(
  {
    cat("\nProject root:\n")
    print(project_root)
    cat("\nCoexistence model directory:\n")
    print(coex_dir)
    cat("\nStage-0 selection file:\n")
    print(selection_info_file)
    cat("\nStage-0 data file:\n")
    print(stage0_data_file)
    cat("\nGermination priors file:\n")
    print(germination_priors_file)
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
## 8. Family, priors, settings
############################
model_family <- gaussian()

## Match the log-Gaussian priors used in stage 0.
## etalam is on log(lambda); etaai and etaaj are on log(alpha).
priors_common <- c(
  prior(normal(1.5, 1), class = "b", nlpar = "etalam"),
  prior(normal(-3, 1),  class = "b", nlpar = "etaai"),
  prior(normal(-3, 1),  class = "b", nlpar = "etaaj"),
  prior(exponential(1), class = "sd", nlpar = "etalam"),
  prior(exponential(1), class = "sigma")
)

ctrl <- list(adapt_delta = 0.995, max_treedepth = 14)

############################
## 9. Model formulas
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
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m2 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp + focsp:Site + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Site + focsp:Cov,
    nl = TRUE
  )

  bf_m3 <- bf(
    fitness_log ~ log(
      exp(etalam) / (1 + exp(etaai) * intra_sc + exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
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
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp,
    etaaj ~ 0 + focsp,
    nl = TRUE
  )

  bf_m2 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp + focsp:Site + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Site + focsp:Cov,
    nl = TRUE
  )

  bf_m3 <- bf(
    fitness_log ~ log(
      exp(etalam - exp(etaai) * intra_sc - exp(etaaj) * inter_sc) + 1
    ),
    etalam ~ 0 + focsp + focsp:Site + focsp:Cov + focsp:fac_sc +
      focsp:Site:fac_sc + focsp:Cov:fac_sc + (1 | SiteBlock),
    etaai ~ 0 + focsp + focsp:Cov,
    etaaj ~ 0 + focsp + focsp:Cov,
    nl = TRUE
  )
}

############################
## 10. Helper functions
############################
fit_or_load_brms <- function(formula_obj, data, family, prior, seed, file_base, control, base_dir,
                             force_refit = FALSE) {
  rds_file <- file.path(base_dir, "fits", paste0(file_base, ".rds"))
  fit_file <- file.path(base_dir, "fits", file_base)

  if (file.exists(rds_file) && !isTRUE(force_refit)) {
    cat("\nLoading existing model:", rds_file, "\n")
    fit <- readRDS(rds_file)

    if (!inherits(fit, "brmsfit")) {
      stop("Loaded object is not a brmsfit: ", rds_file, call. = FALSE)
    }

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
      cores = min(4, max(1, parallel::detectCores() - 1)),
      refresh = 50,
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
    png_file <- file.path(base_dir, "figures", paste0(file_base, "_ppcheck.png"))
    png(png_file, width = 1400, height = 1000, res = 160)
    on.exit({
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    }, add = TRUE)
    print(pp_check(fit))
    grDevices::dev.off()
  }, silent = TRUE)

  invisible(fit)
}

loo_or_load <- function(fit, loo_file, use_reloo = FALSE, base_dir, force_recompute = FALSE) {
  loo_path <- file.path(base_dir, "loo", loo_file)

  if (file.exists(loo_path) && !isTRUE(force_recompute)) {
    cat("\nLoading existing LOO object:", loo_path, "\n")
    loo_obj <- readRDS(loo_path)

  } else {
    cat("\nComputing LOO and saving to:", loo_path, "\n")

    if (use_reloo) {
      loo_obj <- loo(
        fit,
        reloo = TRUE,
        cores = getOption("mc.cores", 1)
      )
    } else {
      loo_obj <- loo(
        fit,
        moment_match = TRUE,
        recompile = TRUE,
        cores = getOption("mc.cores", 1)
      )
    }

    saveRDS(loo_obj, loo_path)
  }

  loo_obj
}

############################
## 11. File tags
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
## 12. Fit or load models
############################
fit_m1 <- fit_or_load_brms(
  formula_obj  = bf_m1,
  data         = dat_stage1,
  family       = model_family,
  prior        = priors_common,
  seed         = 1234,
  file_base    = fit_base_m1,
  control      = ctrl,
  base_dir     = base_dir,
  force_refit  = force_refit_models
)

fit_m2 <- fit_or_load_brms(
  formula_obj  = bf_m2,
  data         = dat_stage1,
  family       = model_family,
  prior        = priors_common,
  seed         = 1235,
  file_base    = fit_base_m2,
  control      = ctrl,
  base_dir     = base_dir,
  force_refit  = force_refit_models
)

fit_m3 <- fit_or_load_brms(
  formula_obj  = bf_m3,
  data         = dat_stage1,
  family       = model_family,
  prior        = priors_common,
  seed         = 1236,
  file_base    = fit_base_m3,
  control      = ctrl,
  base_dir     = base_dir,
  force_refit  = force_refit_models
)

############################
## 13. LOO: compute or load
############################
loo_m1 <- loo_or_load(
  fit             = fit_m1,
  loo_file        = loo_file_m1,
  use_reloo       = FALSE,
  base_dir        = base_dir,
  force_recompute = force_recompute_loo
)

loo_m2 <- loo_or_load(
  fit             = fit_m2,
  loo_file        = loo_file_m2,
  use_reloo       = FALSE,
  base_dir        = base_dir,
  force_recompute = force_recompute_loo
)

loo_m3 <- loo_or_load(
  fit             = fit_m3,
  loo_file        = loo_file_m3,
  use_reloo       = FALSE,
  base_dir        = base_dir,
  force_recompute = force_recompute_loo
)

cat("\n====================\nLOO summaries\n====================\n")
print(loo_m1)
print(loo_m2)
print(loo_m3)

############################
## 14. LOO comparison and final model choice
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

## Choose best by LOO automatically.
## Check elpd_diff and se_diff manually before treating the top model
## as meaningfully better.
selected_model <- rownames(loo_tab)[1]

selected_fit_file <- switch(
  selected_model,
  m1_species  = file.path(base_dir, "fits", paste0(fit_base_m1, ".rds")),
  m2_site_cov = file.path(base_dir, "fits", paste0(fit_base_m2, ".rds")),
  m3_cov      = file.path(base_dir, "fits", paste0(fit_base_m3, ".rds"))
)

alpha_structure_selection_info <- list(
  selected_structure_model      = selected_model,
  selected_fit_file             = selected_fit_file,
  selected_growth_family        = family_tag,
  selected_error_model          = error_to_use,
  site_order                    = site_levels,
  source_stage0_dir             = error_dir,
  source_stage0_selection_rds   = selection_info_file,
  source_stage0_best_model      = if ("best_model" %in% names(selection_info)) selection_info$best_model else NA_character_,
  source_stage0_data_rds        = stage0_data_file,
  stage1_data_rds               = file.path(base_dir, "data", "analysis_data_stage1_with_germination.rds"),
  germination_priors_rds        = germination_priors_file,
  germination_site_priors_rds   = file.path(base_dir, "data", "germination_site_priors_used.rds"),
  note = paste(
    "Germination priors are loaded and validated in stage 1 but are not used in the alpha-structure likelihood.",
    "Use them downstream in coexistence calculations via lambda_eff = g * lambda."
  )
)

## Family-specific save
selection_rds_file <- file.path(
  base_dir,
  "selection",
  paste0("alpha_structure_selection_info_", family_tag, "_", error_tag, ".rds")
)

saveRDS(
  alpha_structure_selection_info,
  selection_rds_file
)

## Also save a generic alias if this matches the primary stage-0 family.
if (identical(family_tag, selection_info$best_growth_family)) {
  saveRDS(
    alpha_structure_selection_info,
    file.path(base_dir, "selection", "alpha_structure_selection_info.rds")
  )
}

############################
## 15. Save summaries
############################
summary_txt <- file.path(
  base_dir,
  "summaries",
  paste0(family_tag, "_", error_tag, "_structure_model_summaries.txt")
)

sink(summary_txt)
cat("\n====================\nInput paths\n====================\n")
cat("project_root:", project_root, "\n")
cat("coex_dir:", coex_dir, "\n")
cat("selection_info_file:", selection_info_file, "\n")
cat("stage0_data_file:", stage0_data_file, "\n")
cat("germination_priors_file:", germination_priors_file, "\n")

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
## 16. Save session info
############################
capture.output(
  sessionInfo(),
  file = file.path(base_dir, "selection", paste0("sessionInfo_", family_tag, "_", error_tag, ".txt"))
)

############################
## 17. Final objects in console
############################
cat("\nDone.\n")
cat("Family branch used:\n", family_tag, "\n")
cat("Selected structure model:\n", selected_model, "\n")
cat("Selection object saved to:\n", normalizePath(selection_rds_file, winslash = "/", mustWork = FALSE), "\n")

loo_tab
alpha_structure_selection_info
