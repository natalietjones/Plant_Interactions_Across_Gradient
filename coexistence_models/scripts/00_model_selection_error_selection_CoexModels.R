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
##
## Key modelling decisions:
## - Conspecific and heterospecific densities are raw counts.
##   This preserves alpha as competition strength per competitor individual.
## - Facilitator density is raw NO_GRASS.
##   It is not centred and not z-scaled, so NO_GRASS = 0 means no facilitator.
## - Raw facilitator effects enter etalam = log(lambda), so the prior on
##   facilitator slopes must be scaled to the raw NO_GRASS unit.
## - Germination priors are loaded and validated here but are not used
##   in this stage-0 likelihood.
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
## 1. Set output folders
############################

base_dir <- here::here("coexistence_models", "00_family_error_selection")

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

data_file <- here::here("data", "field", "fullWAdata.csv")

stopifnot(file.exists(data_file))

dat0 <- readr::read_csv(
  data_file,
  show_col_types = FALSE
)

site_order <- c("BEN", "GH", "NAM", "PJ", "CS")
cov_order  <- c("Shade", "Sun")
sp_order   <- c("TRCY", "TROR")

dat_comp <- dat0 %>%
  mutate(
    Site = case_when(
      Site %in% c("Ben", "BEN", "Bendering") ~ "BEN",
      Site %in% c("Nam", "NAM")              ~ "NAM",
      TRUE                                   ~ as.character(Site)
    ),
    Cov = case_when(
      Cov %in% c("SH", "Shade", "shade") ~ "Shade",
      Cov %in% c("SUN", "Sun", "sun")    ~ "Sun",
      TRUE                               ~ as.character(Cov)
    ),
    Site     = factor(Site, levels = site_order),
    Cov      = factor(Cov, levels = cov_order),
    BLK      = factor(BLK),
    focsp    = factor(focsp, levels = sp_order),
    dens_trt = factor(dens_trt),
    germ     = stringr::str_to_upper(as.character(FOCAL_GERM)) %in% c("Y", "1", "TRUE")
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
    ## Raw competitor densities.
    ## These are not centred and not z-scaled.
    ## Alpha is therefore interpretable as the effect per competitor individual.
    intra_n = intra,
    inter_n = inter,
    
    ## Raw facilitator density.
    ## This is not centred and not z-scaled, so 0 = no facilitator.
    fac_raw = fac,
    
    ## Legacy aliases retained so older downstream scripts do not break immediately.
    ## These are raw values, not scaled values.
    intra_sc = intra_n,
    inter_sc = inter_n,
    fac_sc   = fac_raw,
    
    fitness_raw   = fitness,
    fitness_log   = log(fitness + 1),
    fitness_count = pmax(0, round(fitness))
  ) %>%
  droplevels() %>%
  mutate(
    SiteBlock = interaction(Site, BLK, drop = TRUE)
  )

## Hard check for count models.
## If this fails, do not fit Poisson or negative-binomial models to rounded values
## unless you explicitly justify that decision.
fitness_is_integer <- all(abs(dat_comp$fitness_raw - round(dat_comp$fitness_raw)) < 1e-8)

if (!fitness_is_integer) {
  stop(
    "fitness is not integer-valued. Do not fit Poisson or negative-binomial models ",
    "to rounded fitness_count unless this rounding is explicitly justified."
  )
}

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

germination_priors_file <- here::here(
  "coexistence_models",
  "00_germination_priors",
  "selection",
  "germination_priors.rds"
)

if (!file.exists(germination_priors_file)) {
  stop(
    "Could not find germination prior file:\n  ",
    germination_priors_file,
    "\nRun coexistence_models/00_germination_priors.R first."
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

required_germ_cols <- c(
  "Site", "focsp",
  "g_alpha_final", "g_beta_final", "g_mean_final"
)

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
  print(germ_check)
  
  stop(
    "Missing germination priors for at least one Site x focsp combination used in fullWAdata.csv.\n",
    "Check site/focsp naming in both fullWAdata.csv and germination_priors.rds."
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
print(fitness_is_integer)

cat("\nSummary of rounded count fitness:\n")
print(summary(dat_comp$fitness_count))

cat("\nSummary of raw facilitator density:\n")
fac_summary <- dat_comp %>%
  summarise(
    min_fac    = min(fac_raw, na.rm = TRUE),
    q25_fac    = as.numeric(quantile(fac_raw, 0.25, na.rm = TRUE)),
    median_fac = median(fac_raw, na.rm = TRUE),
    q75_fac    = as.numeric(quantile(fac_raw, 0.75, na.rm = TRUE)),
    q90_fac    = as.numeric(quantile(fac_raw, 0.90, na.rm = TRUE)),
    max_fac    = max(fac_raw, na.rm = TRUE)
  )

print(fac_summary)

write.csv(
  fac_summary,
  file.path(base_dir, "summaries", "raw_facilitator_density_summary.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("\nSummary of raw fitness:\n")
    print(summary(dat_comp$fitness_raw))
    
    cat("\nCheck whether fitness is integer-valued:\n")
    print(fitness_is_integer)
    
    cat("\nSummary of rounded count fitness:\n")
    print(summary(dat_comp$fitness_count))
    
    cat("\nSummary of raw facilitator density:\n")
    print(fac_summary)
    
    cat("\nCounts by Site, Cov, focsp:\n")
    print(dat_comp %>% count(Site, Cov, focsp))
    
    cat("\nCounts by SiteBlock:\n")
    print(dat_comp %>% count(SiteBlock, sort = TRUE))
    
    cat("\nGermination priors loaded from:\n")
    print(normalizePath(germination_priors_file))
  },
  file = file.path(base_dir, "summaries", "data_checks.txt")
)

############################
## 5. Model formulas
############################

## Beverton-Holt:
##   F = lambda / (1 + alpha_ii * N_i + alpha_ij * N_j)
##
## Ricker:
##   F = lambda * exp(- alpha_ii * N_i - alpha_ij * N_j)
##
## Here:
## - intra_n and inter_n are raw competitor counts.
## - fac_raw is raw NO_GRASS.
## - etalam is log(lambda).
## - etaai and etaaj are log(alpha).

## ---- Beverton-Holt + Gaussian raw
form_bh_gaussian <- bf(
  fitness_raw ~ exp(etalam) / (1 + exp(etaai) * intra_n + exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + Poisson count
form_bh_poisson <- bf(
  fitness_count ~ exp(etalam) / (1 + exp(etaai) * intra_n + exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + Negative binomial count
form_bh_negbin <- bf(
  fitness_count ~ exp(etalam) / (1 + exp(etaai) * intra_n + exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Beverton-Holt + Gaussian logfit
form_bh_loggaussian <- bf(
  fitness_log ~ log(
    exp(etalam) / (1 + exp(etaai) * intra_n + exp(etaaj) * inter_n) + 1
  ),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Gaussian raw
form_ricker_gaussian <- bf(
  fitness_raw ~ exp(etalam - exp(etaai) * intra_n - exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Poisson count
form_ricker_poisson <- bf(
  fitness_count ~ exp(etalam - exp(etaai) * intra_n - exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Negative binomial count
form_ricker_negbin <- bf(
  fitness_count ~ exp(etalam - exp(etaai) * intra_n - exp(etaaj) * inter_n),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

## ---- Ricker + Gaussian logfit
form_ricker_loggaussian <- bf(
  fitness_log ~ log(
    exp(etalam - exp(etaai) * intra_n - exp(etaaj) * inter_n) + 1
  ),
  etalam ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_raw +
    focsp:Site:fac_raw +
    focsp:Cov:fac_raw +
    (1 | SiteBlock),
  etaai ~ 0 + focsp,
  etaaj ~ 0 + focsp,
  nl = TRUE
)

############################
## 6. Priors
############################

## etalam is log(lambda).
##
## Important:
## - Baseline etalam terms need priors suitable for log fecundity.
## - Site/Cov deviations should be centred on 0.
## - Raw facilitator slopes should be centred on 0 and scaled to the raw
##   NO_GRASS unit.
## - etaai and etaaj are log(alpha), where alpha is competition strength
##   per competitor individual.
##
## The facilitator-slope prior below is scaled using the observed 90th
## percentile of raw NO_GRASS. This keeps raw units while avoiding an
## implausibly broad prior per one grass individual.

combine_priors <- function(...) {
  prior_list <- Filter(Negate(is.null), list(...))
  do.call(c, prior_list)
}

coef_priors <- function(coefs, prior_text, nlpar = "etalam") {
  coefs <- unique(coefs[!is.na(coefs) & nzchar(coefs)])
  
  if (length(coefs) == 0) {
    return(NULL)
  }
  
  do.call(
    c,
    lapply(coefs, function(cf) {
      brms::prior_string(
        prior_text,
        class = "b",
        nlpar = nlpar,
        coef = cf
      )
    })
  )
}

## Extract exact coefficient names from a representative formula.
## etalam coefficient names are the same across candidate model families.
prior_template <- brms::get_prior(
  formula = form_bh_gaussian,
  data    = dat_comp,
  family  = gaussian()
)

write.csv(
  as.data.frame(prior_template),
  file.path(base_dir, "summaries", "get_prior_template_bh_gaussian.csv"),
  row.names = FALSE
)

etalam_coefs <- prior_template %>%
  dplyr::filter(class == "b", nlpar == "etalam", coef != "") %>%
  dplyr::pull(coef)

## Species-specific baseline log(lambda) terms.
baseline_coefs <- intersect(etalam_coefs, paste0("focsp", sp_order))

## Raw facilitator slope terms, including species-specific,
## site-specific, and cover-specific facilitator effects.
fac_coefs <- etalam_coefs[stringr::str_detect(etalam_coefs, "fac_raw")]

## Remaining etalam coefficients are site/cover deviations.
other_etalam_coefs <- setdiff(etalam_coefs, c(baseline_coefs, fac_coefs))

## Set raw facilitator prior scale.
## Interpretation:
## - The prior SD is chosen so that a no-to-90th-percentile facilitator change
##   corresponds to roughly log(1.5) for one facilitator-slope component.
## - It is capped at 0.05 to avoid very broad per-grass priors if q90 is small.
fac_q90 <- as.numeric(quantile(dat_comp$fac_raw, 0.90, na.rm = TRUE))
fac_q90_for_prior <- max(fac_q90, 1)

fac_slope_prior_sd <- min(0.05, log(1.5) / fac_q90_for_prior)
fac_slope_prior_text <- sprintf("normal(0, %.8f)", fac_slope_prior_sd)

prior_record <- tibble::tibble(
  coef = etalam_coefs,
  prior_group = dplyr::case_when(
    coef %in% baseline_coefs     ~ "baseline_log_lambda",
    coef %in% fac_coefs          ~ "raw_facilitator_slope",
    coef %in% other_etalam_coefs ~ "site_cover_deviation",
    TRUE                         ~ "other"
  ),
  prior = dplyr::case_when(
    coef %in% baseline_coefs     ~ "normal(3, 2)",
    coef %in% fac_coefs          ~ fac_slope_prior_text,
    coef %in% other_etalam_coefs ~ "normal(0, 1)",
    TRUE                         ~ NA_character_
  )
)

write.csv(
  prior_record,
  file.path(base_dir, "summaries", "etalam_prior_coefficient_groups.csv"),
  row.names = FALSE
)

prior_scale_info <- list(
  baseline_log_lambda_prior = "normal(3, 2)",
  site_cover_deviation_prior = "normal(0, 1)",
  raw_facilitator_slope_prior = fac_slope_prior_text,
  raw_facilitator_q90 = fac_q90,
  raw_facilitator_q90_for_prior = fac_q90_for_prior,
  raw_facilitator_slope_prior_sd = fac_slope_prior_sd,
  alpha_log_prior = "normal(-3, 1)",
  note = paste(
    "fac_raw is raw NO_GRASS. It is not centred and not z-scaled.",
    "The facilitator prior is on the effect per one grass individual on log(lambda)."
  )
)

saveRDS(
  prior_scale_info,
  file.path(base_dir, "summaries", "prior_scale_info.rds")
)

capture.output(
  {
    cat("\nPrior scale information:\n")
    print(prior_scale_info)
    
    cat("\nEtalam prior coefficient groups:\n")
    print(prior_record)
  },
  file = file.path(base_dir, "summaries", "prior_scale_info.txt")
)

base_priors <- combine_priors(
  ## Default for site/cover deviations on log(lambda).
  ## A 1-unit coefficient is a ~2.7x multiplicative effect on lambda.
  brms::prior(normal(0, 1), class = "b", nlpar = "etalam"),
  
  ## Species-specific baseline log(lambda).
  ## normal(3, 2) is broad: exp(3) is ~20, with substantial mass much lower/higher.
  coef_priors(
    baseline_coefs,
    "normal(3, 2)",
    nlpar = "etalam"
  ),
  
  ## Raw facilitator effect per one grass individual.
  ## This is centred on zero because facilitation is being tested, not assumed.
  coef_priors(
    fac_coefs,
    fac_slope_prior_text,
    nlpar = "etalam"
  ),
  
  ## Competition coefficients.
  ## exp(-3) is ~0.05 per competitor individual.
  brms::prior(normal(-3, 1), class = "b", nlpar = "etaai"),
  brms::prior(normal(-3, 1), class = "b", nlpar = "etaaj"),
  
  ## Random block variation in etalam.
  brms::prior(exponential(1), class = "sd", nlpar = "etalam")
)

priors_gaussian <- c(
  base_priors,
  brms::prior(exponential(1), class = "sigma")
)

priors_poisson <- base_priors

priors_negbin <- c(
  base_priors,
  brms::prior(exponential(1), class = "shape")
)

priors_loggaussian <- c(
  base_priors,
  brms::prior(exponential(1), class = "sigma")
)

############################
## 7. Fit settings
############################

## Because priors and variable names changed, force a full refit now.
## After successful completion, set both to FALSE for normal reruns.
force_refit <- FALSE
force_recompute_loo <- FALSE

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
  
  if (file.exists(fit_rds) && !force_refit) {
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
  
  if (file.exists(loo_rds) && !force_recompute_loo) {
    cat("\nLoading existing LOO:", loo_name, "\n")
    loo_obj <- readRDS(loo_rds)
  } else {
    cat("\nComputing LOO:", loo_name, "\n")
    
    loo_obj <- loo(
      fit,
      moment_match = TRUE,
      recompile = TRUE,
      cores = fit_args$cores
    )
    
    saveRDS(loo_obj, loo_rds)
  }
  
  pareto_summary <- tibble(
    loo_name = loo_name,
    max_pareto_k = max(loo_obj$diagnostics$pareto_k, na.rm = TRUE),
    n_pareto_k_gt_0.7 = sum(loo_obj$diagnostics$pareto_k > 0.7, na.rm = TRUE),
    n_pareto_k_gt_1.0 = sum(loo_obj$diagnostics$pareto_k > 1.0, na.rm = TRUE)
  )
  
  write.csv(
    pareto_summary,
    file.path(base_dir, "loo", paste0(loo_name, "_pareto_summary.csv")),
    row.names = FALSE
  )
  
  capture.output(
    {
      cat("\nLOO summary for:", loo_name, "\n")
      print(loo_obj)
      cat("\nPareto-k summary:\n")
      print(pareto_summary)
    },
    file = file.path(base_dir, "loo", paste0(loo_name, "_summary.txt"))
  )
  
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

spec_lookup <- setNames(
  model_specs,
  vapply(model_specs, `[[`, character(1), "key")
)

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

loo_df <- as.data.frame(loo_tab) %>%
  tibble::rownames_to_column("model") %>%
  as_tibble()

write.csv(
  loo_df,
  file.path(base_dir, "selection", "loo_compare_family_error_models.csv"),
  row.names = FALSE
)

loo_close_4 <- loo_df %>%
  filter(elpd_diff >= -4)

write.csv(
  loo_close_4,
  file.path(base_dir, "selection", "loo_models_within_4_elpd.csv"),
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

pareto_all <- purrr::map_dfr(
  names(loo_list),
  function(nm) {
    tibble(
      model = nm,
      max_pareto_k = max(loo_list[[nm]]$diagnostics$pareto_k, na.rm = TRUE),
      n_pareto_k_gt_0.7 = sum(loo_list[[nm]]$diagnostics$pareto_k > 0.7, na.rm = TRUE),
      n_pareto_k_gt_1.0 = sum(loo_list[[nm]]$diagnostics$pareto_k > 1.0, na.rm = TRUE)
    )
  }
)

write.csv(
  pareto_all,
  file.path(base_dir, "selection", "pareto_k_summary_all_models.csv"),
  row.names = FALSE
)

############################
## 14. Record selected family + error model
############################

best_model <- rownames(loo_tab)[1]
best_spec  <- spec_lookup[[best_model]]

if (is.null(best_spec)) {
  stop("Could not match best_model from loo_compare to model_specs.")
}

selection_info <- list(
  best_model = best_model,
  best_growth_family = best_spec$growth_family,
  best_error_model = best_spec$error_model,
  best_fit_rds = file.path(base_dir, "fits", paste0(best_spec$fit_name, ".rds")),
  best_fit_name = best_spec$fit_name,
  comparison_csv = file.path(base_dir, "selection", "loo_compare_family_error_models.csv"),
  models_within_4_elpd_csv = file.path(base_dir, "selection", "loo_models_within_4_elpd.csv"),
  weights_csv = file.path(base_dir, "selection", "model_weights_family_error_models.csv"),
  pareto_summary_csv = file.path(base_dir, "selection", "pareto_k_summary_all_models.csv"),
  candidate_models_csv = file.path(base_dir, "selection", "candidate_family_error_models.csv"),
  data_rds = file.path(base_dir, "data", "analysis_data_family_error_compare.rds"),
  germination_priors_rds = normalizePath(germination_priors_file),
  germination_site_priors_rds = file.path(base_dir, "data", "germination_site_priors_used.rds"),
  prior_scale_info_rds = file.path(base_dir, "summaries", "prior_scale_info.rds"),
  raw_facilitator_rule = "fac_raw = NO_GRASS; not centred; not z-scaled; NO_GRASS = 0 means no facilitator.",
  competitor_density_rule = "intra_n and inter_n are raw competitor counts; alpha is per competitor individual.",
  germination_note = "Germination priors are loaded and validated here but are NOT used in stage-0 fitting; use them downstream for lambda_eff = g * lambda.",
  loo_caution = "LOO comparisons are used as a guide to predictive performance. Raw, count, and log-transformed response models are not perfectly comparable on a common likelihood scale."
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
  "1. Use the LOO comparison table as a guide to predictive performance, recognising that raw/count/log-transformed response models are not perfectly comparable on a common likelihood scale.",
  "2. Use pp_check plots, Pareto-k diagnostics, and summaries as complementary diagnostics.",
  "3. If models are very close, prefer the biologically appropriate and stable model rather than relying mechanically on the top-ranked row.",
  "4. Count-based models here are fit to fitness_count, which is only valid because fitness was checked as integer-valued.",
  "5. Germination priors were loaded and validated here, but not used in the stage-0 likelihood.",
  "6. Use germination downstream in coexistence calculations via lambda_eff = g * lambda.",
  "7. Competitor densities are raw counts, so alpha is interpretable as competition strength per competitor individual.",
  "8. Facilitator density is raw NO_GRASS, not centred and not z-scaled, so NO_GRASS = 0 means no facilitator.",
  "9. Facilitator effects enter etalam = log(lambda), so the facilitator slope is the effect of one additional grass individual on log(lambda).",
  "10. Downstream prediction scripts must use raw NO_GRASS values for no-facilitator and high-facilitator contrasts."
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
cat("\nImportant: force_refit and force_recompute_loo are currently TRUE.\n")
cat("After this full refit is complete, set them to FALSE for normal reruns.\n")