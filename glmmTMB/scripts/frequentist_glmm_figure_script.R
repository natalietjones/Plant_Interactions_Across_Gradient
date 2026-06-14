############################################################
## FINAL FREQUENTIST GLMM FIGURE SCRIPT
##
## Produces:
##   Figure 1: density-response plot with raw observations
##   Figure 2: predicted fecundity under no facilitator versus
##             high facilitator across all a priori contexts
##   Figure 3: selected-model summary of facilitation strength
##
## Required inputs:
##   - fullWAdata.csv in data/field/ preferred, data/ fallback
##   - final glmmTMB model objects saved in MODEL_DIR
##   - 15_final_model_selection_summary.csv in OUTDIR
##
## Key contrast used throughout:
##   - No facilitator = NO_GRASS = 0
##   - High facilitator = species-specific 90th percentile of NO_GRASS
##
## Notes:
##   - Fecundity remains a raw count response in fitted GLMMs.
##   - Model predictions use the final selected glmmTMB models.
##   - Figure 2 uses predicted seed number on the response scale.
##   - Figure 3 uses the link-scale log response ratio:
##       log(predicted high facilitator / predicted no facilitator)
##   - Density transformations/scaling are reconstructed from the
##     species-specific data and the density_form recorded in the model
##     selection table.
############################################################

############################################################
## 0. PACKAGES AND USER SETTINGS
############################################################

pkgs <- c(
  "tidyverse", "glmmTMB", "lme4", "ggh4x", "ggnewscale", "MASS", "scales", "patchwork", "car", "here", "readr"
)

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

## Match model-fitting script contrasts.
options(contrasts = c("contr.sum", "contr.poly"))

############################################################
## 0b. PROJECT PATHS
############################################################

## Match the fecundity model-fitting workflow path logic.
## Expected folder layout:
##   <project root>/data/field/fullWAdata.csv   preferred
##   <project root>/data/fullWAdata.csv         fallback
##   <project root>/glmmTMB/scripts/this_script.R
##   <project root>/glmmTMB/glmm_outputs/
##   <project root>/glmmTMB/figures_final_glmm/
##
## Robust to running from:
##   1. the main project root,
##   2. the glmmTMB folder, or
##   3. the glmmTMB/scripts folder.

normalise_existing <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

find_glmm_dir <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  here_root <- normalise_existing(here::here())

  candidates <- unique(normalise_existing(c(
    file.path(here_root, "glmmTMB"),
    here_root,
    wd,
    dirname(wd),
    dirname(dirname(wd))
  )))

  candidates <- candidates[dir.exists(candidates)]

  is_glmm_dir <- vapply(
    candidates,
    function(x) {
      basename(x) == "glmmTMB" ||
        dir.exists(file.path(x, "scripts")) ||
        dir.exists(file.path(x, "glmm_outputs")) ||
        dir.exists(file.path(x, "figures_final_glmm"))
    },
    logical(1)
  )

  hits <- candidates[is_glmm_dir]

  if (length(hits) == 0) {
    stop(
      "Could not locate the glmmTMB folder.\n",
      "Current working directory: ", wd, "\n",
      "here::here(): ", here_root, "\n",
      "Open the project from the parent folder or from glmmTMB, then rerun."
    )
  }

  hits[[1]]
}

resolve_first_existing <- function(paths, label) {
  paths <- unique(normalise_existing(paths))
  hit <- paths[file.exists(paths)]

  if (length(hit) == 0) {
    stop(
      "Could not find ", label, ". Checked:\n",
      paste(paths, collapse = "\n")
    )
  }

  hit[[1]]
}

glmm_dir <- find_glmm_dir()
project_root <- dirname(glmm_dir)

outdir <- file.path(glmm_dir, "glmm_outputs")
figdir <- file.path(glmm_dir, "figures_final_glmm")
model_dir <- file.path(outdir, "model_objects")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

## Prefer the current data convention first: data/field/fullWAdata.csv.
## Keep data/fullWAdata.csv as a fallback because some older scripts used it.
data_file <- resolve_first_existing(
  c(
    file.path(project_root, "data", "field", "fullWAdata.csv"),
    file.path(project_root, "data", "fullWAdata.csv"),
    file.path(glmm_dir, "data", "field", "fullWAdata.csv"),
    file.path(glmm_dir, "data", "fullWAdata.csv"),
    here::here("data", "field", "fullWAdata.csv"),
    here::here("data", "fullWAdata.csv"),
    here::here("fullWAdata.csv")
  ),
  label = "fullWAdata.csv"
)

## Uppercase aliases are used throughout this figure script.
DATA_FILE <- data_file
OUTDIR <- outdir
MODEL_DIR <- model_dir
FIGDIR <- figdir

cat("\nPath setup for figure script:\n")
cat("  glmm_dir:     ", normalizePath(glmm_dir, winslash = "/"), "\n", sep = "")
cat("  project_root: ", normalizePath(project_root, winslash = "/"), "\n", sep = "")
cat("  data_file:    ", normalizePath(DATA_FILE, winslash = "/"), "\n", sep = "")
cat("  output dir:   ", normalizePath(OUTDIR, winslash = "/"), "\n", sep = "")
cat("  model dir:    ", normalizePath(MODEL_DIR, winslash = "/"), "\n", sep = "")
cat("  figure dir:   ", normalizePath(FIGDIR, winslash = "/"), "\n\n", sep = "")

MODEL_SELECTION_FILE <- file.path(OUTDIR, "15_final_model_selection_summary.csv")
ANALYSIS_LABEL_TO_PLOT <- "primary_full_data_auto_density"
SPECIES_TO_PLOT <- c("TRCY", "TROR")
SITE_LEVELS <- c("Ben", "GH", "Nam", "PJ", "CS")
COV_LEVELS <- c("Sun", "Shade")

## Prediction settings.
N_GRID <- 160
FAC_HIGH_QUANTILE <- 0.90
LOW_DENSITY_QUANTILE <- 0.25
HIGH_DENSITY_QUANTILE <- 0.90
N_PARAM_DRAWS <- 2000
SEED <- 123
set.seed(SEED)

## Figure 1 x-axis setting.
## TRUE makes both density curves run across the full panel x-range.
PREDICT_TO_PANEL_MAX <- TRUE
X_PAD_MULT <- 1.02

## Figure sizes.
FIG_WIDTH <- 8
FIG_HEIGHT <- 6
FIG_DPI <- 600

## Graphic settings.
BASE_FAMILY <- "Arial"
BASE_SIZE <- 11
POINT_ALPHA <- 0.28
RAW_POINT_SIZE <- 0.9
RIBBON_ALPHA <- 0.18
LINE_WIDTH <- 0.75

## Raw-data overlay for Figure 2 predicted-fecundity contrasts.
## These points are observational fecundity values assigned to the plotted
## a priori neighbour/facilitator categories. Positive facilitator densities
## are assigned to the plotted "High facilitator" category because the
## model contrast uses 0 versus the positive 90th percentile.
RAW_OVERLAY_ALPHA <- 0.25
RAW_OVERLAY_SIZE <- 0.45
ERRORBAR_LINEWIDTH <- 0.48
ERRORBAR_WIDTH <- 0.05
RAW_OVERLAY_JITTER_WIDTH <- 0.025

fac_cols <- c(
  "No/low facilitator" = "#FFC43D",
  "High facilitator" = "#2C7FB8"
)

cover_cols <- c(
  "Sun" = "#D55E00",
  "Shade" = "#0072B2"
)

competition_cols <- c(
  "No neighbours" = "#170A1C",
  "Low conspecific" = "#067BC2",
  "Low heterospecific" = "#84BCDA",
  "High conspecific" = "#F37748",
  "High heterospecific" = "#D56062"
)

species_labs_chr <- c(
  "TRCY" = "italic('T. cyanopetala')",
  "TROR" = "italic('T. ornata')"
)

context_levels <- c(
  "No neighbours",
  "Low conspecific",
  "Low heterospecific",
  "High conspecific",
  "High heterospecific"
)

fac_labs <- c("No facilitator", "High facilitator")

############################################################
## 1. DATA PREP HELPERS
############################################################

z_safe <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

add_density_transforms <- function(dat) {
  dat %>%
    mutate(
      cons_log1p = log1p(conspecific_density),
      hetero_log1p = log1p(heterospecific_density),
      fac_log1p = log1p(facilitator_density)
    ) %>%
    group_by(focsp) %>%
    mutate(
      cons_raw_z = z_safe(conspecific_density),
      hetero_raw_z = z_safe(heterospecific_density),
      fac_raw_z = z_safe(facilitator_density),
      cons_log1p_z = z_safe(cons_log1p),
      hetero_log1p_z = z_safe(hetero_log1p),
      fac_log1p_z = z_safe(fac_log1p)
    ) %>%
    ungroup()
}

assign_density_form <- function(d, density_form = c("raw", "log1p")) {
  density_form <- match.arg(density_form)
  if (density_form == "raw") {
    d %>%
      mutate(
        cons_x = cons_raw_z,
        hetero_x = hetero_raw_z,
        fac_x = fac_raw_z
      )
  } else {
    d %>%
      mutate(
        cons_x = cons_log1p_z,
        hetero_x = hetero_log1p_z,
        fac_x = fac_log1p_z
      )
  }
}

prep_glmm_data <- function(file = DATA_FILE) {
  raw <- read.csv(file, stringsAsFactors = FALSE) %>%
    mutate(obs_id = row_number())
  
  required_cols <- c(
    "Site", "Cov", "BLK", "TROR", "TRCY", "NO_GRASS",
    "FOCAL_GERM", "focsp", "fitness"
  )
  
  missing_cols <- setdiff(required_cols, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  raw %>%
    mutate(
      FOCAL_GERM_BIN = case_when(
        FOCAL_GERM %in% c("Y", "Yes", "YES", "y", "yes", 1, "1") ~ 1,
        FOCAL_GERM %in% c("N", "No", "NO", "n", "no", 0, "0") ~ 0,
        TRUE ~ NA_real_
      ),
      Site = case_when(
        Site %in% c("BEN", "Ben", "ben") ~ "Ben",
        Site %in% c("GH", "Gh", "gh") ~ "GH",
        Site %in% c("NAM", "Nam", "nam") ~ "Nam",
        Site %in% c("PJ", "Pj", "pj") ~ "PJ",
        Site %in% c("CS", "Cs", "cs") ~ "CS",
        TRUE ~ as.character(Site)
      ),
      Site = factor(Site, levels = SITE_LEVELS),
      Cov = case_when(
        Cov %in% c("SUN", "Sun", "sun") ~ "Sun",
        Cov %in% c("SH", "SHADE", "Shade", "shade", "sh") ~ "Shade",
        TRUE ~ as.character(Cov)
      ),
      Cov = factor(Cov, levels = COV_LEVELS),
      focsp = factor(focsp, levels = SPECIES_TO_PLOT),
      BLK = factor(BLK),
      fecundity = as.numeric(fitness),
      conspecific_density = case_when(
        focsp == "TRCY" ~ as.numeric(TRCY),
        focsp == "TROR" ~ as.numeric(TROR),
        TRUE ~ NA_real_
      ),
      heterospecific_density = case_when(
        focsp == "TRCY" ~ as.numeric(TROR),
        focsp == "TROR" ~ as.numeric(TRCY),
        TRUE ~ NA_real_
      ),
      facilitator_density = as.numeric(NO_GRASS),
      SiteBlock = interaction(Site, BLK, drop = TRUE)
    ) %>%
    filter(
      FOCAL_GERM_BIN == 1,
      !is.na(fecundity),
      fecundity >= 0,
      !is.na(conspecific_density),
      !is.na(heterospecific_density),
      !is.na(facilitator_density),
      conspecific_density >= 0,
      heterospecific_density >= 0,
      facilitator_density >= 0,
      !is.na(Site),
      !is.na(Cov),
      !is.na(focsp),
      !is.na(SiteBlock)
    ) %>%
    mutate(zero_fecundity = fecundity == 0) %>%
    droplevels() %>%
    add_density_transforms()
}

############################################################
## 2. LOAD MODELS AND MODEL-SELECTION METADATA
############################################################

load_model_for_species <- function(species_code) {
  candidate_paths <- c(
    file.path(MODEL_DIR, paste0(species_code, "_primary_final_model.rds")),
    file.path(MODEL_DIR, paste0("fit_", species_code, "_final_model.rds")),
    file.path(MODEL_DIR, paste0(species_code, "_", ANALYSIS_LABEL_TO_PLOT, "_final_model.rds")),
    file.path(OUTDIR, paste0(species_code, "_primary_final_model.rds")),
    file.path(OUTDIR, paste0("fit_", species_code, "_final_model.rds"))
  )
  
  existing <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing) == 0) {
    stop(
      "Could not find a saved final model for ", species_code, ".\n",
      "Expected one of:\n  ", paste(candidate_paths, collapse = "\n  "), "\n\n",
      "Add model-saving code to the end of your fitting workflow, rerun once, then rerun this script."
    )
  }
  
  message("Loading ", species_code, " model from: ", existing[[1]])
  model <- readRDS(existing[[1]])
  if (!inherits(model, "glmmTMB")) {
    stop("Saved object for ", species_code, " is not a glmmTMB model: ", existing[[1]])
  }
  model
}

if (!file.exists(MODEL_SELECTION_FILE)) {
  stop("Cannot find model-selection summary: ", MODEL_SELECTION_FILE)
}

model_selection <- readr::read_csv(MODEL_SELECTION_FILE, show_col_types = FALSE) %>%
  filter(analysis_label == ANALYSIS_LABEL_TO_PLOT, species %in% SPECIES_TO_PLOT) %>%
  mutate(
    species = as.character(species),
    density_form = as.character(density_form)
  )

if (nrow(model_selection) != length(SPECIES_TO_PLOT)) {
  stop(
    "Model-selection table does not contain exactly one row per species for analysis_label = ",
    ANALYSIS_LABEL_TO_PLOT, ". Check ", MODEL_SELECTION_FILE
  )
}

models <- setNames(
  lapply(SPECIES_TO_PLOT, load_model_for_species),
  SPECIES_TO_PLOT
)

cat("\nLoaded model formulas:\n")
purrr::iwalk(models, function(mod, sp) {
  cat("\n", sp, ":\n", sep = "")
  print(formula(mod))
})

df_glmm <- prep_glmm_data(DATA_FILE)

get_species_meta <- function(species_code) {
  model_selection %>%
    filter(.data$species == species_code) %>%
    slice(1)
}

get_species_data <- function(species_code) {
  df_glmm %>%
    filter(focsp == species_code) %>%
    droplevels() %>%
    dplyr::select(
      -dplyr::any_of(c(
        "cons_raw_z", "hetero_raw_z", "fac_raw_z",
        "cons_log1p", "hetero_log1p", "fac_log1p",
        "cons_log1p_z", "hetero_log1p_z", "fac_log1p_z",
        "cons_x", "hetero_x", "fac_x"
      ))
    ) %>%
    add_density_transforms() %>%
    assign_density_form(get_species_meta(species_code)$density_form[[1]])
}

############################################################
## 3. NEWDATA SCALING AND PREDICTION HELPERS
############################################################

get_scalers <- function(d_sp, density_form) {
  density_form <- match.arg(density_form, c("raw", "log1p"))
  
  if (density_form == "raw") {
    list(
      cons_mu = mean(d_sp$conspecific_density, na.rm = TRUE),
      cons_sd = sd(d_sp$conspecific_density, na.rm = TRUE),
      hetero_mu = mean(d_sp$heterospecific_density, na.rm = TRUE),
      hetero_sd = sd(d_sp$heterospecific_density, na.rm = TRUE),
      fac_mu = mean(d_sp$facilitator_density, na.rm = TRUE),
      fac_sd = sd(d_sp$facilitator_density, na.rm = TRUE),
      transform = "raw"
    )
  } else {
    list(
      cons_mu = mean(log1p(d_sp$conspecific_density), na.rm = TRUE),
      cons_sd = sd(log1p(d_sp$conspecific_density), na.rm = TRUE),
      hetero_mu = mean(log1p(d_sp$heterospecific_density), na.rm = TRUE),
      hetero_sd = sd(log1p(d_sp$heterospecific_density), na.rm = TRUE),
      fac_mu = mean(log1p(d_sp$facilitator_density), na.rm = TRUE),
      fac_sd = sd(log1p(d_sp$facilitator_density), na.rm = TRUE),
      transform = "log1p"
    )
  }
}

scale_one <- function(x, mu, sd, transform = c("raw", "log1p")) {
  transform <- match.arg(transform)
  xx <- if (transform == "raw") x else log1p(x)
  if (is.na(sd) || sd == 0) return(rep(0, length(x)))
  (xx - mu) / sd
}

add_scaled_predictors <- function(nd, scalers) {
  nd %>%
    mutate(
      cons_x = scale_one(conspecific_density, scalers$cons_mu, scalers$cons_sd, scalers$transform),
      hetero_x = scale_one(heterospecific_density, scalers$hetero_mu, scalers$hetero_sd, scalers$transform),
      fac_x = scale_one(facilitator_density, scalers$fac_mu, scalers$fac_sd, scalers$transform),
      Site = factor(Site, levels = SITE_LEVELS),
      Cov = factor(Cov, levels = COV_LEVELS)
    )
}

predict_link_response <- function(model, newdata) {
  pp <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  if (is.list(pp)) {
    fit_link <- as.numeric(pp$fit)
    se_link <- as.numeric(pp$se.fit)
  } else {
    fit_link <- as.numeric(pp)
    se_link <- rep(NA_real_, length(fit_link))
  }
  
  tibble(
    fit_link = fit_link,
    se_link = se_link,
    pred = exp(fit_link),
    pred_low = exp(fit_link - 1.96 * se_link),
    pred_high = exp(fit_link + 1.96 * se_link)
  )
}

predict_response_ci <- function(model, newdata) {
  ## glmmTMB versions differ in whether predict(..., type = "response", se.fit = TRUE)
  ## returns a list or an atomic vector. Use the link scale for a stable CI,
  ## then transform back to predicted seed number.
  pp <- predict(
    model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  if (!is.list(pp)) {
    stop(
      "predict(..., se.fit = TRUE) did not return standard errors. ",
      "Update glmmTMB or use the parametric-draw prediction code instead."
    )
  }
  
  fit_link <- as.numeric(pp$fit)
  se_link <- as.numeric(pp$se.fit)
  
  tibble(
    fit_link = fit_link,
    se_link = se_link,
    pred = exp(fit_link),
    pred_low = exp(fit_link - 1.96 * se_link),
    pred_high = exp(fit_link + 1.96 * se_link)
  )
}

align_model_matrix <- function(model, newdata) {
  fixed_form <- lme4::nobars(formula(model))
  trm <- delete.response(terms(fixed_form))
  X <- model.matrix(trm, newdata)
  
  beta_names <- names(fixef(model)$cond)
  missing_cols <- setdiff(beta_names, colnames(X))
  extra_cols <- setdiff(colnames(X), beta_names)
  
  if (length(missing_cols) > 0) {
    X_missing <- matrix(0, nrow = nrow(X), ncol = length(missing_cols))
    colnames(X_missing) <- missing_cols
    X <- cbind(X, X_missing)
  }
  
  if (length(extra_cols) > 0) {
    X <- X[, setdiff(colnames(X), extra_cols), drop = FALSE]
  }
  
  X[, beta_names, drop = FALSE]
}

draw_fixed_response <- function(model, newdata, n_draws = N_PARAM_DRAWS) {
  beta <- fixef(model)$cond
  V <- as.matrix(vcov(model)$cond)
  X <- align_model_matrix(model, newdata)
  
  draws <- MASS::mvrnorm(n = n_draws, mu = beta, Sigma = V)
  eta <- X %*% t(draws)
  exp(eta)
}

draw_fixed_expected_response <- function(model, newdata, n_draws = N_PARAM_DRAWS) {
  ## Draw predicted expected seed number from the fixed-effect uncertainty.
  ## For zero-inflated models with ziformula = ~1, this returns:
  ##   E[y] = conditional_mean * (1 - zero_inflation_probability)
  ## If the zero-inflation component is more complex than an intercept-only model,
  ## the function stops rather than silently giving the wrong expected response.
  
  cond_mu_draws <- draw_fixed_response(model, newdata, n_draws = n_draws)
  
  zi_beta <- fixef(model)$zi
  
  if (length(zi_beta) == 0) {
    return(cond_mu_draws)
  }
  
  if (!(length(zi_beta) == 1 && names(zi_beta)[1] == "(Intercept)")) {
    stop(
      "draw_fixed_expected_response() currently handles intercept-only zero inflation. ",
      "The fitted model has a non-intercept-only zero-inflation component."
    )
  }
  
  zi_V <- as.matrix(vcov(model)$zi)
  zi_draws <- MASS::mvrnorm(n = n_draws, mu = zi_beta, Sigma = zi_V)
  zi_eta <- as.numeric(zi_draws[, 1])
  zi_prob <- plogis(zi_eta)
  
  sweep(cond_mu_draws, 2, 1 - zi_prob, `*`)
}

summarise_prediction_draws <- function(pred_draws) {
  ## Summary of response-scale prediction draws.
  ## For Figure 2, plot mean ± 1 SE so the error bars are symmetric
  ## around the plotted point. Quantile intervals are retained in the
  ## output CSV because they are often preferable for GLMM uncertainty,
  ## but they are asymmetric on the response scale after the log-link
  ## back-transformation.
  
  pred_mean <- apply(pred_draws, 1, mean, na.rm = TRUE)
  pred_median <- apply(pred_draws, 1, median, na.rm = TRUE)
  pred_se_draw <- apply(pred_draws, 1, sd, na.rm = TRUE)
  pred_q025 <- apply(pred_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  pred_q975 <- apply(pred_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  tibble(
    pred = pred_mean,
    pred_mean = pred_mean,
    pred_median = pred_median,
    pred_se_draw = pred_se_draw,
    pred_low = pmax(0, pred_mean - pred_se_draw),
    pred_high = pred_mean + pred_se_draw,
    pred_q025 = pred_q025,
    pred_q975 = pred_q975,
    pred_ci_width = pred_high - pred_low,
    pred_ci_ratio = pred_high / pmax(pred, 1e-12)
  )
}

############################################################
## 4. COMMON FACILITATOR AND DENSITY CONTEXT HELPERS
############################################################

q_positive <- function(x, prob = 0.25) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(0)
  as.numeric(quantile(x, probs = prob, na.rm = TRUE))
}

make_fac_table <- function(d_sp) {
  ## Use species-specific high facilitator density for final-model plots.
  ## No facilitator is always exactly NO_GRASS = 0.
  tmp <- d_sp %>%
    summarise(
      fac_low = 0,
      fac_high = q_positive(facilitator_density, FAC_HIGH_QUANTILE),
      fac_max = max(facilitator_density, na.rm = TRUE)
    ) %>%
    mutate(fac_high = if_else(fac_high <= fac_low, fac_max, fac_high))
  
  tibble(Site = factor(SITE_LEVELS, levels = SITE_LEVELS)) %>%
    bind_cols(tmp[rep(1, length(SITE_LEVELS)), ])
}

make_context_table_by_site <- function(d_sp) {
  d_sp %>%
    group_by(Site) %>%
    summarise(
      cons_low = q_positive(conspecific_density, LOW_DENSITY_QUANTILE),
      hetero_low = q_positive(heterospecific_density, LOW_DENSITY_QUANTILE),
      cons_high = q_positive(conspecific_density, HIGH_DENSITY_QUANTILE),
      hetero_high = q_positive(heterospecific_density, HIGH_DENSITY_QUANTILE),
      .groups = "drop"
    )
}

make_context_table_global <- function(d_sp) {
  tibble(
    competition_context = factor(context_levels, levels = context_levels),
    conspecific_density = c(
      0,
      q_positive(d_sp$conspecific_density, LOW_DENSITY_QUANTILE),
      0,
      q_positive(d_sp$conspecific_density, HIGH_DENSITY_QUANTILE),
      0
    ),
    heterospecific_density = c(
      0,
      0,
      q_positive(d_sp$heterospecific_density, LOW_DENSITY_QUANTILE),
      0,
      q_positive(d_sp$heterospecific_density, HIGH_DENSITY_QUANTILE)
    )
  )
}

############################################################
## 5. FIGURE 1: DENSITY-RESPONSE PLOT
############################################################

make_density_prediction_grid <- function(species_code, density_axis = c("conspecific", "heterospecific")) {
  density_axis <- match.arg(density_axis)
  d_sp <- get_species_data(species_code)
  meta <- get_species_meta(species_code)
  scalers <- get_scalers(d_sp, meta$density_form[[1]])
  fac_table <- make_fac_table(d_sp)
  
  density_limits <- d_sp %>%
    group_by(Site, Cov) %>%
    summarise(
      cons_data_max = max(conspecific_density, na.rm = TRUE),
      hetero_data_max = max(heterospecific_density, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      panel_density_max = if (PREDICT_TO_PANEL_MAX) {
        pmax(cons_data_max, hetero_data_max, na.rm = TRUE)
      } else if (density_axis == "conspecific") {
        cons_data_max
      } else {
        hetero_data_max
      },
      panel_density_max = if_else(
        !is.finite(panel_density_max) | panel_density_max <= 0,
        1,
        panel_density_max * X_PAD_MULT
      )
    )
  
  density_axis_label <- if (density_axis == "conspecific") {
    "Conspecific density"
  } else {
    "Heterospecific density"
  }
  
  grid <- tidyr::crossing(
    Site = factor(SITE_LEVELS, levels = SITE_LEVELS),
    Cov = factor(COV_LEVELS, levels = COV_LEVELS),
    grid_id = seq_len(N_GRID),
    fac_level = c("No/low facilitator", "High facilitator")
  ) %>%
    left_join(
      density_limits %>% dplyr::select(Site, Cov, panel_density_max),
      by = c("Site", "Cov")
    ) %>%
    left_join(
      fac_table %>% dplyr::select(Site, fac_low, fac_high),
      by = "Site"
    ) %>%
    mutate(
      density_raw = (grid_id - 1) / (N_GRID - 1) * panel_density_max,
      conspecific_density = if (density_axis == "conspecific") density_raw else 0,
      heterospecific_density = if (density_axis == "heterospecific") density_raw else 0,
      facilitator_density = if_else(fac_level == "High facilitator", fac_high, fac_low),
      density_axis = density_axis_label,
      species = species_code,
      focsp = factor(species_code, levels = SPECIES_TO_PLOT),
      BLK = factor(levels(d_sp$BLK)[1], levels = levels(d_sp$BLK)),
      SiteBlock = factor(levels(d_sp$SiteBlock)[1], levels = levels(d_sp$SiteBlock))
    ) %>%
    add_scaled_predictors(scalers)
  
  pred <- predict_link_response(models[[species_code]], grid)
  bind_cols(grid, pred)
}

density_pred <- purrr::map_dfr(SPECIES_TO_PLOT, function(sp) {
  bind_rows(
    make_density_prediction_grid(sp, "conspecific"),
    make_density_prediction_grid(sp, "heterospecific")
  )
}) %>%
  mutate(
    species = factor(species, levels = SPECIES_TO_PLOT),
    density_axis = factor(density_axis, levels = c("Conspecific density", "Heterospecific density")),
    fac_level = factor(fac_level, levels = c("No/low facilitator", "High facilitator")),
    panel_row = factor(
      paste(species, Cov, sep = " | "),
      levels = c("TRCY | Sun", "TROR | Sun", "TRCY | Shade", "TROR | Shade")
    )
  )

raw_long <- df_glmm %>%
  filter(focsp %in% SPECIES_TO_PLOT) %>%
  transmute(
    species = factor(as.character(focsp), levels = SPECIES_TO_PLOT),
    Site = factor(Site, levels = SITE_LEVELS),
    Cov = factor(Cov, levels = COV_LEVELS),
    fecundity,
    conspecific_density,
    heterospecific_density
  ) %>%
  pivot_longer(
    cols = c(conspecific_density, heterospecific_density),
    names_to = "density_axis",
    values_to = "density_raw"
  ) %>%
  mutate(
    density_axis = dplyr::recode(
      density_axis,
      conspecific_density = "Conspecific density",
      heterospecific_density = "Heterospecific density"
    ),
    density_axis = factor(density_axis, levels = c("Conspecific density", "Heterospecific density")),
    panel_row = factor(
      paste(species, Cov, sep = " | "),
      levels = c("TRCY | Sun", "TROR | Sun", "TRCY | Shade", "TROR | Shade")
    )
  )

fig1_combined <- ggplot() +
  geom_point(
    data = raw_long,
    aes(x = density_raw, y = fecundity),
    alpha = POINT_ALPHA,
    size = RAW_POINT_SIZE,
    colour = "grey20",
    stroke = 0
  ) +
  geom_ribbon(
    data = density_pred,
    aes(
      x = density_raw,
      ymin = pred_low,
      ymax = pred_high,
      fill = fac_level,
      group = interaction(fac_level, density_axis)
    ),
    alpha = RIBBON_ALPHA,
    colour = NA
  ) +
  geom_line(
    data = density_pred,
    aes(
      x = density_raw,
      y = pred,
      linetype = density_axis,
      colour = fac_level,
      group = interaction(fac_level, density_axis)
    ),
    linewidth = LINE_WIDTH
  ) +
  ggh4x::facet_grid2(
    rows = vars(panel_row),
    cols = vars(Site),
    scales = "free",
    independent = "all",
    axes = "all",
    remove_labels = "none"
  ) +
  scale_colour_manual(values = fac_cols) +
  scale_fill_manual(values = fac_cols) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(
    x = "Neighbour density",
    y = "Predicted fecundity",
    colour = "Facilitator level",
    fill = "Facilitator level",
    linetype = "Density varied"
  ) +
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 8.5, colour = "black"),
    axis.text.y = element_text(size = 8.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.2, "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.35, colour = "black", fill = NA),
    panel.spacing.x = unit(0.08, "lines"),
    panel.spacing.y = unit(0.08, "lines"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text.x = element_text(size = 8.5, face = "bold", margin = margin(1, 1, 1, 1)),
    strip.text.y = element_text(size = 8.0, face = "bold", margin = margin(1, 1, 1, 1)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9.5),
    legend.key.height = unit(0.35, "lines"),
    legend.key.width = unit(0.9, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-4, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 0.35)),
    colour = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

print(fig1_combined)

# ggsave(
#   file.path(FIGDIR, "figure1_density_responses_combined_A4.pdf"),
#   fig1_combined,
#   width = FIG_WIDTH,
#   height = FIG_HEIGHT,
#   device = cairo_pdf
# )

ggsave(
  file.path(FIGDIR, "figure1_density_responses_combined_A4.png"),
  fig1_combined,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = FIG_DPI
)

#write_csv(
#  density_pred,
#  file.path(FIGDIR, "figure1_density_response_predictions_A4.csv")
#)

############################################################
## 6. FIGURE 2: PREDICTED FECUNDITY ACROSS ALL A PRIORI
## FACILITATOR CONTRASTS
############################################################

make_predicted_fecundity_grid <- function(species_code) {
  d_sp <- get_species_data(species_code)
  meta <- get_species_meta(species_code)
  scalers <- get_scalers(d_sp, meta$density_form[[1]])
  fac_table <- make_fac_table(d_sp)
  context_table <- make_context_table_by_site(d_sp)
  
  grid_base <- tidyr::crossing(
    Site = factor(SITE_LEVELS, levels = SITE_LEVELS),
    Cov = factor(COV_LEVELS, levels = COV_LEVELS),
    competition_context = factor(context_levels, levels = context_levels),
    fac_level = factor(fac_labs, levels = fac_labs)
  ) %>%
    left_join(context_table, by = "Site") %>%
    left_join(fac_table %>% dplyr::select(Site, fac_low, fac_high), by = "Site") %>%
    mutate(
      conspecific_density = case_when(
        competition_context == "Low conspecific" ~ cons_low,
        competition_context == "High conspecific" ~ cons_high,
        TRUE ~ 0
      ),
      heterospecific_density = case_when(
        competition_context == "Low heterospecific" ~ hetero_low,
        competition_context == "High heterospecific" ~ hetero_high,
        TRUE ~ 0
      ),
      facilitator_density = if_else(
        fac_level == "High facilitator",
        fac_high,
        fac_low
      ),
      species = species_code,
      species_label = species_labs_chr[[species_code]],
      focsp = factor(species_code, levels = SPECIES_TO_PLOT),
      BLK = factor(levels(d_sp$BLK)[1], levels = levels(d_sp$BLK)),
      SiteBlock = factor(levels(d_sp$SiteBlock)[1], levels = levels(d_sp$SiteBlock))
    ) %>%
    add_scaled_predictors(scalers)
  
  pred_draws <- draw_fixed_expected_response(
    models[[species_code]],
    grid_base,
    n_draws = N_PARAM_DRAWS
  )
  
  pred <- summarise_prediction_draws(pred_draws)
  bind_cols(grid_base, pred)
}

fecundity_pred <- purrr::map_dfr(
  SPECIES_TO_PLOT,
  make_predicted_fecundity_grid
) %>%
  mutate(
    species = factor(species, levels = SPECIES_TO_PLOT),
    species_label = factor(species_label, levels = species_labs_chr[SPECIES_TO_PLOT]),
    Site = factor(Site, levels = SITE_LEVELS),
    Cov = factor(Cov, levels = COV_LEVELS),
    competition_context = factor(competition_context, levels = context_levels),
    fac_level = factor(fac_level, levels = fac_labs)
  )

# write_csv(
#   fecundity_pred,
#   file.path(FIGDIR, "figure2_predicted_fecundity_all_contrasts.csv")
# )

figure2_interval_diagnostics <- fecundity_pred %>%
  arrange(desc(pred_ci_ratio), desc(pred_ci_width)) %>%
  dplyr::select(
    species, Site, Cov, competition_context, fac_level,
    conspecific_density, heterospecific_density, facilitator_density,
    pred, pred_low, pred_high, pred_se_draw, pred_q025, pred_q975, pred_ci_width, pred_ci_ratio
  )

print(figure2_interval_diagnostics, n = 30)

# write_csv(
#   figure2_interval_diagnostics,
#   file.path(FIGDIR, "figure2_prediction_interval_diagnostics.csv")
# )

############################################################
## Plot Figure 2: four-row version
## Rows:
##   (a) TRCY Sun
##   (b) TRCY Shade
##   (c) TROR Sun
##   (d) TROR Shade
############################################################

context_offsets <- tibble(
  competition_context = factor(context_levels, levels = context_levels),
  context_offset = seq(-0.36, 0.36, length.out = length(context_levels))
)

plot_offsets <- tidyr::crossing(
  Cov = factor(COV_LEVELS, levels = COV_LEVELS),
  competition_context = factor(context_levels, levels = context_levels),
  fac_level = factor(fac_labs, levels = fac_labs)
) %>%
  left_join(context_offsets, by = "competition_context") %>%
  mutate(
    fac_offset = if_else(fac_level == "No facilitator", -0.045, 0.045),
    x_offset = context_offset + fac_offset
  ) %>%
  dplyr::select(Cov, competition_context, fac_level, x_offset)

fecundity_pred_plot <- fecundity_pred %>%
  left_join(plot_offsets, by = c("Cov", "competition_context", "fac_level")) %>%
  mutate(
    x_num = as.numeric(Site) + x_offset,
    segment_id = interaction(species, Site, Cov, competition_context, drop = TRUE),
    pred_plot = pred,
    pred_low_plot = pmax(pred_low, 0),
    pred_high_plot = pred_high,
    panel_label = case_when(
      species == "TRCY" & Cov == "Sun"   ~ "'(a)'~italic('T. cyanopetala')~'Sun'",
      species == "TRCY" & Cov == "Shade" ~ "'(b)'~italic('T. cyanopetala')~'Shade'",
      species == "TROR" & Cov == "Sun"   ~ "'(c)'~italic('T. ornata')~'Sun'",
      species == "TROR" & Cov == "Shade" ~ "'(d)'~italic('T. ornata')~'Shade'",
      TRUE ~ as.character(species)
    ),
    panel_label = factor(
      panel_label,
      levels = c(
        "'(a)'~italic('T. cyanopetala')~'Sun'",
        "'(b)'~italic('T. cyanopetala')~'Shade'",
        "'(c)'~italic('T. ornata')~'Sun'",
        "'(d)'~italic('T. ornata')~'Shade'"
      )
    )
  )

############################################################
## Raw fecundity observations assigned to Figure 2 conditions
##
## These are used only as a faded data overlay. The model predictions are
## defined at exact a priori contrast values: no facilitator = 0, high
## facilitator = positive 90th percentile, and low/high neighbour densities
## = positive 25th/90th percentiles. Raw observations rarely fall exactly
## on all of those values, so each observation is assigned to the closest
## plotted single-neighbour context within its species x site x cover cell.
## Observations with both conspecific and heterospecific density > 0 are
## not shown in the overlay because Figure 2 varies only one neighbour
## axis at a time.
############################################################

assign_raw_context_for_plot <- function(species_code) {
  d_sp <- get_species_data(species_code)

  context_table <- make_context_table_by_site(d_sp)
  fac_table <- make_fac_table(d_sp)

  d_sp %>%
    left_join(context_table, by = "Site") %>%
    left_join(
      fac_table %>% dplyr::select(Site, fac_low, fac_high),
      by = "Site"
    ) %>%
    mutate(
      species = as.character(focsp),
      species_label = species_labs_chr[[species_code]],
      fac_level = case_when(
        facilitator_density == 0 ~ "No facilitator",
        facilitator_density > 0 ~ "High facilitator",
        TRUE ~ NA_character_
      ),
      competition_context = case_when(
        conspecific_density == 0 & heterospecific_density == 0 ~ "No neighbours",
        heterospecific_density == 0 &
          conspecific_density > 0 &
          abs(conspecific_density - cons_low) <= abs(conspecific_density - cons_high) ~ "Low conspecific",
        heterospecific_density == 0 & conspecific_density > 0 ~ "High conspecific",
        conspecific_density == 0 &
          heterospecific_density > 0 &
          abs(heterospecific_density - hetero_low) <= abs(heterospecific_density - hetero_high) ~ "Low heterospecific",
        conspecific_density == 0 & heterospecific_density > 0 ~ "High heterospecific",
        TRUE ~ NA_character_
      ),
      raw_overlay_status = case_when(
        is.na(competition_context) ~ "not_plotted_both_neighbour_axes_positive_or_unclassified",
        is.na(fac_level) ~ "not_plotted_facilitator_unclassified",
        TRUE ~ "plotted"
      ),
      raw_overlay_note = case_when(
        raw_overlay_status != "plotted" ~ "Observation not shown in Figure 2 raw overlay",
        facilitator_density > 0 ~ "Positive facilitator observation assigned to High facilitator plotting category",
        TRUE ~ "No facilitator observation"
      )
    )
}

raw_fecundity_context_all <- purrr::map_dfr(
  SPECIES_TO_PLOT,
  assign_raw_context_for_plot
) %>%
  mutate(
    species = factor(species, levels = SPECIES_TO_PLOT),
    species_label = factor(species_label, levels = species_labs_chr[SPECIES_TO_PLOT]),
    Site = factor(Site, levels = SITE_LEVELS),
    Cov = factor(Cov, levels = COV_LEVELS),
    competition_context = factor(competition_context, levels = context_levels),
    fac_level = factor(fac_level, levels = fac_labs)
  )

raw_fecundity_context_plot <- raw_fecundity_context_all %>%
  filter(raw_overlay_status == "plotted") %>%
  left_join(
    plot_offsets,
    by = c("Cov", "competition_context", "fac_level")
  ) %>%
  mutate(
    x_num = as.numeric(Site) + x_offset,
    x_num_raw = x_num + stats::runif(
      dplyr::n(),
      min = -RAW_OVERLAY_JITTER_WIDTH,
      max = RAW_OVERLAY_JITTER_WIDTH
    ),
    panel_label = case_when(
      species == "TRCY" & Cov == "Sun"   ~ "'(a)'~italic('T. cyanopetala')~'Sun'",
      species == "TRCY" & Cov == "Shade" ~ "'(b)'~italic('T. cyanopetala')~'Shade'",
      species == "TROR" & Cov == "Sun"   ~ "'(c)'~italic('T. ornata')~'Sun'",
      species == "TROR" & Cov == "Shade" ~ "'(d)'~italic('T. ornata')~'Shade'",
      TRUE ~ as.character(species)
    ),
    panel_label = factor(panel_label, levels = levels(fecundity_pred_plot$panel_label))
  )

raw_fecundity_context_omitted <- raw_fecundity_context_all %>%
  filter(raw_overlay_status != "plotted")

raw_fecundity_overlay_counts <- raw_fecundity_context_all %>%
  count(
    species, Site, Cov, competition_context, fac_level,
    raw_overlay_status,
    name = "n_raw_observations"
  ) %>%
  arrange(species, Site, Cov, raw_overlay_status, competition_context, fac_level)

print(raw_fecundity_overlay_counts, n = Inf)

# write_csv(
#   raw_fecundity_context_plot,
#   file.path(FIGDIR, "figure2_raw_fecundity_overlay_points.csv")
# )

# write_csv(
#   raw_fecundity_context_omitted,
#   file.path(FIGDIR, "figure2_raw_fecundity_overlay_omitted_observations.csv")
# )

# write_csv(
#   raw_fecundity_overlay_counts,
#   file.path(FIGDIR, "figure2_raw_fecundity_overlay_counts.csv")
# )

figure2_y_axis_diagnostics <- fecundity_pred_plot %>%
  group_by(species, Cov, panel_label) %>%
  summarise(
    n_points = n(),
    max_prediction = max(pred, na.rm = TRUE),
    max_upper_interval = max(pred_high, na.rm = TRUE),
    min_lower_interval = min(pred_low, na.rm = TRUE),
    max_upper_quantile_interval = max(pred_q975, na.rm = TRUE),
    min_lower_quantile_interval = min(pred_q025, na.rm = TRUE),
    .groups = "drop"
  )

print(figure2_y_axis_diagnostics)

# write_csv(
#   figure2_y_axis_diagnostics,
#   file.path(FIGDIR, "figure2_y_axis_diagnostics_no_capping.csv")
# )

figure2_se_symmetry_check <- fecundity_pred_plot %>%
  mutate(
    se_above = pred_high_plot - pred_plot,
    se_below = pred_plot - pred_low_plot,
    se_difference = se_above - se_below,
    lower_truncated_at_zero = pred_low < 0
  ) %>%
  arrange(desc(abs(se_difference))) %>%
  dplyr::select(
    species, Cov, Site, competition_context, fac_level,
    pred, pred_low, pred_high, pred_se_draw,
    se_above, se_below, se_difference, lower_truncated_at_zero
  )

print(figure2_se_symmetry_check, n = 30)

# write_csv(
#   figure2_se_symmetry_check,
#   file.path(FIGDIR, "figure2_se_symmetry_check.csv")
# )

fecundity_pred_plot %>%
  mutate(
    err_height = pred_high_plot - pred_low_plot,
    err_small = err_height < 1
  ) %>%
  arrange(err_height) %>%
  dplyr::select(
    species, Cov, Site, competition_context, fac_level,
    pred_plot, pred_low_plot, pred_high_plot, err_height, err_small
  ) %>%
  print(n = 40)

fig2_predicted_fecundity <- ggplot(
  fecundity_pred_plot,
  aes(
    x = x_num,
    y = pred_plot,
    ymin = pred_low_plot,
    ymax = pred_high_plot,
    colour = competition_context
  )
) +
  ## Raw observations: use the same treatment visual language as model estimates.
  ## No facilitator = open symbol; high facilitator category = filled symbol.
  ## Colour follows neighbour-context colours used for the model estimates.
  geom_line(
    aes(group = segment_id, colour = competition_context),
    linewidth = 0.48,
    alpha = 0.85
  ) +
  geom_point(
    data = fecundity_pred_plot %>% filter(fac_level == "No facilitator"),
    aes(colour = competition_context),
    shape = 21,
    fill = "white",
    size = 1.5,
    stroke = 0.45,
    alpha = 0.80
  ) +
  geom_point(
    data = fecundity_pred_plot %>% filter(fac_level == "High facilitator"),
    aes(colour = competition_context, fill = competition_context),
    shape = 21,
    size = 1.8,
    stroke = 0.45,
    alpha = 0.80
  ) +
  ## Coloured vertical intervals, no horizontal caps, drawn in front of points
  geom_linerange(
    data = fecundity_pred_plot,
    aes(
      x = x_num,
      ymin = pred_low_plot,
      ymax = pred_high_plot,
      colour = competition_context
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    alpha = 1
  ) +
  scale_fill_manual(
    values = competition_cols,
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = tibble(
      panel_label = factor(
        rep("'(a)'~italic('T. cyanopetala')~'Sun'", 2),
        levels = levels(fecundity_pred_plot$panel_label)
      ),
      x_num = 1,
      pred_plot = 1,
      fac_level = factor(c("No facilitator", "High facilitator"), levels = fac_labs)
    ),
    inherit.aes = FALSE,
    aes(x = x_num, y = pred_plot, fill = fac_level),
    shape = 21,
    size = 1.8,
    stroke = 0.45,
    colour = "black",
    alpha = 0,
    show.legend = TRUE
  ) +
  facet_wrap(
    vars(panel_label),
    ncol = 1,
    scales = "free_y",
    labeller = label_parsed,
    strip.position = "top"
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      panel_label %in% c(
        "'(b)'~italic('T. cyanopetala')~'Shade'",
        "'(d)'~italic('T. ornata')~'Shade'"
      ) ~ scale_y_continuous(
        limits = c(0, 60),
        breaks = seq(0, 60, by = 20),
        minor_breaks = seq(10, 60, by = 10),
        expand = expansion(add = c(3, 2))
      ),
      TRUE ~ scale_y_continuous(
        limits = c(0, 125),
        breaks = seq(0, 120, by = 20),
        minor_breaks = seq(10, 120, by = 10),
        expand = expansion(add = c(3, 2))
      )
    )
  ) +
  scale_colour_manual(
    values = competition_cols,
    name = "Neighbour context"
  ) +
  scale_fill_manual(
    name = "Facilitator level",
    values = c(
      "No facilitator" = "white",
      "High facilitator" = "grey25"
    ),
    labels = c(
      "No facilitator" = "- facilitator",
      "High facilitator" = "+ facilitator"
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(SITE_LEVELS),
    labels = SITE_LEVELS,
    limits = c(0.50, length(SITE_LEVELS) + 0.50),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Site",
    y = "Predicted seed number"
  ) +
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY, size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.2, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey82", linewidth = 0.25, linetype = "dashed"),
    panel.grid.minor.y = element_line(colour = "grey88", linewidth = 0.20, linetype = "dashed"),
    panel.border = element_rect(linewidth = 0.35, colour = "black", fill = NA),
    panel.spacing.y = unit(0.35, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", hjust = 0),
    strip.placement = "outside",
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "top",
    legend.justification = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.50, "lines"),
    legend.key.width = unit(0.75, "lines"),
    legend.spacing.y = unit(4, "pt"),
    legend.box.spacing = unit(4, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      ncol = 1,
      override.aes = list(
        shape = 21,
        size = 2.3,
        stroke = 0.45,
        linewidth = 0,
        linetype = 0,
        alpha = 1,
        fill = unname(competition_cols)
      )
    ),
    fill = guide_legend(
      order = 2,
      ncol = 1,
      override.aes = list(
        shape = 21,
        size = 2.3,
        stroke = 0.45,
        colour = "black",
        alpha = 1
      )
    )
  )

############################################################
## Add right-side legend plus final-model effect annotation
############################################################

fmt_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("not retained")
  if (p < 0.001) return("< 0.001")
  paste0("= ", formatC(p, format = "f", digits = 3))
}

get_cond_p <- function(model, possible_terms) {
  coefs <- summary(model)$coefficients$cond
  term <- intersect(possible_terms, rownames(coefs))
  if (length(term) == 0) return(NA_real_)
  unname(coefs[term[1], "Pr(>|z|)"])
}

get_fac_interaction_terms <- function(model) {
  rn <- rownames(summary(model)$coefficients$cond)
  rn[grepl("fac_x", rn) & grepl(":", rn)]
}

get_anova_p <- function(model, term) {
  aa <- tryCatch(
    car::Anova(model, type = 3),
    error = function(e) NULL
  )
  
  if (is.null(aa)) return(NA_real_)
  
  aa_df <- as.data.frame(aa)
  rn <- rownames(aa_df)
  
  if (!term %in% rn) return(NA_real_)
  
  p_col <- grep("Pr", names(aa_df), value = TRUE)[1]
  if (is.na(p_col)) return(NA_real_)
  
  unname(aa_df[term, p_col])
}

trcy_fac_int_terms <- get_fac_interaction_terms(models[["TRCY"]])
tror_fac_int_terms <- get_fac_interaction_terms(models[["TROR"]])

trcy_fac_int_text <- if (length(trcy_fac_int_terms) == 0) {
  "No facilitator interactions retained"
} else {
  paste("Facilitator interactions retained:", paste(trcy_fac_int_terms, collapse = ", "))
}

tror_env_fac_int_text <- if (!any(grepl("Site|Cov", tror_fac_int_terms))) {
  "No Site/Cover × Facilitator retained"
} else {
  paste("Environment × Facilitator retained:", paste(tror_fac_int_terms[grepl("Site|Cov", tror_fac_int_terms)], collapse = ", "))
}

figure2_model_summary_text <- paste0(
  "TRCY final model\n",
  "Site P ", fmt_p(get_anova_p(models[["TRCY"]], "Site")),
  "; Cover P ", fmt_p(get_anova_p(models[["TRCY"]], "Cov")), "\n",
  "Facilitator P ", fmt_p(get_cond_p(models[["TRCY"]], "fac_x")), "\n",
  "Conspecific P ", fmt_p(get_cond_p(models[["TRCY"]], "cons_x")), "\n",
  "Heterospecific P ", fmt_p(get_cond_p(models[["TRCY"]], "hetero_x")), "\n",
  trcy_fac_int_text, "\n\n",
  "TROR final model\n",
  "Site P ", fmt_p(get_anova_p(models[["TROR"]], "Site")),
  "; Cover P ", fmt_p(get_anova_p(models[["TROR"]], "Cov")), "\n",
  "Site × Cover P ", fmt_p(get_anova_p(models[["TROR"]], "Site:Cov")), "\n",
  "Facilitator P ", fmt_p(get_cond_p(models[["TROR"]], "fac_x")), "\n",
  "Fac × conspecific P ", fmt_p(get_cond_p(models[["TROR"]], c("cons_x:fac_x", "fac_x:cons_x"))), "\n",
  "Fac × heterospecific P ", fmt_p(get_cond_p(models[["TROR"]], c("hetero_x:fac_x", "fac_x:hetero_x"))), "\n",
  tror_env_fac_int_text
)

legend_neighbour_df <- tibble(
  competition_context = factor(context_levels, levels = context_levels),
  label = context_levels,
  x_point = 0.08,
  x_text = 0.17,
  y = seq(0.94, 0.74, length.out = length(context_levels))
)

legend_facilitator_df <- tibble(
  fac_level = factor(c("No facilitator", "High facilitator"), levels = fac_labs),
  label = c("- facilitator", "+ facilitator"),
  x_point = 0.08,
  x_text = 0.17,
  y = c(0.595, 0.555)
)

fig2_right_panel <- ggplot() +
  annotate(
    "text",
    x = 0.00,
    y = 0.995,
    label = "Neighbour context",
    hjust = 0,
    vjust = 1,
    family = BASE_FAMILY,
    fontface = "bold",
    size = 10 / .pt
  ) +
  geom_point(
    data = legend_neighbour_df,
    aes(x = x_point, y = y, fill = competition_context),
    shape = 21,
    colour = "black",
    size = 2.3,
    stroke = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = legend_neighbour_df,
    aes(x = x_text, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    family = BASE_FAMILY,
    size = 10 / .pt
  ) +
  scale_fill_manual(values = competition_cols, guide = "none") +
  annotate(
    "text",
    x = 0.00,
    y = 0.67,
    label = "Facilitator level",
    hjust = 0,
    vjust = 1,
    family = BASE_FAMILY,
    fontface = "bold",
    size = 10 / .pt
  ) +
  geom_point(
    data = legend_facilitator_df %>% filter(fac_level == "No facilitator"),
    aes(x = x_point, y = y),
    shape = 21,
    fill = "white",
    colour = "black",
    size = 2.3,
    stroke = 0.45,
    show.legend = FALSE
  ) +
  geom_point(
    data = legend_facilitator_df %>% filter(fac_level == "High facilitator"),
    aes(x = x_point, y = y),
    shape = 21,
    fill = "grey25",
    colour = "black",
    size = 2.3,
    stroke = 0.45,
    show.legend = FALSE
  ) +
  geom_text(
    data = legend_facilitator_df,
    aes(x = x_text, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    family = BASE_FAMILY,
    size = 10 / .pt
  ) +
  annotate(
    "text",
    x = 0.00,
    y = 0.46,
    label = "Model effects",
    hjust = 0,
    vjust = 1,
    family = BASE_FAMILY,
    fontface = "bold",
    size = 10 / .pt
  ) +
  annotate(
    "text",
    x = 0.00,
    y = 0.425,
    label = figure2_model_summary_text,
    hjust = 0,
    vjust = 1,
    family = BASE_FAMILY,
    size = 8.0 / .pt,
    lineheight = 0.88
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = BASE_FAMILY) +
  theme(
    plot.margin = margin(2, 2, 2, 2)
  )

fig2_predicted_fecundity <-
  (fig2_predicted_fecundity + theme(legend.position = "none")) +
  fig2_right_panel +
  patchwork::plot_layout(widths = c(1, 0.34))

print(fig2_predicted_fecundity)

ggsave(
  file.path(FIGDIR, "figure2_predicted_fecundity_all_contrasts.png"),
  fig2_predicted_fecundity,
  width = FIG_WIDTH * 1.12 + (1.5 / 2.54),
  height = FIG_HEIGHT * 1.15,
  dpi = FIG_DPI
)


############################################################
## Figure 2 alternative inset-legend version removed
## Current preferred Figure 2 uses four rows and right-side legends.
############################################################

############################################################
## 7. FIGURE 3: SELECTED-MODEL FACILITATION-STRENGTH SUMMARY
############################################################

make_selected_facilitation_grid <- function(species_code) {
  d_sp <- get_species_data(species_code)
  meta <- get_species_meta(species_code)
  scalers <- get_scalers(d_sp, meta$density_form[[1]])
  
  fac_low <- 0
  fac_high <- q_positive(d_sp$facilitator_density, FAC_HIGH_QUANTILE)
  
  if (species_code == "TRCY") {
    grid_base <- tibble(
      species = species_code,
      species_label = species_labs_chr[[species_code]],
      effect_label = "Overall",
      Site = factor(SITE_LEVELS[1], levels = SITE_LEVELS),
      Cov = factor(COV_LEVELS[1], levels = COV_LEVELS),
      competition_context = factor("Overall", levels = c("Overall", context_levels)),
      conspecific_density = 0,
      heterospecific_density = 0
    )
  } else {
    grid_base <- make_context_table_global(d_sp) %>%
      mutate(
        species = species_code,
        species_label = species_labs_chr[[species_code]],
        effect_label = as.character(competition_context),
        Site = factor(SITE_LEVELS[1], levels = SITE_LEVELS),
        Cov = factor(COV_LEVELS[1], levels = COV_LEVELS)
      )
  }
  
  nd <- bind_rows(
    grid_base %>%
      mutate(
        fac_level = "No facilitator",
        facilitator_density = fac_low
      ),
    grid_base %>%
      mutate(
        fac_level = "High facilitator",
        facilitator_density = fac_high
      )
  ) %>%
    mutate(
      focsp = factor(species_code, levels = SPECIES_TO_PLOT),
      fac_level = factor(fac_level, levels = fac_labs),
      pair_id = interaction(species, effect_label, drop = TRUE),
      BLK = factor(levels(d_sp$BLK)[1], levels = levels(d_sp$BLK)),
      SiteBlock = factor(levels(d_sp$SiteBlock)[1], levels = levels(d_sp$SiteBlock))
    ) %>%
    add_scaled_predictors(scalers)
  
  nd
}

selected_fac_nd <- purrr::map_dfr(
  SPECIES_TO_PLOT,
  make_selected_facilitation_grid
)

selected_fac_summary <- purrr::map_dfr(SPECIES_TO_PLOT, function(sp) {
  nd_sp <- selected_fac_nd %>%
    filter(species == sp) %>%
    arrange(pair_id, fac_level)
  
  pred_draws <- draw_fixed_response(
    models[[sp]],
    nd_sp,
    n_draws = N_PARAM_DRAWS
  )
  
  low_i <- nd_sp$fac_level == "No facilitator"
  high_i <- nd_sp$fac_level == "High facilitator"
  
  low_nd <- nd_sp[low_i, ] %>% arrange(pair_id)
  high_nd <- nd_sp[high_i, ] %>% arrange(pair_id)
  
  low_draws <- pred_draws[low_i, , drop = FALSE]
  high_draws <- pred_draws[high_i, , drop = FALSE]
  
  low_draws <- low_draws[order(nd_sp$pair_id[low_i]), , drop = FALSE]
  high_draws <- high_draws[order(nd_sp$pair_id[high_i]), , drop = FALSE]
  
  eps <- 1e-12
  logRR_draws <- log(pmax(high_draws, eps) / pmax(low_draws, eps))
  
  out <- low_nd %>%
    transmute(
      species,
      species_label,
      effect_label,
      competition_context,
      pair_id,
      conspecific_density,
      heterospecific_density,
      fac_low = facilitator_density
    ) %>%
    left_join(
      high_nd %>% transmute(pair_id, fac_high = facilitator_density),
      by = "pair_id"
    )
  
  bind_cols(
    out,
    tibble(
      logRR = apply(logRR_draws, 1, median, na.rm = TRUE),
      logRR_low = apply(logRR_draws, 1, quantile, probs = 0.025, na.rm = TRUE),
      logRR_high = apply(logRR_draws, 1, quantile, probs = 0.975, na.rm = TRUE),
      p_fac_positive = rowMeans(logRR_draws > 0, na.rm = TRUE)
    )
  )
}) %>%
  mutate(
    species = factor(species, levels = SPECIES_TO_PLOT),
    species_label = factor(species_label, levels = species_labs_chr[SPECIES_TO_PLOT]),
    effect_label = factor(
      effect_label,
      levels = c("Overall", context_levels)
    ),
    pct_change = 100 * (exp(logRR) - 1),
    pct_low = 100 * (exp(logRR_low) - 1),
    pct_high = 100 * (exp(logRR_high) - 1)
  )

# write_csv(
#   selected_fac_summary,
#   file.path(FIGDIR, "figure3_selected_facilitation_strength_summary.csv")
# )


############################################################
## Plot Figure 2: four-row version as bars + error bars
## Rows:
##   (a) TRCY Sun
##   (b) TRCY Shade
##   (c) TROR Sun
##   (d) TROR Shade
############################################################

BAR_WIDTH <- 0.075
ERRORBAR_WIDTH <- 0.045
BAR_LINEWIDTH <- 0.32
ERRORBAR_LINEWIDTH <- 0.35

fig2_predicted_fecundity <- ggplot(
  fecundity_pred_plot,
  aes(x = x_num)
) +
  ## No facilitator: open bars
  geom_col(
    data = fecundity_pred_plot %>% filter(fac_level == "No facilitator"),
    aes(
      y = pred_plot,
      colour = competition_context
    ),
    fill = "white",
    width = BAR_WIDTH,
    linewidth = BAR_LINEWIDTH,
    alpha = 1
  ) +
  ## High facilitator: filled bars
  geom_col(
    data = fecundity_pred_plot %>% filter(fac_level == "High facilitator"),
    aes(
      y = pred_plot,
      fill = competition_context,
      colour = competition_context
    ),
    width = BAR_WIDTH,
    linewidth = BAR_LINEWIDTH,
    alpha = 1
  ) +
  ## Error bars for both facilitator levels
  geom_errorbar(
    aes(
      ymin = pred_low_plot,
      ymax = pred_high_plot,
      colour = competition_context
    ),
    width = ERRORBAR_WIDTH,
    linewidth = ERRORBAR_LINEWIDTH,
    alpha = 1
  ) +
  ## Raw observations: use the same treatment visual language as model estimates.
  ## No facilitator = open symbol; high facilitator category = filled symbol.
  ## Colour follows neighbour-context colours used for the model estimates.
  geom_point(
    data = raw_fecundity_context_plot %>%
      dplyr::filter(fac_level == "No facilitator"),
    inherit.aes = FALSE,
    aes(
      x = x_num_raw,
      y = fecundity,
      colour = competition_context
    ),
    shape = 21,
    fill = "white",
    alpha = RAW_OVERLAY_ALPHA,
    size = RAW_OVERLAY_SIZE,
    stroke = 0.18,
    show.legend = FALSE
  ) +
  geom_point(
    data = raw_fecundity_context_plot %>%
      dplyr::filter(fac_level == "High facilitator"),
    inherit.aes = FALSE,
    aes(
      x = x_num_raw,
      y = fecundity,
      colour = competition_context,
      fill = competition_context
    ),
    shape = 21,
    alpha = RAW_OVERLAY_ALPHA,
    size = RAW_OVERLAY_SIZE,
    stroke = 0.18,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = competition_cols,
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  ## Invisible dummy layer for facilitator-level legend
  geom_col(
    data = tibble(
      panel_label = factor(
        rep("'(a)'~italic('T. cyanopetala')~'Sun'", 2),
        levels = levels(fecundity_pred_plot$panel_label)
      ),
      x_num = 1,
      pred_plot = 1,
      fac_level = factor(c("No facilitator", "High facilitator"), levels = fac_labs)
    ),
    inherit.aes = FALSE,
    aes(x = x_num, y = pred_plot, fill = fac_level),
    width = BAR_WIDTH,
    colour = "black",
    linewidth = BAR_LINEWIDTH,
    alpha = 0,
    show.legend = TRUE
  ) +
  facet_wrap(
    vars(panel_label),
    ncol = 1,
    scales = "fixed",
    labeller = label_parsed,
    strip.position = "top"
  ) +
  scale_colour_manual(
    values = competition_cols,
    name = "Neighbour context"
  ) +
  scale_fill_manual(
    name = "Facilitator level",
    values = c(
      "No facilitator" = "white",
      "High facilitator" = "grey25"
    ),
    labels = c(
      "No facilitator" = "- facilitator",
      "High facilitator" = "+ facilitator"
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(SITE_LEVELS),
    labels = SITE_LEVELS,
    limits = c(0.50, length(SITE_LEVELS) + 0.50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 125),
    breaks = seq(0, 120, by = 20),
    minor_breaks = seq(10, 120, by = 10),
    expand = expansion(add = c(3, 2))
  ) +
  labs(
    x = "Site",
    y = "Predicted seed number"
  ) +
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY, size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.2, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey82", linewidth = 0.25, linetype = "dashed"),
    panel.grid.minor.y = element_line(colour = "grey88", linewidth = 0.20, linetype = "dashed"),
    panel.border = element_rect(linewidth = 0.35, colour = "black", fill = NA),
    panel.spacing.y = unit(0.35, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", hjust = 0),
    strip.placement = "outside",
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "top",
    legend.justification = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.50, "lines"),
    legend.key.width = unit(0.75, "lines"),
    legend.spacing.y = unit(4, "pt"),
    legend.box.spacing = unit(4, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      ncol = 1,
      override.aes = list(
        fill = unname(competition_cols),
        colour = unname(competition_cols),
        linewidth = 0,
        alpha = 1
      )
    ),
    fill = guide_legend(
      order = 2,
      ncol = 1,
      override.aes = list(
        colour = "black",
        alpha = 1
      )
    )
  )

ggsave(
  file.path(FIGDIR, "figure2_predicted_fecundity_all_contrasts_Bars.png"),
  fig2_predicted_fecundity,
  width = FIG_WIDTH * 1.12 + (1.5 / 2.54),
  height = FIG_HEIGHT * 1.15,
  dpi = FIG_DPI
)


############################################################
## Plot Figure 3: selected-model facilitation-strength summary
## Single-panel version with species legend
############################################################

fig3_x_levels <- c(
  "TRCY overall",
  "No neighbours",
  "Low conspecific",
  "Low heterospecific",
  "High conspecific",
  "High heterospecific"
)

fig3_x_labels <- c(
  "TRCY overall" = "Overall",
  "No neighbours" = "No neighbours",
  "Low conspecific" = "Low conspecific",
  "Low heterospecific" = "Low heterospecific",
  "High conspecific" = "High conspecific",
  "High heterospecific" = "High heterospecific"
)

fig3_cols <- c(
  "TRCY overall" = "grey45",
  "No neighbours" = "black",
  "Low conspecific" = competition_cols[["Low conspecific"]],
  "Low heterospecific" = competition_cols[["Low heterospecific"]],
  "High conspecific" = competition_cols[["High conspecific"]],
  "High heterospecific" = competition_cols[["High heterospecific"]]
)

selected_fac_summary_plot <- selected_fac_summary %>%
  mutate(
    effect_plot = case_when(
      species == "TRCY" ~ "TRCY overall",
      species == "TROR" ~ as.character(effect_label),
      TRUE ~ as.character(effect_label)
    ),
    effect_plot = factor(effect_plot, levels = fig3_x_levels),
    x_num = as.numeric(effect_plot),
    species_group = factor(
      if_else(species == "TRCY", "TRCY", "TROR"),
      levels = c("TRCY", "TROR")
    )
  )

logrr_range <- range(
  c(selected_fac_summary_plot$logRR_low, selected_fac_summary_plot$logRR_high),
  na.rm = TRUE
)

logrr_limit <- max(abs(logrr_range))
logrr_limit <- max(0.5, ceiling(logrr_limit * 2) / 2)

fig3_selected_facilitation <- ggplot(
  selected_fac_summary_plot,
  aes(
    x = x_num,
    y = logRR,
    ymin = logRR_low,
    ymax = logRR_high,
    colour = effect_plot,
    fill = effect_plot,
    shape = species_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35 * 1.15
  ) +
  geom_vline(
    xintercept = 1.5,
    linewidth = 0.25,
    colour = "grey70"
  ) +
  geom_linerange(
    linewidth = 0.45 * 1.15,
    alpha = 1
  ) +
  geom_point(
    size = 2.2,
    stroke = 0.45,
    alpha = 1
  ) +
  scale_colour_manual(values = fig3_cols, guide = "none") +
  scale_fill_manual(values = fig3_cols, guide = "none") +
  scale_shape_manual(
    name = "Species",
    values = c("TRCY" = 23, "TROR" = 21),
    labels = c(
      "TRCY" = expression(italic("T. cyanopetala")),
      "TROR" = expression(italic("T. ornata"))
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(fig3_x_levels),
    labels = fig3_x_labels,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  coord_cartesian(ylim = c(-logrr_limit, logrr_limit)) +
  labs(
    x = NULL,
    y = "Facilitation effect"
  ) +
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY, size = 11),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 10, colour = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.2, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.35, colour = "black", fill = NA),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.45, "lines"),
    legend.key.width = unit(0.65, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -2),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        colour = c("grey45", "black"),
        fill = c("grey45", "black"),
        size = 2.4,
        stroke = 0.45
      )
    )
  )

print(fig3_selected_facilitation)

ggsave(
  file.path(FIGDIR, "figure3_selected_facilitation_strength.png"),
  fig3_selected_facilitation,
  width = 4,
  height = 3,
  dpi = FIG_DPI
)

############################################################
## 8. MODEL OUTPUT FILES FOR ANNOTATION DECISIONS
############################################################

sink(file.path(FIGDIR, "final_model_summaries_for_figure_annotations.txt"))
cat("Final model summaries used for figures\n")
cat("Contrast: no facilitator = NO_GRASS 0; high facilitator = species-specific 90th percentile of positive NO_GRASS values\n\n")

for (sp in SPECIES_TO_PLOT) {
  cat("\n============================================================\n")
  cat(sp, " final model\n")
  cat("============================================================\n\n")
  print(summary(models[[sp]]))
}

sink()

############################################################
## 9. QUICK DIAGNOSTIC PRINTS
############################################################

cat("\nFacilitator contrast values used:\n")
purrr::walk(SPECIES_TO_PLOT, function(sp) {
  d_sp <- get_species_data(sp)
  vals <- d_sp %>%
    summarise(
      no_facilitator = 0,
      high_facilitator_90pct_positive_values = q_positive(facilitator_density, FAC_HIGH_QUANTILE),
      max_facilitator = max(facilitator_density, na.rm = TRUE)
    )
  cat("\n", sp, ":\n", sep = "")
  print(vals)
})

cat("\nFiles written to: ", FIGDIR, "\n", sep = "")
