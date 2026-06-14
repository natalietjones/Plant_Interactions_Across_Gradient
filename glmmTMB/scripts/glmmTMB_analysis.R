############################################################
## FREQUENTIST GLMM
## 
##
## Purpose:
##   1. Check response and predictor distributions before model fitting.
##   2. Fit species-specific GLMMs for TRCY and TROR.
##   3. Select count family, density transformation, and fixed-effect structure by AIC.
##   4. Include targeted a priori facilitation interactions in the model-selection set.
##   5. Use manual likelihood-ratio tests for targeted nested facilitation contrasts.
##   6. Record non-estimable/data-deficient models.
##
## Main decisions:
##   - Fecundity is kept as a raw count response. It is NOT log-transformed
##     because the GLMM uses a count distribution with a log link.
##   - Density predictors are evaluated as:
##       raw_z   = z-score(raw density)
##       log1p_z = z-score(log1p(raw density))
##     The log1p transform reduces leverage from right-skewed density values.
##   - z-scaling is recalculated within each species-specific analysis dataset,
##     including outlier sensitivity datasets.
##   - Outlier removal is sensitivity analysis only, unless observations are
##     known data-entry or measurement errors.
##   - A priori facilitation interactions are tested using targeted models;
##     they are retained in the final model only if they are estimable and
##     improve model support from AIC
############################################################

############################################################
## 0. PACKAGES AND SETUP
############################################################

pkgs <- c(
  "tidyverse", "glmmTMB", "DHARMa", "car",
  "broom.mixed", "performance", "here", "readr",
  "emmeans", "multcomp", "multcompView"
)

missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

## Sum-to-zero contrasts are required for interpretable Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

############################################################
## 0a. OUTPUT OPTIONS
############################################################

## Keep CSV writing focused on the interpretation tables.
## Set TRUE if you want all extra diagnostic tables written as CSVs.
WRITE_EXTENDED_TABLES <- FALSE

## Set TRUE to save every species-level emmeans table.
## Combined posthoc interpretation tables are saved either way.
SAVE_DETAILED_POSTHOC <- FALSE

############################################################
## 0b. PROJECT PATHS
############################################################

## Expected folder layout:
##   <project root>/data/field/fullWAdata.csv   preferred
##   <project root>/data/fullWAdata.csv         fallback
##   <project root>/glmmTMB/scripts/this_script.R
##   <project root>/glmmTMB/glmm_outputs/
##   <project root>/glmmTMB/figures_final_glmm/
##
## This block is deliberately robust to running the script from either:
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

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, path)
  if (!file.exists(path)) stop("File was not created: ", path)
  message("Saved: ", normalizePath(path, winslash = "/", mustWork = TRUE))
  invisible(path)
}

safe_save_rds <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(x, path)
  if (!file.exists(path)) stop("File was not created: ", path)
  message("Saved: ", normalizePath(path, winslash = "/", mustWork = TRUE))
  invisible(path)
}

cat("\nPath setup:\n")
cat("  glmm_dir:     ", normalizePath(glmm_dir, winslash = "/"), "\n", sep = "")
cat("  project_root: ", normalizePath(project_root, winslash = "/"), "\n", sep = "")
cat("  data_file:    ", normalizePath(data_file, winslash = "/"), "\n", sep = "")
cat("  output dir:   ", normalizePath(outdir, winslash = "/"), "\n", sep = "")
cat("  figure dir:   ", normalizePath(figdir, winslash = "/"), "\n\n", sep = "")

decision_log <- character(0)

log_decision <- function(...) {
  msg <- paste0(...)
  decision_log <<- c(decision_log, msg)
  message(msg)
}

print_tbl <- function(x, n = Inf) {
  ## Avoid print.data.frame(..., n = Inf), where `n` can be partially
  ## matched to `na.print` and halt the script with:
  ## "invalid 'na.print' specification".
  if (inherits(x, "data.frame")) {
    print(tibble::as_tibble(x), n = n, width = Inf)
  } else {
    print(x)
  }
}

log_decision("DECISION: Fecundity is analysed as a raw count response in count GLMMs with a log link.")
log_decision("DECISION: Density predictors are compared as raw z-scores versus log1p-transformed z-scores.")
log_decision("DECISION: z-scaling is recalculated within each species-specific analysis dataset, including outlier sensitivity datasets.")
log_decision("DECISION: Outlier removal is treated as sensitivity analysis, not the primary analysis.")
log_decision("DECISION: Targeted a priori facilitation interactions are part of model selection and are retained only if estimable and supported.")
log_decision("DECISION: Random intercepts use Site x Cover x Block, so Sun block 1 and Shade block 1 are not treated as the same random-effect level.")

############################################################
## Robust print helper
############################################################

print_posthoc <- function(x, n = Inf) {
  if (inherits(x, "data.frame")) {
    print(tibble::as_tibble(x), n = n, width = Inf)
  } else {
    print(x)
  }
}

############################################################
## Robust emmeans output helpers
############################################################

get_ci_cols <- function(df) {
  lower_candidates <- c("asymp.LCL", "lower.CL", "LCL", "lower.HPD")
  upper_candidates <- c("asymp.UCL", "upper.CL", "UCL", "upper.HPD")
  
  lower_col <- intersect(lower_candidates, names(df))[1]
  upper_col <- intersect(upper_candidates, names(df))[1]
  
  if (is.na(lower_col) || is.na(upper_col)) {
    stop(
      "Could not find confidence interval columns. Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  list(lower = lower_col, upper = upper_col)
}

get_ratio_col <- function(df) {
  ratio_candidates <- c("ratio", "response", "emmean")
  ratio_col <- intersect(ratio_candidates, names(df))[1]
  
  if (is.na(ratio_col)) {
    stop(
      "Could not find ratio/response column. Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  ratio_col
}

as_percent_contrast <- function(x) {
  df <- tibble::as_tibble(as.data.frame(x))
  ci <- get_ci_cols(df)
  ratio_col <- get_ratio_col(df)
  
  df %>%
    dplyr::mutate(
      percent_change = (.data[[ratio_col]] - 1) * 100,
      lower_percent = (.data[[ci$lower]] - 1) * 100,
      upper_percent = (.data[[ci$upper]] - 1) * 100
    )
}

as_percent_trend <- function(x, trend_col) {
  df <- tibble::as_tibble(as.data.frame(x))
  ci <- get_ci_cols(df)
  
  if (!trend_col %in% names(df)) {
    stop(
      "Could not find trend column '", trend_col, "'. Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  df %>%
    dplyr::mutate(
      estimate_exp = exp(.data[[trend_col]]),
      lower_exp = exp(.data[[ci$lower]]),
      upper_exp = exp(.data[[ci$upper]]),
      percent_change_per_1SD = (estimate_exp - 1) * 100,
      lower_percent = (lower_exp - 1) * 100,
      upper_percent = (upper_exp - 1) * 100
    )
}

############################################################
## Manual compact-letter display from pairwise p-values
############################################################

make_letters_from_pairs <- function(emm_resp, pair_resp, group_col = "Site", alpha = 0.05) {
  
  emm_df <- tibble::as_tibble(as.data.frame(emm_resp))
  pair_df <- tibble::as_tibble(as.data.frame(pair_resp))
  
  if (!group_col %in% names(emm_df)) {
    stop("Group column not found in emmeans object: ", group_col)
  }
  
  contrast_col <- intersect(c("contrast", "comparison"), names(pair_df))[1]
  p_col <- intersect(c("p.value", "p.value."), names(pair_df))[1]
  
  if (is.na(contrast_col) || is.na(p_col)) {
    warning(
      "Could not compute letters. Pairwise table columns were:\n",
      paste(names(pair_df), collapse = ", ")
    )
    return(
      emm_df %>%
        dplyr::mutate(.group = NA_character_)
    )
  }
  
  contrasts_clean <- pair_df[[contrast_col]] %>%
    as.character() %>%
    stringr::str_replace_all("\\s+", "") %>%
    stringr::str_replace_all("/", "-") %>%
    stringr::str_replace_all("−", "-")
  
  pvals <- pair_df[[p_col]]
  names(pvals) <- contrasts_clean
  
  letter_obj <- multcompView::multcompLetters(
    pvals,
    threshold = alpha
  )
  
  letters_df <- tibble::tibble(
    !!group_col := names(letter_obj$Letters),
    .group = unname(letter_obj$Letters)
  )
  
  emm_df %>%
    dplyr::left_join(letters_df, by = group_col) %>%
    dplyr::mutate(.group = stringr::str_squish(.group))
}


############################################################
## Core modelling helpers
############################################################

z_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

skewness_simple <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3)
}

add_density_transforms <- function(dat) {
  dat %>%
    dplyr::mutate(
      cons_log1p = log1p(conspecific_density),
      hetero_log1p = log1p(heterospecific_density),
      fac_log1p = log1p(facilitator_density)
    ) %>%
    dplyr::group_by(focsp) %>%
    dplyr::mutate(
      cons_raw_z = z_safe(conspecific_density),
      hetero_raw_z = z_safe(heterospecific_density),
      fac_raw_z = z_safe(facilitator_density),
      cons_log1p_z = z_safe(cons_log1p),
      hetero_log1p_z = z_safe(hetero_log1p),
      fac_log1p_z = z_safe(fac_log1p)
    ) %>%
    dplyr::ungroup()
}

assign_density_form <- function(dat, density_form = c("raw", "log1p")) {
  density_form <- match.arg(density_form)
  
  if (density_form == "raw") {
    dat %>%
      dplyr::mutate(
        cons_x = cons_raw_z,
        hetero_x = hetero_raw_z,
        fac_x = fac_raw_z
      )
  } else {
    dat %>%
      dplyr::mutate(
        cons_x = cons_log1p_z,
        hetero_x = hetero_log1p_z,
        fac_x = fac_log1p_z
      )
  }
}

## Random structure deliberately uses Site x Cover x Block.
## This means Sun block 1 and Shade block 1 are NOT linked by the random effect.
make_formula <- function(rhs) {
  stats::as.formula(
    paste0("fecundity ~ ", rhs, " + (1 | SiteCovBlock)")
  )
}

ctrl <- glmmTMB::glmmTMBControl(
  optimizer = stats::optim,
  optArgs = list(method = "BFGS"),
  profile = TRUE
)

family_specs <- list(
  poisson_count = list(
    family = stats::poisson(link = "log"),
    ziformula = ~0
  ),
  nb1_count = list(
    family = glmmTMB::nbinom1(link = "log"),
    ziformula = ~0
  ),
  nb2_count = list(
    family = glmmTMB::nbinom2(link = "log"),
    ziformula = ~0
  ),
  zip_count = list(
    family = stats::poisson(link = "log"),
    ziformula = ~1
  ),
  zinb1_count = list(
    family = glmmTMB::nbinom1(link = "log"),
    ziformula = ~1
  ),
  zinb2_count = list(
    family = glmmTMB::nbinom2(link = "log"),
    ziformula = ~1
  )
)

safe_fit <- function(formula, data, family, ziformula = ~0, control = ctrl) {
  tryCatch(
    suppressWarnings(
      glmmTMB::glmmTMB(
        formula = formula,
        data = data,
        family = family,
        ziformula = ziformula,
        control = control
      )
    ),
    error = function(e) {
      paste0("FIT_ERROR: ", conditionMessage(e))
    }
  )
}

model_row <- function(model, model_name, n_obs) {
  
  base_row <- tibble::tibble(
    model = model_name,
    usable = FALSE,
    n = n_obs,
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    df = NA_real_,
    dropped_reason = NA_character_,
    nonestimable_terms = NA_character_,
    n_NA_fixed = NA_integer_,
    n_NA_vcov = NA_integer_
  )
  
  if (!inherits(model, "glmmTMB")) {
    base_row$dropped_reason <- as.character(model)[1]
    return(base_row)
  }
  
  fixed <- tryCatch(glmmTMB::fixef(model)$cond, error = function(e) NA)
  vc <- tryCatch(as.matrix(stats::vcov(model)$cond), error = function(e) matrix(NA_real_, 1, 1))
  ll_obj <- tryCatch(stats::logLik(model), error = function(e) NA)
  
  n_na_fixed <- sum(is.na(fixed))
  n_na_vcov <- sum(is.na(vc))
  nonestimable_terms <- if (n_na_fixed > 0) {
    paste(names(fixed)[is.na(fixed)], collapse = "; ")
  } else {
    NA_character_
  }
  
  convergence_ok <- isTRUE(model$fit$convergence == 0)
  pdhess_ok <- isTRUE(model$sdr$pdHess)
  ll <- suppressWarnings(as.numeric(ll_obj))
  aic_val <- suppressWarnings(tryCatch(stats::AIC(model), error = function(e) NA_real_))
  bic_val <- suppressWarnings(tryCatch(stats::BIC(model), error = function(e) NA_real_))
  df_val <- suppressWarnings(tryCatch(attr(ll_obj, "df"), error = function(e) NA_real_))
  
  reason <- character(0)
  if (!convergence_ok) reason <- c(reason, paste0("nonzero convergence code: ", model$fit$convergence))
  if (!pdhess_ok) reason <- c(reason, "non-positive-definite Hessian")
  if (n_na_fixed > 0) reason <- c(reason, "NA/non-estimable fixed effects")
  if (n_na_vcov > 0) reason <- c(reason, "NA fixed-effect covariance entries")
  if (!is.finite(ll)) reason <- c(reason, "non-finite logLik")
  if (!is.finite(aic_val)) reason <- c(reason, "non-finite AIC")
  
  usable <- length(reason) == 0
  
  tibble::tibble(
    model = model_name,
    usable = usable,
    n = n_obs,
    logLik = ll,
    AIC = aic_val,
    BIC = bic_val,
    df = as.numeric(df_val),
    dropped_reason = ifelse(usable, NA_character_, paste(reason, collapse = "; ")),
    nonestimable_terms = nonestimable_terms,
    n_NA_fixed = n_na_fixed,
    n_NA_vcov = n_na_vcov
  )
}

aic_table <- function(models, n_obs) {
  out <- purrr::imap_dfr(
    models,
    ~ model_row(.x, .y, n_obs)
  )
  
  if (any(out$usable, na.rm = TRUE)) {
    min_aic <- min(out$AIC[out$usable], na.rm = TRUE)
    out <- out %>%
      dplyr::mutate(
        delta_AIC = dplyr::if_else(
          usable,
          AIC - min_aic,
          NA_real_
        )
      )
  } else {
    out <- out %>% dplyr::mutate(delta_AIC = NA_real_)
  }
  
  out %>%
    dplyr::arrange(delta_AIC, AIC, df)
}

choose_best <- function(models, tab, delta_rule = 2) {
  usable <- tab %>%
    dplyr::filter(usable, is.finite(AIC), is.finite(delta_AIC))
  
  if (nrow(usable) == 0) {
    stop("No usable models available for selection.")
  }
  
  ## Parsimony rule: among models within delta_rule AIC of the best model,
  ## choose the one with the fewest parameters, then lowest AIC.
  selected <- usable %>%
    dplyr::filter(delta_AIC <= delta_rule) %>%
    dplyr::arrange(df, AIC) %>%
    dplyr::slice(1)
  
  selected_name <- selected$model[[1]]
  list(
    name = selected_name,
    model = models[[selected_name]],
    selected_table = selected
  )
}

manual_lrt_row <- function(models, base_name, test_name, species, n_obs) {
  
  base_model <- models[[base_name]]
  test_model <- models[[test_name]]
  
  base_row <- model_row(base_model, base_name, n_obs)
  test_row <- model_row(test_model, test_name, n_obs)
  
  if (!isTRUE(base_row$usable) || !isTRUE(test_row$usable)) {
    return(tibble::tibble(
      species = species,
      base_model = base_name,
      test_model = test_name,
      base_usable = base_row$usable,
      test_usable = test_row$usable,
      base_AIC = base_row$AIC,
      test_AIC = test_row$AIC,
      delta_AIC_test_minus_base = test_row$AIC - base_row$AIC,
      LRT = NA_real_,
      df_diff = NA_real_,
      p_value = NA_real_,
      note = paste(
        "LRT not run:",
        paste(na.omit(c(base_row$dropped_reason, test_row$dropped_reason)), collapse = "; ")
      )
    ))
  }
  
  df_diff <- test_row$df - base_row$df
  LRT <- 2 * (test_row$logLik - base_row$logLik)
  
  if (!is.finite(df_diff) || df_diff <= 0 || !is.finite(LRT)) {
    p_value <- NA_real_
    note <- "LRT not valid because df difference or likelihood difference was not positive/finite."
  } else {
    p_value <- stats::pchisq(LRT, df = df_diff, lower.tail = FALSE)
    note <- NA_character_
  }
  
  tibble::tibble(
    species = species,
    base_model = base_name,
    test_model = test_name,
    base_usable = TRUE,
    test_usable = TRUE,
    base_AIC = base_row$AIC,
    test_AIC = test_row$AIC,
    delta_AIC_test_minus_base = test_row$AIC - base_row$AIC,
    LRT = LRT,
    df_diff = df_diff,
    p_value = p_value,
    note = note
  )
}

extract_fixed_effects <- function(model) {
  coefs <- summary(model)$coefficients$cond
  
  tibble::as_tibble(coefs, rownames = "term") %>%
    dplyr::rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      z_value = `z value`,
      p_value = `Pr(>|z|)`
    )
}


############################################################
## 2. READ AND PREPARE DATA
############################################################
raw <- readr::read_csv(data_file, show_col_types = FALSE) %>%
  dplyr::mutate(obs_id = dplyr::row_number())

required_cols <- c(
  "Site", "Cov", "BLK", "TROR", "TRCY", "NO_GRASS",
  "FOCAL_GERM", "focsp", "fitness"
)

missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_glmm <- raw %>%
  dplyr::mutate(
    FOCAL_GERM_BIN = dplyr::case_when(
      FOCAL_GERM %in% c("Y", "Yes", "YES", "y", "yes", 1, "1") ~ 1,
      FOCAL_GERM %in% c("N", "No", "NO", "n", "no", 0, "0") ~ 0,
      TRUE ~ NA_real_
    ),
    Site = dplyr::case_when(
      Site %in% c("BEN", "Ben", "ben") ~ "Ben",
      Site %in% c("GH", "Gh", "gh") ~ "GH",
      Site %in% c("NAM", "Nam", "nam") ~ "Nam",
      Site %in% c("PJ", "Pj", "pj") ~ "PJ",
      Site %in% c("CS", "Cs", "cs") ~ "CS",
      TRUE ~ as.character(Site)
    ),
    Site = factor(Site, levels = c("Ben", "GH", "Nam", "PJ", "CS")),
    Cov = dplyr::case_when(
      Cov %in% c("SUN", "Sun", "sun") ~ "Sun",
      Cov %in% c("SH", "SHADE", "Shade", "shade", "sh") ~ "Shade",
      TRUE ~ as.character(Cov)
    ),
    Cov = factor(Cov, levels = c("Sun", "Shade")),
    focsp = factor(focsp, levels = c("TRCY", "TROR")),
    BLK = factor(BLK),
    fecundity = as.numeric(fitness),
    conspecific_density = dplyr::case_when(
      focsp == "TRCY" ~ as.numeric(TRCY),
      focsp == "TROR" ~ as.numeric(TROR),
      TRUE ~ NA_real_
    ),
    heterospecific_density = dplyr::case_when(
      focsp == "TRCY" ~ as.numeric(TROR),
      focsp == "TROR" ~ as.numeric(TRCY),
      TRUE ~ NA_real_
    ),
    facilitator_density = as.numeric(NO_GRASS),
    SiteCovBlock = interaction(Site, Cov, BLK, drop = TRUE, sep = "_")
  ) %>%
  dplyr::filter(
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
    !is.na(SiteCovBlock)
  ) %>%
  dplyr::mutate(zero_fecundity = fecundity == 0) %>%
  droplevels() %>%
  add_density_transforms()

############################################################
## 3. PRE-MODEL DIAGNOSTICS: RESPONSE, DENSITIES, OUTLIERS
############################################################

response_summary <- df_glmm %>%
  dplyr::group_by(focsp) %>%
  dplyr::summarise(
    n = n(),
    zeros = sum(fecundity == 0),
    zero_prop = mean(fecundity == 0),
    mean_fecundity = mean(fecundity),
    var_fecundity = var(fecundity),
    var_to_mean = var_fecundity / mean_fecundity,
    median_fecundity = median(fecundity),
    q95_fecundity = quantile(fecundity, 0.95),
    q99_fecundity = quantile(fecundity, 0.99),
    max_fecundity = max(fecundity),
    skew_fecundity_raw = skewness_simple(fecundity),
    skew_fecundity_log1p = skewness_simple(log1p(fecundity)),
    .groups = "drop"
  )

density_transform_summary <- df_glmm %>%
  dplyr::group_by(focsp) %>%
  dplyr::summarise(
    n = n(),
    cons_raw_skew = skewness_simple(conspecific_density),
    cons_log1p_skew = skewness_simple(cons_log1p),
    hetero_raw_skew = skewness_simple(heterospecific_density),
    hetero_log1p_skew = skewness_simple(hetero_log1p),
    fac_raw_skew = skewness_simple(facilitator_density),
    fac_log1p_skew = skewness_simple(fac_log1p),
    .groups = "drop"
  )

z_check <- df_glmm %>%
  dplyr::group_by(focsp) %>%
  dplyr::summarise(
    cons_raw_z_mean = mean(cons_raw_z),
    cons_raw_z_sd = sd(cons_raw_z),
    cons_log1p_z_mean = mean(cons_log1p_z),
    cons_log1p_z_sd = sd(cons_log1p_z),
    hetero_raw_z_mean = mean(hetero_raw_z),
    hetero_raw_z_sd = sd(hetero_raw_z),
    hetero_log1p_z_mean = mean(hetero_log1p_z),
    hetero_log1p_z_sd = sd(hetero_log1p_z),
    fac_raw_z_mean = mean(fac_raw_z),
    fac_raw_z_sd = sd(fac_raw_z),
    fac_log1p_z_mean = mean(fac_log1p_z),
    fac_log1p_z_sd = sd(fac_log1p_z),
    .groups = "drop"
  )

data_support_by_cell <- df_glmm %>%
  dplyr::group_by(focsp, Site, Cov) %>%
  dplyr::summarise(
    n = n(),
    fecundity_max = max(fecundity),
    cons_nonzero = sum(conspecific_density > 0),
    hetero_nonzero = sum(heterospecific_density > 0),
    fac_nonzero = sum(facilitator_density > 0),
    cons_sd = sd(conspecific_density),
    hetero_sd = sd(heterospecific_density),
    fac_sd = sd(facilitator_density),
    flag_low_n = n < 10,
    flag_no_cons_variation = is.na(cons_sd) | cons_sd == 0,
    flag_no_hetero_variation = is.na(hetero_sd) | hetero_sd == 0,
    flag_no_fac_variation = is.na(fac_sd) | fac_sd == 0,
    .groups = "drop"
  )

outlier_flags <- df_glmm %>%
  dplyr::group_by(focsp) %>%
  dplyr::mutate(
    fec_q1 = quantile(fecundity, 0.25, na.rm = TRUE),
    fec_q3 = quantile(fecundity, 0.75, na.rm = TRUE),
    fec_iqr = fec_q3 - fec_q1,
    fec_upper_iqr3 = fec_q3 + 3 * fec_iqr,
    fec_q99 = quantile(fecundity, 0.99, na.rm = TRUE),
    high_fec_iqr3 = fecundity > fec_upper_iqr3,
    high_fec_q99 = fecundity >= fec_q99,
    high_fec_flag = high_fec_iqr3 | high_fec_q99
  ) %>%
  dplyr::ungroup()

high_fec_report <- outlier_flags %>%
  dplyr::filter(high_fec_flag) %>%
  dplyr::arrange(focsp, desc(fecundity)) %>%
  dplyr::select(
    obs_id, focsp, Site, Cov, BLK,
    fecundity, conspecific_density, heterospecific_density, facilitator_density,
    TRCY, TROR, NO_GRASS, high_fec_iqr3, high_fec_q99
  )

log_decision("\n===== RESPONSE SUMMARY =====")
print_tbl(response_summary)
log_decision("\n===== DENSITY TRANSFORMATION SUMMARY =====")
print_tbl(density_transform_summary)
log_decision("\n===== Z-SCORE CHECK: values should have mean ~0 and SD ~1 within species =====")
print_tbl(z_check)
log_decision("\n===== SITE x COVER DATA SUPPORT CHECK =====")
print_tbl(data_support_by_cell, n = Inf)
log_decision("\n===== HIGH-FECUNDITY OBSERVATIONS FLAGGED FOR SENSITIVITY ONLY =====")
print_tbl(high_fec_report, n = Inf)

############################################################
## 4. DATA SUPPORT FOR SITE x COVER x FACILITATION INTERACTIONS
############################################################

interaction_support_from_data <- function(dat, species_code) {
  dat %>%
    dplyr::filter(focsp == species_code) %>%
    dplyr::group_by(Site, Cov) %>%
    dplyr::summarise(
      species = species_code,
      n = n(),
      fecundity_zeros = sum(fecundity == 0),
      fecundity_zero_prop = mean(fecundity == 0),
      fac_min = min(facilitator_density, na.rm = TRUE),
      fac_median = median(facilitator_density, na.rm = TRUE),
      fac_max = max(facilitator_density, na.rm = TRUE),
      fac_sd = sd(facilitator_density, na.rm = TRUE),
      fac_unique = dplyr::n_distinct(facilitator_density),
      fac_positive_n = sum(facilitator_density > 0, na.rm = TRUE),
      fac_positive_prop = mean(facilitator_density > 0, na.rm = TRUE),
      cons_unique = dplyr::n_distinct(conspecific_density),
      hetero_unique = dplyr::n_distinct(heterospecific_density),
      both_cons_hetero_positive_n = sum(
        conspecific_density > 0 & heterospecific_density > 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      weak_fac_support = fac_unique < 3 | fac_sd == 0 | fac_positive_n < 5,
      weak_density_pair_support = both_cons_hetero_positive_n < 5
    ) %>%
    dplyr::relocate(species)
}

support_table <- dplyr::bind_rows(
  interaction_support_from_data(df_glmm, "TRCY"),
  interaction_support_from_data(df_glmm, "TROR")
) %>%
  dplyr::arrange(species, Site, Cov)

log_decision("\n===== DATA SUPPORT FOR INTERACTION TERMS =====")
print_tbl(support_table, n = Inf)

############################################################
## 5. MODEL FORMULAS
############################################################

candidate_formulas <- list(
  ## Baseline environmental structures.
  additive = make_formula(
    "Site + Cov + cons_x + hetero_x + fac_x"
  ),
  site_cover = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x"
  ),
  
  ## Targeted a priori facilitation interactions.
  fac_by_site = make_formula(
    "Site + Cov + cons_x + hetero_x + fac_x + fac_x:Site"
  ),
  fac_by_cover = make_formula(
    "Site + Cov + cons_x + hetero_x + fac_x + fac_x:Cov"
  ),
  fac_by_site_cover = make_formula(
    "Site + Cov + cons_x + hetero_x + fac_x + fac_x:Site + fac_x:Cov"
  ),
  fac_by_density = make_formula(
    "Site + Cov + cons_x + hetero_x + fac_x + fac_x:cons_x + fac_x:hetero_x"
  ),
  site_cover_fac_by_density = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + fac_x:cons_x + fac_x:hetero_x"
  ),
  site_cover_fac_context = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + fac_x:Site + fac_x:Cov"
  ),
  broad_facilitation_context = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + fac_x:Site + fac_x:Cov + fac_x:cons_x + fac_x:hetero_x"
  ),
  
  ## Other ecological/diagnostic structures.
  site_density = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + Site:(cons_x + hetero_x + fac_x)"
  ),
  cover_density = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + Cov:(cons_x + hetero_x + fac_x)"
  ),
  quadratic_additive = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + I(cons_x^2) + I(hetero_x^2) + I(fac_x^2)"
  ),
  density_pairwise_diagnostic = make_formula(
    "Site * Cov + cons_x + hetero_x + fac_x + cons_x:hetero_x + cons_x:fac_x + hetero_x:fac_x"
  ),
  saturated_diagnostic = make_formula(
    "Site * Cov * (cons_x + hetero_x + fac_x) + cons_x:hetero_x + cons_x:fac_x + hetero_x:fac_x"
  )
)

apriori_lrt_pairs <- tibble::tribble(
  ~base_model, ~test_model,
  "additive", "fac_by_site",
  "additive", "fac_by_cover",
  "additive", "fac_by_site_cover",
  "additive", "fac_by_density",
  "site_cover", "site_cover_fac_by_density",
  "site_cover", "site_cover_fac_context",
  "site_cover", "broad_facilitation_context"
)

apriori_model_names <- unique(c(apriori_lrt_pairs$base_model, apriori_lrt_pairs$test_model))

############################################################
## 6. SPECIES-SPECIFIC PIPELINE
############################################################

fit_species_pipeline <- function(
    sp_name,
    dat,
    analysis_label = "primary",
    density_form_mode = c("auto", "raw", "log1p"),
    force_family = NULL,
    run_dharma = TRUE
) {
  
  density_form_mode <- match.arg(density_form_mode)
  
  log_decision("\n\n============================================================")
  log_decision("SPECIES: ", sp_name, " | ANALYSIS: ", analysis_label)
  log_decision("============================================================")
  
  ## Recalculate transformations within the exact species-specific dataset.
  ## This also applies to sensitivity datasets after outlier removal.
  d0 <- dat %>%
    dplyr::filter(focsp == sp_name) %>%
    droplevels() %>%
    dplyr::select(
      -any_of(c(
        "cons_raw_z", "hetero_raw_z", "fac_raw_z",
        "cons_log1p", "hetero_log1p", "fac_log1p",
        "cons_log1p_z", "hetero_log1p_z", "fac_log1p_z",
        "cons_x", "hetero_x", "fac_x"
      ))
    ) %>%
    add_density_transforms()
  
  density_forms_to_test <- if (density_form_mode == "auto") {
    c("raw", "log1p")
  } else {
    density_form_mode
  }
  
  ##########################################################
  ## 6A. Joint density-form and family screen
  ##########################################################
  
  family_density_models <- list()
  
  for (density_form in density_forms_to_test) {
    d_form <- assign_density_form(d0, density_form)
    f_screen <- fecundity ~ Site * Cov + cons_x + hetero_x + fac_x + (1 | SiteCovBlock)
    
    for (family_name in names(family_specs)) {
      if (!is.null(force_family) && family_name != force_family) next
      model_name <- paste(density_form, family_name, sep = "__")
      
      family_density_models[[model_name]] <- safe_fit(
        formula = f_screen,
        data = d_form,
        family = family_specs[[family_name]]$family,
        ziformula = family_specs[[family_name]]$ziformula,
        control = ctrl
      )
    }
  }
  
  family_density_tab <- aic_table(family_density_models, nrow(d0)) %>%
    tidyr::separate(
      model,
      into = c("density_form", "family_name"),
      sep = "__",
      remove = FALSE
    )
  
  log_decision("\n--- AIC SCREEN: density transformation and count family ---")
  print_tbl(family_density_tab, n = Inf)
  
  family_density_choice <- choose_best(family_density_models, family_density_tab, delta_rule = 2)
  
  selected_density_form <- family_density_choice$selected_table$density_form[[1]]
  selected_family_name <- family_density_choice$selected_table$family_name[[1]]
  
  log_decision(
    "SELECTED DENSITY FORM: ", selected_density_form,
    " | SELECTED FAMILY: ", selected_family_name
  )
  
  ##########################################################
  ## 6B. Fixed-effect and a priori interaction model selection
  ##########################################################
  
  d <- assign_density_form(d0, selected_density_form)
  selected_family <- family_specs[[selected_family_name]]
  
  candidate_models <- purrr::imap(
    candidate_formulas,
    ~ safe_fit(
      formula = .x,
      data = d,
      family = selected_family$family,
      ziformula = selected_family$ziformula,
      control = ctrl
    )
  )
  
  fixed_tab <- aic_table(candidate_models, nrow(d))
  
  log_decision("\n--- AIC SCREEN: fixed-effect and a priori interaction structure ---")
  print_tbl(fixed_tab, n = Inf)
  
  dropped_models <- fixed_tab %>%
    dplyr::filter(!usable) %>%
    dplyr::select(model, dropped_reason, nonestimable_terms, n_NA_fixed, n_NA_vcov)
  
  log_decision("\n--- DATA-DEFICIENT / DROPPED MODELS ---")
  if (nrow(dropped_models) == 0) {
    print_tbl(tibble::tibble(note = "No models were dropped for rank deficiency, convergence failure, or non-estimable coefficients."))
  } else {
    print_tbl(dropped_models, n = Inf)
  }
  
  fixed_choice <- choose_best(candidate_models, fixed_tab, delta_rule = 2)
  final_model <- fixed_choice$model
  final_model_name <- fixed_choice$name
  
  log_decision(
    "FINAL SELECTED MODEL: ", final_model_name,
    " | density_form = ", selected_density_form,
    " | family = ", selected_family_name
  )
  print(formula(final_model))
  print(summary(final_model))
  
  ##########################################################
  ## 6C. Manual LRTs for targeted a priori facilitation models
  ##########################################################
  
  apriori_model_table <- fixed_tab %>%
    dplyr::filter(model %in% apriori_model_names) %>%
    dplyr::mutate(
      species = sp_name,
      analysis_label = analysis_label,
      .before = model
    )
  
  apriori_lrt_table <- purrr::pmap_dfr(
    apriori_lrt_pairs,
    function(base_model, test_model) {
      manual_lrt_row(
        models = candidate_models,
        base_name = base_model,
        test_name = test_model,
        species = sp_name,
        n_obs = nrow(d)
      )
    }
  ) %>%
    dplyr::mutate(analysis_label = analysis_label, .after = species) %>%
    dplyr::arrange(delta_AIC_test_minus_base)
  
  log_decision("\n--- TARGETED A PRIORI FACILITATION MODELS: AIC/logLik ---")
  print_tbl(apriori_model_table, n = Inf)
  
  log_decision("\n--- MANUAL LRTs FOR TARGETED A PRIORI FACILITATION MODELS ---")
  print_tbl(apriori_lrt_table, n = Inf)
  
  broad_model <- candidate_models$broad_facilitation_context
  broad_row <- model_row(broad_model, "broad_facilitation_context", nrow(d))
  broad_type3 <- if (isTRUE(broad_row$usable)) {
    tryCatch(car::Anova(broad_model, type = 3), error = function(e) e)
  } else {
    paste("Broad targeted model was not usable; Type III tests not run.", broad_row$dropped_reason)
  }
  
  log_decision("\n--- TYPE III TESTS FOR BROAD TARGETED FACILITATION MODEL ---")
  print(broad_type3)
  
  ##########################################################
  ## 6D. Final outputs and diagnostics
  ##########################################################
  
  final_effects <- extract_fixed_effects(final_model) %>%
    dplyr::mutate(
      species = sp_name,
      analysis_label = analysis_label,
      density_form = selected_density_form,
      family_name = selected_family_name,
      final_model_name = final_model_name,
      n = nrow(d),
      AIC = AIC(final_model)
    )
  
  anova_type3 <- tryCatch(car::Anova(final_model, type = 3), error = function(e) e)
  
  log_decision("\n--- TYPE III WALD TESTS FOR FINAL MODEL ---")
  print(anova_type3)
  
  dharma_tests <- NULL
  dharma_flagged <- tibble::tibble()
  sim <- NULL
  
  if (run_dharma) {
    sim <- DHARMa::simulateResiduals(final_model, n = 1000, plot = FALSE)
    
    dharma_tests <- list(
      dispersion = DHARMa::testDispersion(sim),
      zero_inflation = DHARMa::testZeroInflation(sim),
      uniformity = DHARMa::testUniformity(sim),
      outliers = DHARMa::testOutliers(sim)
    )
    
    log_decision("\n--- DHARMa TESTS FOR FINAL MODEL ---")
    print(dharma_tests)
    
    scaled <- sim$scaledResiduals
    
    dharma_flagged <- d %>%
      dplyr::mutate(
        dharma_scaled_residual = scaled[seq_len(n())],
        dharma_lower_tail_025 = dharma_scaled_residual < 0.025,
        dharma_upper_tail_975 = dharma_scaled_residual > 0.975,
        dharma_tail_flag = dharma_lower_tail_025 | dharma_upper_tail_975,
        dharma_exact_flag = dharma_scaled_residual %in% c(0, 1)
      ) %>%
      dplyr::filter(dharma_tail_flag | dharma_exact_flag) %>%
      dplyr::arrange(desc(fecundity)) %>%
      dplyr::select(
        obs_id, focsp, Site, Cov, BLK,
        fecundity, conspecific_density, heterospecific_density, facilitator_density,
        TRCY, TROR, NO_GRASS,
        dharma_scaled_residual, dharma_tail_flag, dharma_exact_flag
      )
    
    log_decision("\n--- DHARMa-FLAGGED OBSERVATIONS FOR SENSITIVITY ONLY ---")
    if (nrow(dharma_flagged) == 0) {
      print_tbl(tibble::tibble(note = "No observations flagged using residual tails <0.025 or >0.975, or exact 0/1 scaled residuals."))
    } else {
      print_tbl(dharma_flagged, n = Inf)
    }
  }
  
  list(
    species = sp_name,
    analysis_label = analysis_label,
    data = d,
    density_family_AIC = family_density_tab,
    fixed_AIC = fixed_tab,
    dropped_models = dropped_models,
    selected_density_form = selected_density_form,
    selected_family_name = selected_family_name,
    final_model_name = final_model_name,
    final_model = final_model,
    final_effects = final_effects,
    type3 = anova_type3,
    candidate_models = candidate_models,
    apriori_model_table = apriori_model_table,
    apriori_lrt_table = apriori_lrt_table,
    broad_facilitation_type3 = broad_type3,
    dharma_tests = dharma_tests,
    dharma_flagged = dharma_flagged,
    dharma_sim = sim
  )
}

############################################################
## 7. PRIMARY ANALYSES: FULL DATA, AUTO RAW VS LOG1P CHOICE
############################################################

fit_TRCY <- fit_species_pipeline(
  sp_name = "TRCY",
  dat = df_glmm,
  analysis_label = "primary_full_data_auto_density",
  density_form_mode = "auto",
  run_dharma = TRUE
)

fit_TROR <- fit_species_pipeline(
  sp_name = "TROR",
  dat = df_glmm,
  analysis_label = "primary_full_data_auto_density",
  density_form_mode = "auto",
  run_dharma = TRUE
)

############################################################
## 8. TROR OUTLIER SENSITIVITY ANALYSES
############################################################

tror_dharma_ids <- fit_TROR$dharma_flagged %>%
  dplyr::pull(obs_id) %>%
  unique()

tror_high_fec_ids <- high_fec_report %>%
  dplyr::filter(focsp == "TROR") %>%
  dplyr::pull(obs_id) %>%
  unique()

tror_sensitivity_ids <- if (length(tror_dharma_ids) > 0) {
  log_decision("OUTLIER SENSITIVITY: using DHARMa-flagged TROR obs_id values.")
  tror_dharma_ids
} else {
  log_decision("OUTLIER SENSITIVITY: no DHARMa flags; using high-fecundity TROR obs_id values.")
  tror_high_fec_ids
}

tror_outlier_record <- df_glmm %>%
  dplyr::filter(obs_id %in% tror_sensitivity_ids) %>%
  dplyr::arrange(desc(fecundity)) %>%
  dplyr::select(
    obs_id, focsp, Site, Cov, BLK,
    fecundity, conspecific_density, heterospecific_density, facilitator_density,
    TRCY, TROR, NO_GRASS
  )

log_decision("\n===== TROR OBSERVATIONS REMOVED IN OUTLIER SENSITIVITY ANALYSES =====")
print_tbl(tror_outlier_record, n = Inf)

df_glmm_no_TROR_outliers <- df_glmm %>%
  dplyr::filter(!(focsp == "TROR" & obs_id %in% tror_sensitivity_ids))

fit_TROR_no_outliers_auto <- NULL
fit_TROR_no_outliers_raw <- NULL

if (length(tror_sensitivity_ids) > 0) {
  
  fit_TROR_no_outliers_auto <- fit_species_pipeline(
    sp_name = "TROR",
    dat = df_glmm_no_TROR_outliers,
    analysis_label = "TROR_no_outliers_auto_density",
    density_form_mode = "auto",
    run_dharma = TRUE
  )
  
  fit_TROR_no_outliers_raw <- fit_species_pipeline(
    sp_name = "TROR",
    dat = df_glmm_no_TROR_outliers,
    analysis_label = "TROR_no_outliers_raw_density_forced",
    density_form_mode = "raw",
    run_dharma = TRUE
  )
  
} else {
  log_decision("No TROR sensitivity models run because no outlier observations were identified.")
}

############################################################
## 9. COMPACT SAVED OUTPUTS
############################################################

optional_bind <- function(...) {
  dplyr::bind_rows(purrr::compact(list(...)))
}

aic_tables <- optional_bind(
  fit_TRCY$density_family_AIC %>%
    dplyr::mutate(species = "TRCY", analysis_label = fit_TRCY$analysis_label, selection_stage = "density_family"),
  fit_TROR$density_family_AIC %>%
    dplyr::mutate(species = "TROR", analysis_label = fit_TROR$analysis_label, selection_stage = "density_family"),
  fit_TRCY$fixed_AIC %>%
    dplyr::mutate(species = "TRCY", analysis_label = fit_TRCY$analysis_label, selection_stage = "fixed_structure"),
  fit_TROR$fixed_AIC %>%
    dplyr::mutate(species = "TROR", analysis_label = fit_TROR$analysis_label, selection_stage = "fixed_structure"),
  if (!is.null(fit_TROR_no_outliers_auto)) {
    fit_TROR_no_outliers_auto$density_family_AIC %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_auto$analysis_label, selection_stage = "density_family")
  },
  if (!is.null(fit_TROR_no_outliers_auto)) {
    fit_TROR_no_outliers_auto$fixed_AIC %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_auto$analysis_label, selection_stage = "fixed_structure")
  },
  if (!is.null(fit_TROR_no_outliers_raw)) {
    fit_TROR_no_outliers_raw$density_family_AIC %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_raw$analysis_label, selection_stage = "density_family")
  },
  if (!is.null(fit_TROR_no_outliers_raw)) {
    fit_TROR_no_outliers_raw$fixed_AIC %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_raw$analysis_label, selection_stage = "fixed_structure")
  }
)

final_model_summary <- optional_bind(
  fit_TRCY$final_effects,
  fit_TROR$final_effects,
  if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$final_effects,
  if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$final_effects
)

dropped_model_summary <- optional_bind(
  fit_TRCY$dropped_models %>%
    dplyr::mutate(species = "TRCY", analysis_label = fit_TRCY$analysis_label),
  fit_TROR$dropped_models %>%
    dplyr::mutate(species = "TROR", analysis_label = fit_TROR$analysis_label),
  if (!is.null(fit_TROR_no_outliers_auto)) {
    fit_TROR_no_outliers_auto$dropped_models %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_auto$analysis_label)
  },
  if (!is.null(fit_TROR_no_outliers_raw)) {
    fit_TROR_no_outliers_raw$dropped_models %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_raw$analysis_label)
  }
)

dharma_flag_summary <- optional_bind(
  fit_TRCY$dharma_flagged %>%
    dplyr::mutate(species = "TRCY", analysis_label = fit_TRCY$analysis_label),
  fit_TROR$dharma_flagged %>%
    dplyr::mutate(species = "TROR", analysis_label = fit_TROR$analysis_label),
  if (!is.null(fit_TROR_no_outliers_auto)) {
    fit_TROR_no_outliers_auto$dharma_flagged %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_auto$analysis_label)
  },
  if (!is.null(fit_TROR_no_outliers_raw)) {
    fit_TROR_no_outliers_raw$dharma_flagged %>%
      dplyr::mutate(species = "TROR", analysis_label = fit_TROR_no_outliers_raw$analysis_label)
  }
)

apriori_model_summary <- optional_bind(
  fit_TRCY$apriori_model_table,
  fit_TROR$apriori_model_table,
  if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$apriori_model_table,
  if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$apriori_model_table
)

apriori_lrt_summary <- optional_bind(
  fit_TRCY$apriori_lrt_table,
  fit_TROR$apriori_lrt_table,
  if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$apriori_lrt_table,
  if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$apriori_lrt_table
)

safe_write_csv(response_summary, file.path(outdir, "01_response_summary.csv"))
safe_write_csv(density_transform_summary, file.path(outdir, "02_density_transformation_summary.csv"))
safe_write_csv(data_support_by_cell, file.path(outdir, "04_site_cover_data_support.csv"))
safe_write_csv(aic_tables, file.path(outdir, "07_AIC_selection_tables.csv"))
safe_write_csv(dropped_model_summary, file.path(outdir, "08_dropped_or_rank_deficient_models.csv"))
safe_write_csv(final_model_summary, file.path(outdir, "09_final_model_effects.csv"))
safe_write_csv(support_table, file.path(outdir, "11_interaction_data_support_by_species_site_cover.csv"))
safe_write_csv(apriori_model_summary, file.path(outdir, "12_apriori_facilitation_AIC_logLik_table.csv"))
safe_write_csv(apriori_lrt_summary, file.path(outdir, "13_apriori_facilitation_manual_LRT_table.csv"))

if (isTRUE(WRITE_EXTENDED_TABLES)) {
  safe_write_csv(z_check, file.path(outdir, "03_z_transformation_check.csv"))
  safe_write_csv(high_fec_report, file.path(outdir, "05_high_fecundity_observations_for_sensitivity.csv"))
  safe_write_csv(tror_outlier_record, file.path(outdir, "06_TROR_removed_observations_sensitivity.csv"))
  safe_write_csv(dharma_flag_summary, file.path(outdir, "10_DHARMa_flagged_observations.csv"))
}

writeLines(decision_log, file.path(outdir, "00_decision_log.txt"))

############################################################
## 10. FINAL CONSOLE SUMMARY
############################################################

cat("\n\n============================================================\n")
cat("FINAL MODEL SELECTION SUMMARY\n")
cat("============================================================\n")

final_selection_table <- tibble::tibble(
  species = c("TRCY", "TROR",
              if (!is.null(fit_TROR_no_outliers_auto)) "TROR" else NULL,
              if (!is.null(fit_TROR_no_outliers_raw)) "TROR" else NULL),
  analysis_label = c(
    fit_TRCY$analysis_label,
    fit_TROR$analysis_label,
    if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$analysis_label else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$analysis_label else NULL
  ),
  n = c(
    nrow(fit_TRCY$data),
    nrow(fit_TROR$data),
    if (!is.null(fit_TROR_no_outliers_auto)) nrow(fit_TROR_no_outliers_auto$data) else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) nrow(fit_TROR_no_outliers_raw$data) else NULL
  ),
  density_form = c(
    fit_TRCY$selected_density_form,
    fit_TROR$selected_density_form,
    if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$selected_density_form else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$selected_density_form else NULL
  ),
  family = c(
    fit_TRCY$selected_family_name,
    fit_TROR$selected_family_name,
    if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$selected_family_name else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$selected_family_name else NULL
  ),
  final_model = c(
    fit_TRCY$final_model_name,
    fit_TROR$final_model_name,
    if (!is.null(fit_TROR_no_outliers_auto)) fit_TROR_no_outliers_auto$final_model_name else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) fit_TROR_no_outliers_raw$final_model_name else NULL
  ),
  AIC = c(
    AIC(fit_TRCY$final_model),
    AIC(fit_TROR$final_model),
    if (!is.null(fit_TROR_no_outliers_auto)) AIC(fit_TROR_no_outliers_auto$final_model) else NULL,
    if (!is.null(fit_TROR_no_outliers_raw)) AIC(fit_TROR_no_outliers_raw$final_model) else NULL
  )
)

print_tbl(final_selection_table, n = Inf)
safe_write_csv(final_selection_table, file.path(outdir, "15_final_model_selection_summary.csv"))

cat("\n\n============================================================\n")
cat("DROPPED / DATA-DEFICIENT MODELS\n")
cat("============================================================\n")
if (nrow(dropped_model_summary) == 0) {
  print_tbl(tibble::tibble(note = "No models dropped."))
} else {
  print_tbl(dropped_model_summary, n = Inf)
}

cat("\n\n============================================================\n")
cat("A PRIORI FACILITATION MANUAL LRT SUMMARY\n")
cat("============================================================\n")
print_tbl(apriori_lrt_summary, n = Inf)

cat("\n\nEssential output files saved to: ", outdir, "\n", sep = "")
cat("Figure folder for downstream plots: ", figdir, "\n", sep = "")

############################################################
## SAVE FINAL GLMM MODEL OBJECTS FOR PLOTTING
############################################################

model_dir <- file.path(outdir, "model_objects")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

safe_save_rds(
  fit_TRCY$final_model,
  file.path(model_dir, "TRCY_primary_final_model.rds")
)

safe_save_rds(
  fit_TROR$final_model,
  file.path(model_dir, "TROR_primary_final_model.rds")
)

safe_save_rds(
  list(TRCY = fit_TRCY, TROR = fit_TROR),
  file.path(model_dir, "primary_fit_objects_full.rds")
)

## Optional: save TROR outlier-sensitivity models if they exist.
if (exists("fit_TROR_no_outliers_auto") && !is.null(fit_TROR_no_outliers_auto)) {
  safe_save_rds(
    fit_TROR_no_outliers_auto$final_model,
    file.path(model_dir, "TROR_no_outliers_auto_final_model.rds")
  )
}

if (exists("fit_TROR_no_outliers_raw") && !is.null(fit_TROR_no_outliers_raw)) {
  safe_save_rds(
    fit_TROR_no_outliers_raw$final_model,
    file.path(model_dir, "TROR_no_outliers_raw_final_model.rds")
  )
}

############################################################
## PRINT FINAL MODEL RESULTS
############################################################
summary(fit_TRCY$final_model)
summary(fit_TROR$final_model)

############################################################
## POST HOC CONTRASTS USING EMMEANS
##
## Run after:
##   fit_TRCY <- fit_species_pipeline(...)
##   fit_TROR <- fit_species_pipeline(...)
##
## Notes:
##   - Models use a log link.
##   - Site and cover contrasts are summarised as response ratios.
##   - Percent change = (ratio - 1) * 100.
##   - Site contrasts are compared against Ben.
##   - Cover contrast is Sun vs Shade.
##   - Density/facilitation trends are on the log expected fecundity scale
##     and are converted to percent change per one-SD increase.
############################################################

pkgs_posthoc <- c("emmeans", "multcomp", "multcompView")
missing_posthoc <- pkgs_posthoc[!pkgs_posthoc %in% rownames(installed.packages())]
if (length(missing_posthoc) > 0) install.packages(missing_posthoc)
invisible(lapply(pkgs_posthoc, library, character.only = TRUE))

emmeans::emm_options(rg.limit = 50000)

posthoc_dir <- file.path(outdir, "posthoc_emmeans")
dir.create(posthoc_dir, showWarnings = FALSE, recursive = TRUE)

save_posthoc <- function(x, filename) {
  path <- file.path(posthoc_dir, filename)
  readr::write_csv(as.data.frame(x), path)
  message("Saved: ", normalizePath(path, winslash = "/", mustWork = TRUE))
  invisible(path)
}

############################################################
## Main post hoc function
############################################################

run_emmeans_posthoc <- function(fit_obj, species_code) {
  
  model <- fit_obj$final_model
  
  cat("\n\n============================================================\n")
  cat("EMMEANS POST HOC CONTRASTS:", species_code, "\n")
  cat("============================================================\n")
  
  ##########################################################
  ## 1. Site emmeans and contrasts
  ##########################################################
  
  emm_site_link <- emmeans::emmeans(
    model,
    ~ Site,
    type = "link"
  )
  
  emm_site_resp <- emmeans::emmeans(
    model,
    ~ Site,
    type = "response"
  )
  
  site_all_pairs_link <- pairs(
    emm_site_link,
    adjust = "tukey"
  )
  
  site_all_pairs_resp <- summary(
    site_all_pairs_link,
    type = "response",
    infer = c(TRUE, TRUE)
  )
  
  site_all_pairs_percent <- as_percent_contrast(site_all_pairs_resp)
  
  ben_ref <- which(levels(fit_obj$data$Site) == "Ben")
  
  site_vs_ben_link <- contrast(
    emm_site_link,
    method = "trt.vs.ctrl",
    ref = ben_ref,
    adjust = "dunnettx"
  )
  
  site_vs_ben_resp <- summary(
    site_vs_ben_link,
    type = "response",
    infer = c(TRUE, TRUE)
  )
  
  site_vs_ben_percent <- as_percent_contrast(site_vs_ben_resp)
  
  site_letters <- make_letters_from_pairs(
    emm_resp = emm_site_resp,
    pair_resp = site_all_pairs_resp,
    group_col = "Site",
    alpha = 0.05
  )
  
  cat("\n--- Site estimated marginal means, response scale ---\n")
  print(emm_site_resp)
  
  cat("\n--- Site contrasts versus Ben, response ratios ---\n")
  print(site_vs_ben_resp)
  
  cat("\n--- Site contrasts versus Ben, percent change ---\n")
  print_posthoc(site_vs_ben_percent, n = Inf)
  
  cat("\n--- Site compact letter display ---\n")
  print_posthoc(site_letters, n = Inf)
  
  ##########################################################
  ## 2. Cover emmeans and contrast
  ##########################################################
  
  emm_cov_link <- emmeans::emmeans(
    model,
    ~ Cov,
    type = "link"
  )
  
  emm_cov_resp <- emmeans::emmeans(
    model,
    ~ Cov,
    type = "response"
  )
  
  cov_pair_link <- pairs(
    emm_cov_link,
    adjust = "none"
  )
  
  cov_pair_resp <- summary(
    cov_pair_link,
    type = "response",
    infer = c(TRUE, TRUE)
  )
  
  cov_pair_percent <- as_percent_contrast(cov_pair_resp)
  
  cat("\n--- Cover estimated marginal means, response scale ---\n")
  print(emm_cov_resp)
  
  cat("\n--- Cover contrast, response ratio ---\n")
  print(cov_pair_resp)
  
  cat("\n--- Cover contrast, percent change ---\n")
  print_posthoc(cov_pair_percent, n = Inf)
  
  ##########################################################
  ## 3. Density and facilitator trends
  ##########################################################
  
  trend_cons <- emmeans::emtrends(
    model,
    ~ 1,
    var = "cons_x"
  )
  
  trend_hetero <- emmeans::emtrends(
    model,
    ~ 1,
    var = "hetero_x"
  )
  
  trend_fac <- emmeans::emtrends(
    model,
    ~ 1,
    var = "fac_x"
  )
  
  trend_cons_sum <- summary(trend_cons, infer = c(TRUE, TRUE))
  trend_hetero_sum <- summary(trend_hetero, infer = c(TRUE, TRUE))
  trend_fac_sum <- summary(trend_fac, infer = c(TRUE, TRUE))
  
  trend_cons_percent <- as_percent_trend(trend_cons_sum, "cons_x.trend")
  trend_hetero_percent <- as_percent_trend(trend_hetero_sum, "hetero_x.trend")
  trend_fac_percent <- as_percent_trend(trend_fac_sum, "fac_x.trend")
  
  cat("\n--- Conspecific density trend ---\n")
  print(trend_cons_sum)
  
  cat("\n--- Heterospecific density trend ---\n")
  print(trend_hetero_sum)
  
  cat("\n--- Facilitator density trend ---\n")
  print(trend_fac_sum)
  
  ##########################################################
  ## 4. TROR-specific facilitator-by-density trends
  ##########################################################
  
  fac_trends_by_density <- NULL
  fac_trend_pairwise <- NULL
  
  if (species_code == "TROR") {
    
    fac_trends_by_cons <- emmeans::emtrends(
      model,
      ~ cons_x,
      var = "fac_x",
      at = list(
        cons_x = c(-1, 0, 1),
        hetero_x = 0
      )
    )
    
    fac_trends_by_hetero <- emmeans::emtrends(
      model,
      ~ hetero_x,
      var = "fac_x",
      at = list(
        hetero_x = c(-1, 0, 1),
        cons_x = 0
      )
    )
    
    fac_cons_sum <- summary(fac_trends_by_cons, infer = c(TRUE, TRUE))
    fac_hetero_sum <- summary(fac_trends_by_hetero, infer = c(TRUE, TRUE))
    
    fac_trends_by_density <- dplyr::bind_rows(
      as_percent_trend(fac_cons_sum, "fac_x.trend") %>%
        dplyr::mutate(context_axis = "conspecific_density"),
      as_percent_trend(fac_hetero_sum, "fac_x.trend") %>%
        dplyr::mutate(context_axis = "heterospecific_density")
    )
    
    fac_trend_pairwise <- dplyr::bind_rows(
      as.data.frame(summary(
        pairs(fac_trends_by_cons, adjust = "tukey"),
        infer = c(TRUE, TRUE)
      )) %>%
        dplyr::mutate(context_axis = "conspecific_density"),
      as.data.frame(summary(
        pairs(fac_trends_by_hetero, adjust = "tukey"),
        infer = c(TRUE, TRUE)
      )) %>%
        dplyr::mutate(context_axis = "heterospecific_density")
    )
    
    cat("\n--- TROR facilitator slope across neighbour-density context ---\n")
    print_posthoc(fac_trends_by_density, n = Inf)
    
    cat("\n--- TROR pairwise contrasts among facilitator slopes ---\n")
    print_posthoc(fac_trend_pairwise, n = Inf)  }
  
  ##########################################################
  ## 5. Save outputs
  ##########################################################
  
  if (isTRUE(SAVE_DETAILED_POSTHOC)) {
    save_posthoc(as.data.frame(emm_site_resp), paste0(species_code, "_emmeans_site_response.csv"))
    save_posthoc(site_all_pairs_percent, paste0(species_code, "_site_all_pairwise_response_ratio_percent.csv"))
    save_posthoc(site_vs_ben_percent, paste0(species_code, "_site_vs_Ben_response_ratio_percent.csv"))
    save_posthoc(site_letters, paste0(species_code, "_site_letters_tukey.csv"))
    
    save_posthoc(as.data.frame(emm_cov_resp), paste0(species_code, "_emmeans_cover_response.csv"))
    save_posthoc(cov_pair_percent, paste0(species_code, "_cover_contrast_response_ratio_percent.csv"))
    
    save_posthoc(trend_cons_percent, paste0(species_code, "_trend_conspecific_density_percent.csv"))
    save_posthoc(trend_hetero_percent, paste0(species_code, "_trend_heterospecific_density_percent.csv"))
    save_posthoc(trend_fac_percent, paste0(species_code, "_trend_facilitator_density_percent.csv"))
    
    if (!is.null(fac_trends_by_density)) {
      save_posthoc(
        fac_trends_by_density,
        paste0(species_code, "_facilitator_trends_by_density_context.csv")
      )
    }
    
    if (!is.null(fac_trend_pairwise)) {
      save_posthoc(
        fac_trend_pairwise,
        paste0(species_code, "_facilitator_trend_pairwise_by_density_context.csv")
      )
    }
    
  }
  
  list(
    species = species_code,
    emm_site_response = emm_site_resp,
    site_all_pairs_percent = site_all_pairs_percent,
    site_vs_ben_percent = site_vs_ben_percent,
    site_letters = site_letters,
    emm_cover_response = emm_cov_resp,
    cover_pair_percent = cov_pair_percent,
    trend_cons_percent = trend_cons_percent,
    trend_hetero_percent = trend_hetero_percent,
    trend_fac_percent = trend_fac_percent,
    fac_trends_by_density = fac_trends_by_density,
    fac_trend_pairwise = fac_trend_pairwise
  )
}

############################################################
## Run post hoc contrasts
############################################################

posthoc_TRCY <- run_emmeans_posthoc(fit_TRCY, "TRCY")
posthoc_TROR <- run_emmeans_posthoc(fit_TROR, "TROR")

############################################################
## Combined compact summaries
############################################################

site_vs_ben_posthoc <- dplyr::bind_rows(
  posthoc_TRCY$site_vs_ben_percent %>% dplyr::mutate(species = "TRCY"),
  posthoc_TROR$site_vs_ben_percent %>% dplyr::mutate(species = "TROR")
)

cover_posthoc <- dplyr::bind_rows(
  posthoc_TRCY$cover_pair_percent %>% dplyr::mutate(species = "TRCY"),
  posthoc_TROR$cover_pair_percent %>% dplyr::mutate(species = "TROR")
)

density_trend_posthoc <- dplyr::bind_rows(
  posthoc_TRCY$trend_cons_percent %>%
    dplyr::mutate(species = "TRCY", density = "conspecific"),
  posthoc_TRCY$trend_hetero_percent %>%
    dplyr::mutate(species = "TRCY", density = "heterospecific"),
  posthoc_TRCY$trend_fac_percent %>%
    dplyr::mutate(species = "TRCY", density = "facilitator"),
  posthoc_TROR$trend_cons_percent %>%
    dplyr::mutate(species = "TROR", density = "conspecific"),
  posthoc_TROR$trend_hetero_percent %>%
    dplyr::mutate(species = "TROR", density = "heterospecific"),
  posthoc_TROR$trend_fac_percent %>%
    dplyr::mutate(species = "TROR", density = "facilitator")
)

safe_write_csv(
  site_vs_ben_posthoc,
  file.path(posthoc_dir, "combined_site_vs_Ben_response_ratio_percent.csv")
)

safe_write_csv(
  cover_posthoc,
  file.path(posthoc_dir, "combined_cover_response_ratio_percent.csv")
)

safe_write_csv(
  density_trend_posthoc,
  file.path(posthoc_dir, "combined_density_trends_percent.csv")
)

if (!is.null(posthoc_TROR$fac_trends_by_density)) {
  safe_write_csv(
    posthoc_TROR$fac_trends_by_density,
    file.path(posthoc_dir, "combined_TROR_facilitator_trends_by_density_context.csv")
  )
}

if (!is.null(posthoc_TROR$fac_trend_pairwise)) {
  safe_write_csv(
    posthoc_TROR$fac_trend_pairwise,
    file.path(posthoc_dir, "combined_TROR_facilitator_trend_pairwise_by_density_context.csv")
  )
}

cat("\nPost hoc emmeans outputs saved to: ", posthoc_dir, "\n", sep = "")