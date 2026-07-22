############################################################
## 02_plot_outcomes_BH_Ricker.R
##
## Purpose:
##   Recalculate coexistence outcomes and figures for both final
##   Ricker and Beverton-Holt (BH) Bayesian models
##   - runs Ricker and BH in one script;
##   - writes figures to coexistence_models/02_outcomes/figures;
##   - compares BH versus Ricker posterior outcome probabilities;
##   - uses the previously working contour plotting structure with model-specific BH/Ricker calculations;
##   - prints and saves a summary of how different the estimated
##     outcomes are between model families.
############################################################

rm(list = ls())

############################
## 0. Packages
############################

required_packages <- c(
  "tidyverse",
  "brms",
  "posterior",
  "ggrepel",
  "patchwork",
  "ggnewscale",
  "here",
  "readr",
  "MASS",
  "grid"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

options(mc.cores = max(1, parallel::detectCores() - 1))


############################
## 1. User options
############################

## Leave NULL unless the automatic project finder fails.
## This should be the folder containing 01_alpha_structure_selection
## and 02_outcomes, or the folder containing coexistence_models.
project_dir_override <- NULL

## Run both final model families.
families_to_run <- c("Ricker", "BH")

## Manual fit overrides. Leave as NA unless automatic fit search fails.
## Example:
## fit_file_overrides["BH"] <- "coexistence_models/01_alpha_structure_selection/fits/fit_BH_m1_loggaussian.rds"
fit_file_overrides <- c(
  Ricker = NA_character_,
  BH     = NA_character_
)

## Fit-search terms used only when a manual override is not supplied.
## The search is intentionally broader than the previous script because
## your observed selected Ricker fit was named fit_Ricker_m1_loggaussian.rds,
## not fit_Ricker_m1_species.rds.
fit_search_terms <- c("m1_species", "m1_loggaussian", "m1")

## Main default: no germination correction.
## If TRUE, lambda is multiplied by germination draws from upstream priors.
use_germination <- FALSE

## high facilitator = mean of positive fac_sc values >= this quantile.
high_fac_quantile <- 0.90

## If TRUE, keeps only draws where both species have positive single-species growth
## for niche/fitness-ratio summaries.
restrict_kappa_to_positive_growth <- FALSE

## Figure options
save_png_dpi_main <- 300
save_png_dpi_contours <- 600


############################
## 2. Output/path helpers
############################

sanitize_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

short_error_code <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    identical(x, "gaussian_logfit") ~ "lg",
    identical(x, "loggaussian") ~ "lg",
    grepl("gaussian", x, ignore.case = TRUE) & grepl("log", x, ignore.case = TRUE) ~ "lg",
    grepl("negbin", x, ignore.case = TRUE) ~ "nb",
    grepl("poisson", x, ignore.case = TRUE) ~ "pois",
    TRUE ~ tolower(substr(sanitize_token(x), 1, 8))
  )
}

short_family_code <- function(family_tag) {
  family_tag <- as.character(family_tag)
  dplyr::case_when(
    toupper(family_tag) == "RICKER" ~ "ri",
    toupper(family_tag) == "BH" ~ "bh",
    TRUE ~ tolower(substr(sanitize_token(family_tag), 1, 8))
  )
}

make_file_prefix <- function(family_tag, error_tag, use_germination) {
  paste(
    short_family_code(family_tag),
    short_error_code(error_tag),
    ifelse(use_germination, "g", "ng"),
    sep = "_"
  )
}

table_path <- function(prefix, tag, ext = "csv") {
  file.path(tables_dir, paste0(prefix, "_", tag, ".", ext))
}

draw_path <- function(prefix, tag, ext = "rds") {
  file.path(draws_dir, paste0(prefix, "_", tag, ".", ext))
}

fig_path <- function(prefix, tag, ext = "png") {
  file.path(figures_dir, paste0(prefix, "_", tag, ".", ext))
}

summary_path <- function(prefix, tag, ext = "csv") {
  file.path(summaries_dir, paste0(prefix, "_", tag, ".", ext))
}

safe_path_report <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  path_report <- safe_path_report(path)
  if (nchar(path_report) > 240) {
    warning(
      "Output path is still long and may fail on Windows/OneDrive:\n",
      path_report,
      "\nShorten the project folder if this fails."
    )
  }

  readr::write_csv(x, path)

  if (!file.exists(path)) {
    stop("File was not created: ", path_report)
  }

  cat("Saved table: ", path_report, "\n", sep = "")
  invisible(path)
}

safe_save_rds <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  path_report <- safe_path_report(path)
  if (nchar(path_report) > 240) {
    warning(
      "Output path is still long and may fail on Windows/OneDrive:\n",
      path_report,
      "\nShorten the project folder if this fails."
    )
  }

  saveRDS(x, path)

  if (!file.exists(path)) {
    stop("RDS file was not created: ", path_report)
  }

  cat("Saved RDS: ", path_report, "\n", sep = "")
  invisible(path)
}

safe_save_plot <- function(plot, path, width, height, dpi) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  path_report <- safe_path_report(path)
  if (nchar(path_report) > 240) {
    warning(
      "Figure path is still long and may fail on Windows/OneDrive:\n",
      path_report,
      "\nShorten the project folder if this fails."
    )
  }

  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  if (!file.exists(path)) {
    stop("Figure file was not created: ", path_report)
  }

  cat("Saved figure: ", path_report, "\n", sep = "")
  invisible(path)
}

check_write_dir <- function(folder) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  test_file <- file.path(folder, paste0(".write_test_", Sys.getpid(), ".tmp"))

  ok <- tryCatch(
    {
      writeLines("ok", test_file)
      file.exists(test_file)
    },
    error = function(e) {
      message("Write test failed for: ", safe_path_report(folder))
      message("Error: ", conditionMessage(e))
      FALSE
    }
  )

  if (file.exists(test_file)) {
    unlink(test_file)
  }

  if (!isTRUE(ok)) {
    stop(
      "Output folder is not writable: ", safe_path_report(folder), "\n",
      "This is not a model issue. It is a path/permission/OneDrive issue."
    )
  }

  invisible(TRUE)
}

resolve_existing_file <- function(path,
                                  fallback_dirs = character(),
                                  label = "file",
                                  prefer_current = TRUE,
                                  search_pattern = NULL) {

  fallback_dirs <- unique(fallback_dirs[dir.exists(fallback_dirs)])
  original_path <- path

  ## Prefer files in the current project folders. alpha_info can contain
  ## absolute paths from earlier runs; those should not win over current files.
  if (!is.null(path) && length(path) == 1 && !is.na(path) && nzchar(path)) {

    path_base <- basename(path)

    if (prefer_current && nzchar(path_base)) {
      current_candidates <- file.path(fallback_dirs, path_base)
      current_candidates <- current_candidates[file.exists(current_candidates)]

      if (length(current_candidates) > 0) {
        newest <- current_candidates[which.max(file.info(current_candidates)$mtime)]
        return(normalizePath(newest, winslash = "/", mustWork = TRUE))
      }
    }

    if (!prefer_current && file.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }

    if (!prefer_current && nzchar(path_base)) {
      for (d in fallback_dirs) {
        candidate <- file.path(d, path_base)
        if (file.exists(candidate)) {
          return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
        }
      }
    }

    if (file.exists(path)) {
      warning(
        "Stored path exists but was not used for ", label, " because prefer_current = TRUE.\n",
        "This prevents accidental reuse of stale absolute paths from earlier model runs.\n",
        "Stored path was:\n",
        normalizePath(path, winslash = "/", mustWork = TRUE)
      )
    }
  }

  if (!is.null(search_pattern) && length(search_pattern) == 1 && nzchar(search_pattern)) {
    searched <- unlist(
      lapply(
        fallback_dirs,
        function(d) list.files(d, pattern = search_pattern, full.names = TRUE, recursive = FALSE)
      ),
      use.names = FALSE
    )

    searched <- searched[file.exists(searched)]

    if (length(searched) > 0) {
      searched <- searched[order(file.info(searched)$mtime, decreasing = TRUE)]
      warning(
        "Could not resolve ", label, " by stored path/basename. ",
        "Using most recently modified current-project match:\n",
        searched[1]
      )
      return(normalizePath(searched[1], winslash = "/", mustWork = TRUE))
    }
  }

  stop(
    "Could not resolve ", label, ".\n",
    "Stored path was: ", original_path, "\n",
    "Current-project folders checked:\n",
    paste(fallback_dirs, collapse = "\n"),
    if (!is.null(search_pattern)) paste0("\nSearch pattern used: ", search_pattern) else ""
  )
}

find_project_dir <- function(project_dir_override = NULL) {

  if (!is.null(project_dir_override) &&
      length(project_dir_override) == 1 &&
      !is.na(project_dir_override) &&
      nzchar(project_dir_override)) {

    if (dir.exists(file.path(project_dir_override, "01_alpha_structure_selection"))) {
      return(normalizePath(project_dir_override, winslash = "/", mustWork = TRUE))
    }

    if (dir.exists(file.path(project_dir_override, "coexistence_models", "01_alpha_structure_selection"))) {
      return(normalizePath(file.path(project_dir_override, "coexistence_models"), winslash = "/", mustWork = TRUE))
    }

    stop(
      "project_dir_override was supplied but 01_alpha_structure_selection was not found under it:\n",
      project_dir_override
    )
  }

  if (dir.exists(here::here("coexistence_models", "01_alpha_structure_selection"))) {
    return(normalizePath(here::here("coexistence_models"), winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(here::here("01_alpha_structure_selection"))) {
    return(normalizePath(here::here(), winslash = "/", mustWork = TRUE))
  }

  stop(
    "Could not find 01_alpha_structure_selection. Check your R project root.\n",
    "Current here::here() is: ", here::here(), "\n",
    "Set project_dir_override manually if needed."
  )
}

find_fit_for_family <- function(family_tag, structure_dir, fit_search_terms = c("m1_species", "m1_loggaussian", "m1")) {

  family_tag <- dplyr::case_when(
    toupper(family_tag) == "BH" ~ "BH",
    toupper(family_tag) == "RICKER" ~ "Ricker",
    TRUE ~ as.character(family_tag)
  )

  fits_dir <- file.path(structure_dir, "fits")

  if (!dir.exists(fits_dir)) {
    stop("Cannot search for fit because fits_dir does not exist: ", fits_dir)
  }

  all_fit_files <- list.files(fits_dir, pattern = "\\.rds$", full.names = TRUE, recursive = FALSE)

  if (length(all_fit_files) == 0) {
    stop("No .rds fit files found in: ", fits_dir)
  }

  bn <- basename(all_fit_files)

  not_aux <- !grepl("loo|psis|kfold|compare|summary|draw|posterior", bn, ignore.case = TRUE)

  family_match <- if (family_tag == "BH") {
    grepl("(^|[_-])BH([_-]|\\.|$)", bn, ignore.case = TRUE)
  } else {
    grepl("Ricker", bn, ignore.case = TRUE)
  }

  candidates <- all_fit_files[family_match & not_aux]

  if (length(candidates) == 0) {
    stop(
      "Could not find a candidate fit for family = ", family_tag, ".\n",
      "Searched folder: ", fits_dir, "\n",
      "Available .rds files:\n",
      paste(bn, collapse = "\n")
    )
  }

  cbn <- basename(candidates)

  score <- rep(0, length(candidates))

  for (i in seq_along(fit_search_terms)) {
    term <- fit_search_terms[[i]]
    score <- score + ifelse(grepl(term, cbn, ignore.case = TRUE), 100 - i, 0)
  }

  score <- score + ifelse(grepl("loggaussian|gaussian_logfit", cbn, ignore.case = TRUE), 10, 0)
  score <- score - ifelse(grepl("raw|poisson|negbin|nb|count", cbn, ignore.case = TRUE), 20, 0)

  finfo <- file.info(candidates)

  candidate_table <- tibble::tibble(
    file = candidates,
    basename = cbn,
    score = score,
    modified = finfo$mtime
  ) %>%
    dplyr::arrange(dplyr::desc(score), dplyr::desc(modified))

  selected <- candidate_table$file[[1]]

  cat("\nFit candidates for ", family_tag, ":\n", sep = "")
  print(candidate_table, n = Inf)
  cat("Selected fit for ", family_tag, ": ", selected, "\n", sep = "")

  normalizePath(selected, winslash = "/", mustWork = TRUE)
}


############################
## 3. Input/output folders
############################

project_dir <- find_project_dir(project_dir_override)

structure_dir <- file.path(project_dir, "01_alpha_structure_selection")
outcomes_dir  <- file.path(project_dir, "02_outcomes")
tables_dir    <- file.path(outcomes_dir, "tables")
draws_dir     <- file.path(outcomes_dir, "draws")
figures_dir   <- file.path(outcomes_dir, "figures")
summaries_dir <- file.path(outcomes_dir, "summaries")

dir.create(outcomes_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(draws_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(summaries_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(structure_dir)) stop("structure_dir does not exist: ", structure_dir)

invisible(lapply(c(tables_dir, draws_dir, figures_dir, summaries_dir), check_write_dir))

cat("\nProject dir:   ", safe_path_report(project_dir), "\n", sep = "")
cat("Structure dir: ", safe_path_report(structure_dir), "\n", sep = "")
cat("Outcomes dir:  ", safe_path_report(outcomes_dir), "\n", sep = "")
cat("Tables dir:    ", safe_path_report(tables_dir), "\n", sep = "")
cat("Figures dir:   ", safe_path_report(figures_dir), "\n", sep = "")
cat("Summaries dir: ", safe_path_report(summaries_dir), "\n\n", sep = "")


############################
## 4. Load alpha-selection metadata and prepared data path
############################

alpha_info_file <- file.path(structure_dir, "selection", "alpha_structure_selection_info.rds")

if (!file.exists(alpha_info_file)) {
  stop("Cannot find alpha_structure_selection_info.rds at: ", alpha_info_file)
}

alpha_info <- readRDS(alpha_info_file)

stage1_data_file <- resolve_existing_file(
  alpha_info$stage1_data_rds,
  fallback_dirs = c(
    file.path(structure_dir, "data"),
    file.path(structure_dir, "draws"),
    structure_dir,
    project_dir
  ),
  label = "stage1 data file",
  prefer_current = TRUE,
  search_pattern = "stage1.*\\.rds$|analysis_data.*\\.rds$"
)

selected_growth_family_from_alpha <- as.character(alpha_info$selected_growth_family)

error_tag <- ifelse(
  identical(alpha_info$selected_error_model, "gaussian_logfit"),
  "loggaussian",
  as.character(alpha_info$selected_error_model)
)

prefix_tag <- ifelse(use_germination, "germ", "nogerm")


############################
## 5. Run one model family
############################

run_one_family <- function(family_tag) {

  family_tag <- dplyr::case_when(
    toupper(family_tag) == "BH" ~ "BH",
    toupper(family_tag) == "RICKER" ~ "Ricker",
    TRUE ~ as.character(family_tag)
  )

  fit_override <- fit_file_overrides[[family_tag]]

  if (!is.na(fit_override) && nzchar(fit_override)) {
    fit_file <- resolve_existing_file(
      fit_override,
      fallback_dirs = c(
        file.path(structure_dir, "fits"),
        structure_dir,
        project_dir
      ),
      label = paste0("manual fit override for ", family_tag)
    )
  } else {
    fit_file <- find_fit_for_family(
      family_tag = family_tag,
      structure_dir = structure_dir,
      fit_search_terms = fit_search_terms
    )
  }

  dat_stage1 <- readRDS(stage1_data_file)
  fit_final  <- readRDS(fit_file)

  file_prefix <- make_file_prefix(family_tag, error_tag, use_germination)

  cat("\n============================================================\n")
  cat("Running family: ", family_tag, "\n", sep = "")
  cat("Stage-1 data:   ", stage1_data_file, "\n", sep = "")
  cat("Fit file:       ", fit_file, "\n", sep = "")
  cat("Output prefix:  ", file_prefix, "\n", sep = "")
  cat("============================================================\n\n")

############################
## 3a. Standardise/check data labels
############################

required_cols <- c("focsp", "Site", "Cov", "BLK", "fac_sc")
missing_cols <- setdiff(required_cols, names(dat_stage1))

if (length(missing_cols) > 0) {
  stop("dat_stage1 is missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat_stage1 <- dat_stage1 %>%
  dplyr::mutate(
    Site = toupper(as.character(Site)),
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

if (length(bad_sites) > 0) {
  stop("Unexpected Site values after standardisation: ", paste(bad_sites, collapse = ", "))
}
if (length(bad_cov) > 0) {
  stop("Unexpected Cov values after standardisation: ", paste(bad_cov, collapse = ", "))
}
if (length(bad_focsp) > 0) {
  stop("Unexpected focsp values after standardisation: ", paste(bad_focsp, collapse = ", "))
}

## Site order: use alpha_info if it matches; otherwise wet-to-dry default.
site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

if (!is.null(alpha_info$site_order)) {
  alpha_site_order <- toupper(as.character(alpha_info$site_order))
  if (setequal(alpha_site_order, unique(dat_stage1$Site))) {
    site_order <- alpha_site_order
  }
}

site_order <- site_order[site_order %in% unique(dat_stage1$Site)]

## Keep cover order as Shade, Sun for model consistency.
## Plotting sections can reverse order locally where needed.
cov_order <- c("Shade", "Sun")
cov_order <- cov_order[cov_order %in% unique(dat_stage1$Cov)]

if (length(site_order) == 0) stop("site_order is empty after checking dat_stage1$Site.")
if (length(cov_order) == 0)  stop("cov_order is empty after checking dat_stage1$Cov.")

dat_stage1 <- dat_stage1 %>%
  dplyr::mutate(
    Site = factor(Site, levels = site_order),
    Cov  = factor(Cov, levels = cov_order),
    focsp = factor(focsp, levels = allowed_focsp),
    SiteBlock = factor(SiteBlock)
  )

if (any(is.na(dat_stage1$fac_sc))) {
  stop("dat_stage1$fac_sc contains NA values. Fix these before running coexistence calculations.")
}

if (!any(abs(dat_stage1$fac_sc) < .Machine$double.eps^0.5, na.rm = TRUE)) {
  warning(
    "No fac_sc values equal to 0 were found. This script assumes fac_sc = 0 means no facilitator. ",
    "If fac_sc was centred or z-scaled, the no-facilitator predictions will be wrong."
  )
}


############################
## 3b. Define facilitator levels: none vs high
############################

selected_fac_levels <- bind_rows(
  dat_stage1 %>%
    dplyr::distinct(Site, Cov) %>%
    dplyr::mutate(
      fac_sc = 0,
      fac_level = "None",
      n_positive = NA_integer_
    ),
  dat_stage1 %>%
    dplyr::group_by(Site, Cov) %>%
    dplyr::summarise(
      n_positive = sum(.data$fac_sc > 0, na.rm = TRUE),
      fac_sc = {
        x <- .data$fac_sc[is.finite(.data$fac_sc) & .data$fac_sc > 0]
        if (length(x) == 0) {
          0
        } else {
          q <- stats::quantile(x, high_fac_quantile, na.rm = TRUE)
          mean(x[x >= q], na.rm = TRUE)
        }
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(fac_level = "High")
) %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_level = factor(fac_level, levels = c("None", "High"))
  )

selected_fac_levels_wide <- selected_fac_levels %>%
  dplyr::select(Site, Cov, fac_level, fac_sc) %>%
  tidyr::pivot_wider(
    names_from = fac_level,
    values_from = fac_sc
  ) %>%
  dplyr::left_join(
    selected_fac_levels %>%
      dplyr::filter(fac_level == "High") %>%
      dplyr::select(Site, Cov, n_positive),
    by = c("Site", "Cov")
  ) %>%
  dplyr::mutate(delta = High - None) %>%
  dplyr::arrange(Site, Cov)

print(selected_fac_levels_wide)


############################
## 3c. Define facilitator levels: none/low/high for optional summaries
############################

selected_fac_outcome_levels <- dat_stage1 %>%
  dplyr::group_by(Site, Cov) %>%
  dplyr::summarise(
    n_positive = sum(.data$fac_sc > 0, na.rm = TRUE),
    fac_none = 0,
    fac_low = {
      x <- .data$fac_sc[is.finite(.data$fac_sc) & .data$fac_sc > 0]
      if (length(x) == 0) 0 else as.numeric(stats::quantile(x, 0.25, na.rm = TRUE))
    },
    fac_high = {
      x <- .data$fac_sc[is.finite(.data$fac_sc) & .data$fac_sc > 0]
      if (length(x) == 0) {
        0
      } else {
        q <- stats::quantile(x, high_fac_quantile, na.rm = TRUE)
        mean(x[x >= q], na.rm = TRUE)
      }
    },
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(fac_none, fac_low, fac_high),
    names_to = "fac_key",
    values_to = "fac_sc"
  ) %>%
  dplyr::mutate(
    fac_label = dplyr::case_when(
      fac_key == "fac_none" ~ "No facilitator",
      fac_key == "fac_low"  ~ "Low facilitator",
      fac_key == "fac_high" ~ "High facilitator",
      TRUE ~ fac_key
    ),
    fac_label = factor(
      fac_label,
      levels = c("No facilitator", "Low facilitator", "High facilitator")
    ),
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order)
  ) %>%
  dplyr::select(Site, Cov, fac_label, fac_sc, n_positive)

selected_fac_outcome_levels_wide <- selected_fac_outcome_levels %>%
  dplyr::select(-n_positive) %>%
  tidyr::pivot_wider(
    names_from = fac_label,
    values_from = fac_sc
  ) %>%
  dplyr::arrange(Site, Cov)

print(selected_fac_outcome_levels_wide)


############################
## 4. Helper: matrix -> long data frame
############################

mat_to_long <- function(mat, meta, value_name) {
  meta <- meta %>% dplyr::mutate(row_id = dplyr::row_number())
  colnames(mat) <- paste0("V", seq_len(ncol(mat)))

  tibble::as_tibble(mat, .name_repair = "minimal") %>%
    dplyr::mutate(draw = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = -draw,
      names_to = "col",
      values_to = value_name
    ) %>%
    dplyr::mutate(row_id = as.integer(sub("^V", "", col))) %>%
    dplyr::select(-col) %>%
    dplyr::left_join(meta, by = "row_id") %>%
    dplyr::select(draw, dplyr::everything(), -row_id)
}


############################
## 5. Build prediction grids
############################

fac_values_all <- sort(unique(c(
  dat_stage1$fac_sc,
  selected_fac_levels$fac_sc,
  selected_fac_outcome_levels$fac_sc
)))

fac_values_all <- fac_values_all[is.finite(fac_values_all)]

if (length(fac_values_all) == 0) stop("fac_values_all is empty.")

newdata_defaults <- dat_stage1 %>%
  dplyr::slice(1) %>%
  dplyr::select(BLK, SiteBlock)

lambda_grid <- tidyr::expand_grid(
  focsp  = levels(dat_stage1$focsp),
  Site   = levels(dat_stage1$Site),
  Cov    = levels(dat_stage1$Cov),
  fac_sc = fac_values_all
) %>%
  dplyr::mutate(
    BLK       = newdata_defaults$BLK[1],
    SiteBlock = newdata_defaults$SiteBlock[1],
    intra_sc  = 0,
    inter_sc  = 0,
    focsp = factor(focsp, levels = levels(dat_stage1$focsp)),
    Site = factor(Site, levels = levels(dat_stage1$Site)),
    Cov = factor(Cov, levels = levels(dat_stage1$Cov)),
    BLK = factor(BLK, levels = levels(dat_stage1$BLK)),
    SiteBlock = factor(SiteBlock, levels = levels(dat_stage1$SiteBlock))
  )

alpha_grid <- tidyr::expand_grid(
  focsp = levels(dat_stage1$focsp),
  Site  = levels(dat_stage1$Site),
  Cov   = levels(dat_stage1$Cov)
) %>%
  dplyr::mutate(
    fac_sc    = 0,
    BLK       = newdata_defaults$BLK[1],
    SiteBlock = newdata_defaults$SiteBlock[1],
    intra_sc  = 0,
    inter_sc  = 0,
    focsp = factor(focsp, levels = levels(dat_stage1$focsp)),
    Site = factor(Site, levels = levels(dat_stage1$Site)),
    Cov = factor(Cov, levels = levels(dat_stage1$Cov)),
    BLK = factor(BLK, levels = levels(dat_stage1$BLK)),
    SiteBlock = factor(SiteBlock, levels = levels(dat_stage1$SiteBlock))
  )


############################
## 6. Extract posterior draws of nonlinear predictors
############################

etalam_draws <- brms::posterior_linpred(
  fit_final,
  newdata    = lambda_grid,
  nlpar      = "etalam",
  transform  = FALSE,
  re_formula = NA
)

etaai_draws <- brms::posterior_linpred(
  fit_final,
  newdata    = alpha_grid,
  nlpar      = "etaai",
  transform  = FALSE,
  re_formula = NA
)

etaaj_draws <- brms::posterior_linpred(
  fit_final,
  newdata    = alpha_grid,
  nlpar      = "etaaj",
  transform  = FALSE,
  re_formula = NA
)


############################
## 7. Convert to lambda and alpha on natural scale
############################

lambda_draws <- mat_to_long(etalam_draws, lambda_grid, "etalam") %>%
  dplyr::mutate(lambda = exp(etalam)) %>%
  dplyr::select(draw, focsp, Site, Cov, fac_sc, lambda)

alpha_ii_draws <- mat_to_long(etaai_draws, alpha_grid, "etaai") %>%
  dplyr::mutate(alpha_ii = exp(etaai)) %>%
  dplyr::select(draw, focsp, Site, Cov, alpha_ii)

## alpha_ij_focspTRCY = effect of TROR on TRCY
## alpha_ij_focspTROR = effect of TRCY on TROR
alpha_ij_draws <- mat_to_long(etaaj_draws, alpha_grid, "etaaj") %>%
  dplyr::mutate(alpha_ij = exp(etaaj)) %>%
  dplyr::select(draw, focsp, Site, Cov, alpha_ij)

alpha_draws <- alpha_ii_draws %>%
  dplyr::left_join(alpha_ij_draws, by = c("draw", "focsp", "Site", "Cov"))


############################
## 8. Save summaries of lambda and alpha
############################

lambda_summary <- lambda_draws %>%
  dplyr::group_by(focsp, Site, Cov, fac_sc) %>%
  dplyr::summarise(
    lambda_med = stats::median(lambda),
    lambda_lwr = stats::quantile(lambda, 0.025),
    lambda_upr = stats::quantile(lambda, 0.975),
    .groups = "drop"
  )

alpha_summary <- alpha_draws %>%
  dplyr::group_by(focsp, Site, Cov) %>%
  dplyr::summarise(
    alpha_ii_med = stats::median(alpha_ii),
    alpha_ii_lwr = stats::quantile(alpha_ii, 0.025),
    alpha_ii_upr = stats::quantile(alpha_ii, 0.975),
    alpha_ij_med = stats::median(alpha_ij),
    alpha_ij_lwr = stats::quantile(alpha_ij, 0.025),
    alpha_ij_upr = stats::quantile(alpha_ij, 0.975),
    .groups = "drop"
  )

safe_write_csv(lambda_summary, table_path(file_prefix, "lam"))
safe_write_csv(alpha_summary,  table_path(file_prefix, "alpha"))


############################
## 9. Germination draws OR no-germination default
############################

n_draws <- dplyr::n_distinct(lambda_draws$draw)

if (use_germination) {

  germ_site_priors_file <- NULL

  if (!is.null(alpha_info$germination_site_priors_rds) &&
      length(alpha_info$germination_site_priors_rds) == 1 &&
      !is.na(alpha_info$germination_site_priors_rds)) {

    ## Prefer the current project copy. Do not silently reuse old absolute paths
    ## stored inside alpha_info from previous runs.
    germ_site_priors_file <- tryCatch(
      resolve_existing_file(
        alpha_info$germination_site_priors_rds,
        fallback_dirs = c(file.path(structure_dir, "data"), structure_dir, project_dir),
        label = "germination site-priors file",
        prefer_current = TRUE,
        search_pattern = "germ.*site.*prior.*\\.rds$|site.*prior.*germ.*\\.rds$"
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(germ_site_priors_file)) {

    germ_priors <- readRDS(germ_site_priors_file) %>%
      tibble::as_tibble()

  } else {

    germ_info_file <- resolve_existing_file(
      alpha_info$germination_priors_rds,
      fallback_dirs = c(file.path(structure_dir, "data"), structure_dir, project_dir),
      label = "germination priors file"
    )

    germ_info <- readRDS(germ_info_file)
    germ_priors <- germ_info$site_priors %>%
      tibble::as_tibble()
  }

  germ_required <- c("Site", "focsp", "g_alpha_final", "g_beta_final", "g_mean_final")
  germ_missing <- setdiff(germ_required, names(germ_priors))

  if (length(germ_missing) > 0) {
    stop("germ_priors is missing columns: ", paste(germ_missing, collapse = ", "))
  }

  germ_priors <- germ_priors %>%
    dplyr::mutate(
      Site  = toupper(as.character(Site)),
      focsp = toupper(as.character(focsp))
    )

  set.seed(123)

  germ_draws <- tidyr::expand_grid(
    draw = seq_len(n_draws),
    germ_priors %>% dplyr::select(Site, focsp, g_alpha_final, g_beta_final, g_mean_final)
  ) %>%
    dplyr::mutate(
      g_draw = stats::rbeta(
        n = dplyr::n(),
        shape1 = g_alpha_final,
        shape2 = g_beta_final
      )
    )

  germ_summary <- germ_draws %>%
    dplyr::group_by(Site, focsp) %>%
    dplyr::summarise(
      g_input_mean = dplyr::first(g_mean_final),
      g_draw_med   = stats::median(g_draw),
      g_draw_lwr   = stats::quantile(g_draw, 0.025),
      g_draw_upr   = stats::quantile(g_draw, 0.975),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Site = factor(Site, levels = site_order)) %>%
    dplyr::arrange(focsp, Site)

  germ_wide <- germ_draws %>%
    dplyr::select(draw, Site, focsp, g_draw) %>%
    tidyr::pivot_wider(
      names_from  = focsp,
      values_from = g_draw,
      names_glue  = "g_{focsp}"
    )

} else {

  germ_summary <- tidyr::expand_grid(
    Site = levels(dat_stage1$Site),
    focsp = levels(dat_stage1$focsp)
  ) %>%
    dplyr::mutate(
      g_input_mean = 1,
      g_draw_med   = 1,
      g_draw_lwr   = 1,
      g_draw_upr   = 1
    ) %>%
    dplyr::mutate(Site = factor(Site, levels = site_order)) %>%
    dplyr::arrange(focsp, Site)

  germ_wide <- tidyr::expand_grid(
    draw = seq_len(n_draws),
    Site = levels(dat_stage1$Site)
  ) %>%
    dplyr::mutate(
      g_TRCY = 1,
      g_TROR = 1
    )
}

safe_write_csv(germ_summary, table_path(file_prefix, "germ"))
safe_save_rds(germ_wide, draw_path(file_prefix, "germ"))


############################
## 10. Wide format for coexistence calculations
############################

lambda_wide <- lambda_draws %>%
  tidyr::pivot_wider(
    names_from  = focsp,
    values_from = lambda,
    names_glue  = "lambda_{focsp}"
  )

alpha_wide <- alpha_draws %>%
  tidyr::pivot_wider(
    names_from  = focsp,
    values_from = c(alpha_ii, alpha_ij),
    names_glue  = "{.value}_{focsp}"
  )

coex_draws <- lambda_wide %>%
  dplyr::left_join(alpha_wide, by = c("draw", "Site", "Cov")) %>%
  dplyr::left_join(germ_wide,  by = c("draw", "Site"))

coex_required <- c(
  "lambda_TRCY", "lambda_TROR", "g_TRCY", "g_TROR",
  "alpha_ii_TRCY", "alpha_ii_TROR", "alpha_ij_TRCY", "alpha_ij_TROR"
)

coex_missing <- setdiff(coex_required, names(coex_draws))

if (length(coex_missing) > 0) {
  stop("coex_draws is missing columns after joining: ", paste(coex_missing, collapse = ", "))
}


############################
## 11. Germination / seed survival
############################

## No seed bank in either version.
s_vals <- c(TROR = 0, TRCY = 0)

coex_draws <- coex_draws %>%
  dplyr::mutate(
    eta_TROR = lambda_TROR * g_TROR /
      (1 - s_vals["TROR"] * (1 - g_TROR)),
    eta_TRCY = lambda_TRCY * g_TRCY /
      (1 - s_vals["TRCY"] * (1 - g_TRCY))
  )


############################
## 12. Niche and fitness differences
############################

## Define the density-independent growth term for the fitness ratio.
## BH:     r_i = eta_i - 1
## Ricker: r_i = log(eta_i)

if (family_tag == "BH") {

  coex_draws <- coex_draws %>%
    dplyr::mutate(
      r_TROR = eta_TROR - 1,
      r_TRCY = eta_TRCY - 1
    )

} else if (family_tag == "Ricker") {

  coex_draws <- coex_draws %>%
    dplyr::mutate(
      r_TROR = log(eta_TROR),
      r_TRCY = log(eta_TRCY)
    )

} else {
  stop("Unknown growth family: ", family_tag)
}

coex_draws <- coex_draws %>%
  dplyr::mutate(
    rho = sqrt(
      (alpha_ij_TROR * alpha_ij_TRCY) /
        (alpha_ii_TROR * alpha_ii_TRCY)
    ),

    niche_diff = 1 - rho,

    valid_kappa_base =
      is.finite(r_TRCY) &
      is.finite(r_TROR) &
      r_TROR != 0 &
      is.finite(rho) &
      rho > 0 &
      alpha_ij_TROR > 0 &
      alpha_ii_TROR > 0 &
      alpha_ii_TRCY > 0 &
      alpha_ij_TRCY > 0
  )

if (restrict_kappa_to_positive_growth) {
  coex_draws <- coex_draws %>%
    dplyr::mutate(
      valid_kappa = valid_kappa_base & r_TRCY > 0 & r_TROR > 0
    )
} else {
  coex_draws <- coex_draws %>%
    dplyr::mutate(
      valid_kappa = valid_kappa_base
    )
}

coex_draws <- coex_draws %>%
  dplyr::mutate(
    fitness_ratio_TRCY_over_TROR = dplyr::if_else(
      valid_kappa,
      (r_TRCY / r_TROR) *
        sqrt(
          (alpha_ij_TROR * alpha_ii_TROR) /
            (alpha_ii_TRCY * alpha_ij_TRCY)
        ),
      NA_real_
    ),

    inv_rho = 1 / rho,

    lower_bound = pmin(rho, inv_rho),
    upper_bound = pmax(rho, inv_rho),

    coexist = valid_kappa &
      rho < 1 &
      fitness_ratio_TRCY_over_TROR > lower_bound &
      fitness_ratio_TRCY_over_TROR < upper_bound,

    priority_effect = valid_kappa &
      rho > 1 &
      fitness_ratio_TRCY_over_TROR > lower_bound &
      fitness_ratio_TRCY_over_TROR < upper_bound,

    TRCY_wins = valid_kappa &
      fitness_ratio_TRCY_over_TROR > upper_bound,

    TROR_wins = valid_kappa &
      fitness_ratio_TRCY_over_TROR < lower_bound
  )

eta_validity_summary <- coex_draws %>%
  dplyr::summarise(
    p_eta_TRCY_gt1 = mean(eta_TRCY > 1, na.rm = TRUE),
    p_eta_TROR_gt1 = mean(eta_TROR > 1, na.rm = TRUE),
    p_both_eta_gt1 = mean(eta_TRCY > 1 & eta_TROR > 1, na.rm = TRUE),
    p_valid_kappa = mean(valid_kappa, na.rm = TRUE)
  )

outcome_sum_check <- coex_draws %>%
  dplyr::mutate(
    outcome_sum =
      as.numeric(coexist) +
      as.numeric(TRCY_wins) +
      as.numeric(TROR_wins) +
      as.numeric(priority_effect)
  ) %>%
  dplyr::count(outcome_sum)

print(eta_validity_summary)
print(outcome_sum_check)


############################
## 12b. Stacked posterior outcome probabilities by cover
############################

fac_join_digits <- 10

cov_sun_value <- intersect(c("Sun", "SUN", "sun"), cov_order)[1]
cov_shade_value <- intersect(c("Shade", "SH", "shade", "SHADE"), cov_order)[1]

if (is.na(cov_sun_value) || is.na(cov_shade_value)) {
  stop(
    "Could not identify Sun/Shade cover levels from cov_order. cov_order is: ",
    paste(cov_order, collapse = ", ")
  )
}

## Top row = Sun, bottom row = Shade for the stacked outcome plot.
cover_plot_levels <- c(cov_sun_value, cov_shade_value)

selected_fac_outcome_join <- selected_fac_levels %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits),
    fac_present = dplyr::case_when(
      fac_level == "None" ~ "No",
      fac_level == "High" ~ "Yes",
      TRUE ~ as.character(fac_level)
    ),
    fac_present = factor(fac_present, levels = c("No", "Yes"))
  ) %>%
  dplyr::select(Site, Cov, fac_present, fac_sc, fac_sc_join)

coex_outcome_draws_2fac_cover <- coex_draws %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
  ) %>%
  dplyr::inner_join(
    selected_fac_outcome_join,
    by = c("Site", "Cov", "fac_sc_join")
  )

if (nrow(coex_outcome_draws_2fac_cover) == 0) {
  stop(
    "coex_outcome_draws_2fac_cover has zero rows. Check that selected_fac_levels$fac_sc ",
    "was included in fac_values_all before lambda_grid and coex_draws were created."
  )
}

outcome_bar_summary_cover <- coex_outcome_draws_2fac_cover %>%
  dplyr::group_by(Site, Cov, fac_present) %>%
  dplyr::summarise(
    p_TRCY_wins = mean(TRCY_wins, na.rm = TRUE),
    p_TROR_wins = mean(TROR_wins, na.rm = TRUE),
    p_coexist = mean(coexist, na.rm = TRUE),
    p_priority = mean(priority_effect, na.rm = TRUE),
    p_valid_kappa = mean(valid_kappa, na.rm = TRUE),
    n_rows = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_sum = p_TRCY_wins + p_TROR_wins + p_coexist + p_priority,
    p_unclassified = 1 - p_sum
  ) %>%
  dplyr::arrange(Cov, Site, fac_present)

print(outcome_bar_summary_cover)

safe_write_csv(outcome_bar_summary_cover, table_path(file_prefix, "stack"))

outcome_bar_long_cover <- outcome_bar_summary_cover %>%
  dplyr::select(
    Site, Cov, fac_present,
    p_TRCY_wins, p_TROR_wins, p_coexist, p_priority
  ) %>%
  tidyr::pivot_longer(
    cols = c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority),
    names_to = "outcome",
    values_to = "probability"
  ) %>%
  dplyr::mutate(
    outcome = dplyr::case_when(
      outcome == "p_TRCY_wins" ~ "TRCY wins",
      outcome == "p_TROR_wins" ~ "TROR wins",
      outcome == "p_coexist" ~ "Coexistence",
      outcome == "p_priority" ~ "Priority effect",
      TRUE ~ outcome
    ),
    outcome = factor(
      outcome,
      levels = c("TRCY wins", "TROR wins", "Coexistence", "Priority effect")
    ),
    Site = factor(as.character(Site), levels = site_order),
    Cov = factor(as.character(Cov), levels = cover_plot_levels),
    cover_plot = factor(
      as.character(Cov),
      levels = cover_plot_levels,
      labels = c("Sun", "Shade")
    ),
    fac_present = factor(fac_present, levels = c("No", "Yes"))
  )

outcome_cols <- c(
  "TRCY wins" = "#83E8BA",
  "TROR wins" = "grey",
  "Coexistence" = "#72A1E5",
  "Priority effect" = "#9883E5"
)

p_outcome_stack_fac_present_cover <- ggplot(
  outcome_bar_long_cover,
  aes(x = fac_present, y = probability, fill = outcome)
) +
  geom_col(
    width = 0.78,
    colour = "white",
    linewidth = 0.25
  ) +
  facet_grid(cover_plot ~ Site) +
  scale_fill_manual(
    values = outcome_cols,
    name = "Posterior outcome"
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Facilitators present",
    y = "Posterior probability"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(0.35, "lines"),
    panel.spacing.y = grid::unit(0.25, "lines"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5)
  )

print(p_outcome_stack_fac_present_cover)

safe_save_plot(
  p_outcome_stack_fac_present_cover,
  fig_path(file_prefix, "stack"),
  width = 8.5,
  height = 4.8,
  dpi = save_png_dpi_main
)

############################
## 13. Restrict to observed environments only
############################

obs_env <- dat_stage1 %>%
  dplyr::distinct(Site, Cov, fac_sc)

coex_draws_obs <- coex_draws %>%
  dplyr::semi_join(obs_env, by = c("Site", "Cov", "fac_sc"))

env_counts <- dat_stage1 %>%
  dplyr::count(Site, Cov, fac_sc, name = "n_obs")


############################
## 14. Summarise propagated uncertainty by environment
############################

coex_summary <- coex_draws_obs %>%
  dplyr::group_by(Site, Cov, fac_sc) %>%
  dplyr::summarise(
    niche_med = stats::median(niche_diff, na.rm = TRUE),
    niche_lwr = stats::quantile(niche_diff, 0.025, na.rm = TRUE),
    niche_upr = stats::quantile(niche_diff, 0.975, na.rm = TRUE),

    fit_med = stats::median(fitness_ratio_TRCY_over_TROR, na.rm = TRUE),
    fit_lwr = stats::quantile(fitness_ratio_TRCY_over_TROR, 0.025, na.rm = TRUE),
    fit_upr = stats::quantile(fitness_ratio_TRCY_over_TROR, 0.975, na.rm = TRUE),

    p_coexist = mean(coexist, na.rm = TRUE),
    p_TRCY_wins = mean(TRCY_wins, na.rm = TRUE),
    p_TROR_wins = mean(TROR_wins, na.rm = TRUE),
    p_priority = mean(priority_effect, na.rm = TRUE),

    mean_g_TROR = mean(g_TROR, na.rm = TRUE),
    mean_g_TRCY = mean(g_TRCY, na.rm = TRUE),
    p_valid_kappa = mean(valid_kappa, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::left_join(env_counts, by = c("Site", "Cov", "fac_sc")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p_sum = sum(c(p_coexist, p_TRCY_wins, p_TROR_wins, p_priority), na.rm = TRUE),
    p_unclassified = 1 - p_sum,
    max_support = max(c(p_coexist, p_TRCY_wins, p_TROR_wins, p_priority), na.rm = TRUE),
    most_supported_outcome = c("coexist", "TRCY_wins", "TROR_wins", "priority_effect")[
      which.max(c(p_coexist, p_TRCY_wins, p_TROR_wins, p_priority))
    ],
    outcome_90 = ifelse(max_support >= 0.90, most_supported_outcome, "ambiguous")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order)
  )

safe_write_csv(coex_summary, table_path(file_prefix, "coex_obs"))

safe_save_rds(lambda_draws,   draw_path(file_prefix, "lambda"))
safe_save_rds(alpha_draws,    draw_path(file_prefix, "alpha"))
safe_save_rds(coex_draws_obs, draw_path(file_prefix, "coex_obs"))
safe_save_rds(coex_draws, draw_path(file_prefix, "coex_all"))

############################
## 15. Quick checks
############################

capture.output(
  {
    cat("\nUsing germination:\n")
    print(use_germination)

    cat("\nGrowth family:\n")
    print(family_tag)

    cat("\nSelected family from alpha_info:\n")
    print(selected_growth_family_from_alpha)

    cat("\nFamily run mode:\n")
    print("automatic two-family run")

    cat("\nSelected error model:\n")
    print(alpha_info$selected_error_model)

    cat("\nStage-1 data file:\n")
    print(stage1_data_file)

    cat("\nFit file:\n")
    print(fit_file)

    cat("\nProject dir:\n")
    print(project_dir)

    cat("\nFigures dir:\n")
    print(figures_dir)

    cat("\nSite order used:\n")
    print(site_order)

    cat("\nCov order used:\n")
    print(cov_order)

    cat("\nSelected facilitator levels:\n")
    print(selected_fac_levels_wide)

    cat("\nEta validity summary:\n")
    print(eta_validity_summary)

    cat("\nOutcome-sum check:\n")
    print(outcome_sum_check)

    cat("\nOutcome-90 counts:\n")
    print(coex_summary %>% dplyr::count(outcome_90))

    cat("\nMost-supported-outcome counts:\n")
    print(coex_summary %>% dplyr::count(most_supported_outcome))

    cat("\nSummary of max support:\n")
    print(summary(coex_summary$max_support))
  },
  file = summary_path(file_prefix, "checks", ext = "txt")
)


############################
## 16. LDGR distributions across sites
############################

fac_join_digits <- 10

panel_row_levels <- c(
  paste(cov_sun_value, "None", sep = "-"),
  paste(cov_sun_value, "High", sep = "-"),
  paste(cov_shade_value, "None", sep = "-"),
  paste(cov_shade_value, "High", sep = "-")
)

selected_env_ldgr <- selected_fac_levels %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits),
    panel_row = paste(Cov, fac_level, sep = "-"),
    panel_row = factor(panel_row, levels = panel_row_levels)
  )

ldgr_base <- coex_draws %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
  ) %>%
  dplyr::inner_join(
    selected_env_ldgr %>%
      dplyr::select(Site, Cov, fac_sc_join, fac_level, panel_row),
    by = c("Site", "Cov", "fac_sc_join")
  )

if (nrow(ldgr_base) == 0) {
  debug_selected <- selected_env_ldgr %>%
    dplyr::select(Site, Cov, fac_level, fac_sc, fac_sc_join) %>%
    dplyr::arrange(Site, Cov, fac_level)

  debug_coex <- coex_draws %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    ) %>%
    dplyr::distinct(Site, Cov, fac_sc, fac_sc_join) %>%
    dplyr::arrange(Site, Cov, fac_sc)

  safe_write_csv(debug_selected, summary_path(file_prefix, "dbg_ldgr_env"))
  safe_write_csv(debug_coex,    summary_path(file_prefix, "dbg_coex_env"))

  stop(
    "ldgr_base has zero rows. Debug CSVs written to 02_outcomes/summaries. ",
    "Most likely causes are mismatched Site/Cov labels or fac_sc values not present in lambda_grid/coex_draws."
  )
}

ldgr_draws <- ldgr_base %>%
  dplyr::mutate(
    Nstar_TROR = dplyr::case_when(
      family_tag == "BH" & eta_TROR > 1 & g_TROR > 0 ~
        (eta_TROR - 1) / (alpha_ii_TROR * g_TROR),
      family_tag == "Ricker" & eta_TROR > 1 & g_TROR > 0 ~
        log(eta_TROR) / (alpha_ii_TROR * g_TROR),
      TRUE ~ 0
    ),
    Nstar_TRCY = dplyr::case_when(
      family_tag == "BH" & eta_TRCY > 1 & g_TRCY > 0 ~
        (eta_TRCY - 1) / (alpha_ii_TRCY * g_TRCY),
      family_tag == "Ricker" & eta_TRCY > 1 & g_TRCY > 0 ~
        log(eta_TRCY) / (alpha_ii_TRCY * g_TRCY),
      TRUE ~ 0
    ),

    ## Invasion growth when rare: omit invader self-competition.
    ## TROR invades resident TRCY at TRCY equilibrium.
    ldgr_TROR = dplyr::case_when(
      family_tag == "BH" ~ log(
        (1 - g_TROR) * s_vals["TROR"] +
          (g_TROR * lambda_TROR) /
          (1 + alpha_ij_TROR * g_TRCY * Nstar_TRCY)
      ),
      family_tag == "Ricker" ~ log(
        (1 - g_TROR) * s_vals["TROR"] +
          (g_TROR * lambda_TROR) *
          exp(-alpha_ij_TROR * g_TRCY * Nstar_TRCY)
      ),
      TRUE ~ NA_real_
    ),

    ## TRCY invades resident TROR at TROR equilibrium.
    ldgr_TRCY = dplyr::case_when(
      family_tag == "BH" ~ log(
        (1 - g_TRCY) * s_vals["TRCY"] +
          (g_TRCY * lambda_TRCY) /
          (1 + alpha_ij_TRCY * g_TROR * Nstar_TROR)
      ),
      family_tag == "Ricker" ~ log(
        (1 - g_TRCY) * s_vals["TRCY"] +
          (g_TRCY * lambda_TRCY) *
          exp(-alpha_ij_TRCY * g_TROR * Nstar_TROR)
      ),
      TRUE ~ NA_real_
    )
  )

ldgr_long <- ldgr_draws %>%
  dplyr::select(draw, Site, Cov, fac_sc, fac_level, panel_row, ldgr_TROR, ldgr_TRCY) %>%
  tidyr::pivot_longer(
    cols = c(ldgr_TROR, ldgr_TRCY),
    names_to = "species",
    values_to = "ldgr"
  ) %>%
  dplyr::mutate(
    species = dplyr::recode(
      species,
      ldgr_TROR = "TROR",
      ldgr_TRCY = "TRCY"
    ),
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    panel_row = factor(panel_row, levels = panel_row_levels)
  ) %>%
  dplyr::filter(is.finite(ldgr))

if (nrow(ldgr_long) == 0) {
  stop("ldgr_long has zero finite LDGR rows. Check lambda, alpha, eta, and germination draws.")
}

ldgr_points <- ldgr_long %>%
  dplyr::group_by(Site, panel_row, species) %>%
  dplyr::summarise(ldgr_point = stats::median(ldgr), .groups = "drop") %>%
  dplyr::mutate(Site = factor(as.character(Site), levels = site_order))

p_coexist_ldgr <- ldgr_draws %>%
  dplyr::group_by(Site, panel_row) %>%
  dplyr::summarise(
    p_coexist = mean(ldgr_TROR > 0 & ldgr_TRCY > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    prob_lab = paste0(round(100 * p_coexist), "%"),
    Site = factor(as.character(Site), levels = site_order),
    panel_row = factor(panel_row, levels = panel_row_levels)
  )

ldgr_xpos <- ldgr_long %>%
  dplyr::group_by(Site, panel_row) %>%
  dplyr::summarise(x_lab = stats::quantile(ldgr, 0.005, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    panel_row = factor(panel_row, levels = panel_row_levels)
  )

prob_plot <- p_coexist_ldgr %>%
  dplyr::left_join(ldgr_xpos, by = c("Site", "panel_row")) %>%
  dplyr::mutate(Site = factor(as.character(Site), levels = site_order))

species_cols <- c("TRCY" = "#83E8BA", "TROR" = "grey")

p_ldgr <- ggplot(ldgr_long, aes(x = ldgr, colour = species)) +
  geom_density(linewidth = 1.0, adjust = 1.2, alpha = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_point(
    data = ldgr_points,
    aes(x = ldgr_point, y = 0, colour = species),
    inherit.aes = FALSE,
    size = 1.5,
    alpha = 0.8
  ) +
  geom_text(
    data = prob_plot,
    aes(x = x_lab, y = Inf, label = prob_lab),
    inherit.aes = FALSE,
    hjust = 0.7,
    vjust = 1.4,
    size = 4,
    colour = "black"
  ) +
  facet_grid(panel_row ~ Site, scales = "free_y") +
  scale_colour_manual(values = species_cols, drop = FALSE) +
  scale_x_continuous(limits = c(-5, 5), breaks = c(-5, -3, 0, 3, 5)) +
  labs(
    x = "Low-density growth rate (LDGR)",
    y = "Density",
    colour = "Species"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "top",
    axis.title.x = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 11, colour = "black"),
    axis.text.x  = element_text(size = 11, colour = "black"),
    axis.text.y  = element_text(size = 11, colour = "black"),
    strip.text   = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    legend.text  = element_text(size = 11, colour = "black")
  ) +
  guides(alpha = "none")

safe_save_plot(
  p_ldgr,
  fig_path(file_prefix, "ldgr"),
  width = 8,
  height = 7,
  dpi = save_png_dpi_main
)

print(p_ldgr)


############################
## 17. Save LDGR tables
############################

ldgr_summary <- ldgr_long %>%
  dplyr::group_by(Site, Cov, fac_level, panel_row, species) %>%
  dplyr::summarise(
    ldgr_med = stats::median(ldgr, na.rm = TRUE),
    ldgr_lwr = stats::quantile(ldgr, 0.025, na.rm = TRUE),
    ldgr_upr = stats::quantile(ldgr, 0.975, na.rm = TRUE),
    p_ldgr_positive = mean(ldgr > 0, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(ldgr_summary,   table_path(file_prefix, "ldgr_sum"))
safe_write_csv(p_coexist_ldgr, table_path(file_prefix, "ldgr_pcoex"))
safe_save_rds(ldgr_draws,       draw_path(file_prefix, "ldgr"))


############################
## 18. Coexistence contour helpers
############################

x_min <- -1.0
x_max <- 1.0
y_min <- 0.05
y_max <- 7.00

y_min_log <- log10(y_min)
y_max_log <- log10(y_max)

make_density_grid <- function(dat, kde_n = 180) {

  ## This follows the plotting logic from the older BH script that produced
  ## the working contour panels. It is intentionally less aggressive than the
  ## previous safety-filtered version, because that version could create empty
  ## density layers and trigger ggplot/ggnewscale rendering errors.
  dat <- dat %>%
    dplyr::filter(
      is.finite(niche_diff),
      is.finite(fitness_ratio_log10)
    )

  if (nrow(dat) < 2) {
    return(tibble::tibble(
      niche_diff = numeric(),
      fitness_ratio_log10 = numeric(),
      density = numeric(),
      ndensity = numeric()
    ))
  }

  ## kde2d needs non-zero spread. Add a tiny deterministic perturbation only
  ## when a group is numerically degenerate.
  sdx <- stats::sd(dat$niche_diff, na.rm = TRUE)
  sdy <- stats::sd(dat$fitness_ratio_log10, na.rm = TRUE)

  if (!is.finite(sdx) || sdx == 0) {
    dat$niche_diff <- dat$niche_diff + seq(-1e-8, 1e-8, length.out = nrow(dat))
  }

  if (!is.finite(sdy) || sdy == 0) {
    dat$fitness_ratio_log10 <- dat$fitness_ratio_log10 + seq(-1e-8, 1e-8, length.out = nrow(dat))
  }

  kd <- MASS::kde2d(
    x = dat$niche_diff,
    y = dat$fitness_ratio_log10,
    n = kde_n,
    lims = c(x_min, x_max, y_min_log, y_max_log)
  )

  grid <- expand.grid(
    niche_diff = kd$x,
    fitness_ratio_log10 = kd$y
  )

  grid$density <- as.vector(kd$z)
  max_density <- max(grid$density, na.rm = TRUE)

  if (!is.finite(max_density) || max_density <= 0) {
    grid$ndensity <- NA_real_
  } else {
    grid$ndensity <- grid$density / max_density
  }

  tibble::as_tibble(grid)
}

build_contour_density_data <- function(contour_plot_dat, group_vars, kde_n = 180) {

  out <- contour_plot_dat %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::group_modify(~ make_density_grid(.x, kde_n = kde_n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      density_level = dplyr::case_when(
        ndensity > 6 / 7 ~ "6",
        ndensity > 5 / 7 ~ "5",
        ndensity > 4 / 7 ~ "4",
        ndensity > 3 / 7 ~ "3",
        ndensity > 2 / 7 ~ "2",
        ndensity > 1 / 7 ~ "1",
        TRUE ~ NA_character_
      ),
      density_level = factor(density_level, levels = as.character(1:6))
    ) %>%
    dplyr::filter(!is.na(density_level))

  ## If a group has data but no grid cells above the retained threshold because
  ## of numerical edge cases, retain the highest-density cells as level 6. This
  ## avoids empty ggnewscale layers while preserving the visual intent.
  if (nrow(out) == 0 && nrow(contour_plot_dat) > 0) {
    all_grid <- contour_plot_dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::group_modify(~ make_density_grid(.x, kde_n = kde_n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is.finite(density)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::slice_max(order_by = density, n = max(1, floor(kde_n / 10)), with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(density_level = factor("6", levels = as.character(1:6)))

    out <- all_grid
  }

  out
}

boundary_x_max <- min(x_max, 0.999)

boundary_pos <- tibble::tibble(
  niche_diff = seq(0.0001, boundary_x_max, length.out = 500)
) %>%
  dplyr::mutate(
    rho = 1 - niche_diff,
    lower = log10(rho),
    upper = log10(1 / rho)
  )

boundary_neg <- tibble::tibble(
  niche_diff = seq(x_min, -0.0001, length.out = 300)
) %>%
  dplyr::mutate(
    rho = 1 - niche_diff,
    lower = log10(1 / rho),
    upper = log10(rho)
  )

background_coexist <- boundary_pos %>%
  dplyr::mutate(
    ymin = lower,
    ymax = upper
  )

region_labels <- tibble::tibble(
  niche_diff = c(0.39, 0.39, 0.57, -0.70),
  fitness_ratio_log10 = log10(c(5.5, 0.14, 1.0, 1.0)),
  label = c("TRCY wins", "TROR wins", "Coexistence", "Priority\neffect")
)

expand_contour_layer <- function(layer_dat, facet_dat) {
  ## Ensures non-density layers carry the same facet variables as the plot.
  ## This prevents facet_grid() from failing when a KDE density layer is empty.
  tidyr::crossing(facet_dat, layer_dat)
}

kde_n <- 180

blue_level_cols <- c(
  "6" = "#08306B",
  "5" = "#08519C",
  "4" = "#2C7FB8",
  "3" = "#7FB3D5",
  "2" = "#BFDDF0",
  "1" = "#E5F1FA"
)

orange_level_cols <- c(
  "6" = "#7F2704",
  "5" = "#A63603",
  "4" = "#D94801",
  "3" = "#F18805",
  "2" = "#F6B35A",
  "1" = "#FDE6C9"
)


############################
## 19. Coexistence contour diagrams: combined cover
############################

fac_join_digits <- 10

selected_fac_contour_join <- selected_fac_levels %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits),
    fac_label = dplyr::case_when(
      fac_level == "None" ~ "- Facilitator",
      fac_level == "High" ~ "+ Facilitator",
      TRUE ~ as.character(fac_level)
    ),
    fac_label = factor(fac_label, levels = c("- Facilitator", "+ Facilitator"))
  ) %>%
  dplyr::select(Site, Cov, fac_label, fac_sc, fac_sc_join) %>%
  dplyr::distinct()

coex_contour_draws <- coex_draws %>%
  dplyr::mutate(
    Site = factor(as.character(Site), levels = site_order),
    Cov  = factor(as.character(Cov), levels = cov_order),
    fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
  ) %>%
  dplyr::inner_join(
    selected_fac_contour_join,
    by = c("Site", "Cov", "fac_sc_join"),
    relationship = "many-to-many"
  )

if (nrow(coex_contour_draws) == 0) {
  stop(
    "coex_contour_draws has zero rows. Check that selected_fac_levels$fac_sc ",
    "was included in fac_values_all before lambda_grid and coex_draws were created."
  )
}

contour_plot_dat <- coex_contour_draws %>%
  dplyr::filter(
    is.finite(niche_diff),
    is.finite(fitness_ratio_TRCY_over_TROR),
    fitness_ratio_TRCY_over_TROR > 0
  ) %>%
  dplyr::mutate(
    fitness_ratio_log10 = log10(fitness_ratio_TRCY_over_TROR),
    Site = factor(as.character(Site), levels = site_order),
    Cov = factor(as.character(Cov), levels = cov_order),
    fac_label = factor(fac_label, levels = c("- Facilitator", "+ Facilitator"))
  ) %>%
  dplyr::filter(
    niche_diff >= x_min,
    niche_diff <= x_max,
    fitness_ratio_log10 >= y_min_log,
    fitness_ratio_log10 <= y_max_log
  )

if (nrow(contour_plot_dat) == 0) {
  stop("contour_plot_dat has zero rows after filtering to plot limits.")
}

set.seed(123)
max_draws_per_site_fac <- 80000

contour_plot_dat_combined <- contour_plot_dat %>%
  dplyr::group_by(Site, fac_label) %>%
  dplyr::group_modify(~ {
    n_sample <- min(nrow(.x), max_draws_per_site_fac)
    dplyr::slice_sample(.x, n = n_sample)
  }) %>%
  dplyr::ungroup()

contour_group_check <- contour_plot_dat_combined %>%
  dplyr::count(Site, fac_label, name = "n_draws")

print(contour_group_check, n = Inf)

if (any(contour_group_check$n_draws < 50)) {
  warning(
    "Some Site x facilitator groups have fewer than 50 posterior draws after filtering. ",
    "Contours may be unstable or missing."
  )
}

density_plot_dat <- build_contour_density_data(
  contour_plot_dat_combined,
  group_vars = c("Site", "fac_label"),
  kde_n = kde_n
)

if (nrow(density_plot_dat) == 0) {
  warning(
    "density_plot_dat has zero rows for the combined-cover contour plot. ",
    "The background and boundary regions will still be plotted. ",
    "Check contour_group_check and the x/y limits if you expected density shading."
  )
}

contour_facets_combined <- tidyr::expand_grid(
  fac_label = factor(c("- Facilitator", "+ Facilitator"), levels = c("- Facilitator", "+ Facilitator")),
  Site = factor(site_order, levels = site_order)
)

background_coexist_combined <- expand_contour_layer(background_coexist, contour_facets_combined)
boundary_pos_combined      <- expand_contour_layer(boundary_pos, contour_facets_combined)
boundary_neg_combined      <- expand_contour_layer(boundary_neg, contour_facets_combined)
region_labels_combined    <- expand_contour_layer(region_labels, contour_facets_combined)

p_coex_contours_combined <- ggplot() +
  geom_blank(
    data = contour_facets_combined,
    aes(x = 0, y = 0),
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = background_coexist_combined,
    aes(x = niche_diff, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    colour = NA
  ) +
  geom_raster(
    data = density_plot_dat %>% dplyr::filter(fac_label == "- Facilitator"),
    aes(x = niche_diff, y = fitness_ratio_log10, fill = density_level),
    interpolate = FALSE
  ) +
  scale_fill_manual(
    name = "Levels",
    values = blue_level_cols,
    breaks = as.character(6:1),
    labels = as.character(6:1),
    drop = FALSE,
    guide = guide_legend(title.position = "top", order = 1, override.aes = list(colour = NA))
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = density_plot_dat %>% dplyr::filter(fac_label == "+ Facilitator"),
    aes(x = niche_diff, y = fitness_ratio_log10, fill = density_level),
    interpolate = FALSE
  ) +
  scale_fill_manual(
    name = "Levels",
    values = orange_level_cols,
    breaks = as.character(6:1),
    labels = as.character(6:1),
    drop = FALSE,
    guide = guide_legend(title.position = "top", order = 2, override.aes = list(colour = NA))
  ) +
  geom_line(data = boundary_pos_combined, aes(x = niche_diff, y = lower), inherit.aes = FALSE, colour = "red", linewidth = 0.5) +
  geom_line(data = boundary_pos_combined, aes(x = niche_diff, y = upper), inherit.aes = FALSE, colour = "red", linewidth = 0.5) +
  geom_line(data = boundary_neg_combined, aes(x = niche_diff, y = lower), inherit.aes = FALSE, colour = "grey40", linewidth = 0.5) +
  geom_line(data = boundary_neg_combined, aes(x = niche_diff, y = upper), inherit.aes = FALSE, colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey25") +
  geom_text(
    data = region_labels_combined,
    aes(x = niche_diff, y = fitness_ratio_log10, label = label),
    inherit.aes = FALSE,
    size = 2.6,
    colour = "black",
    fontface = "bold",
    lineheight = 0.9
  ) +
  facet_grid(fac_label ~ Site) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = c(-1.0, -0.50, 0, 0.5, 0.8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(y_min_log, y_max_log),
    breaks = log10(c(0.1, 0.5, 1, 3, 5)),
    labels = c("0.1", "0.5", "1.0", "3.0", "5.0"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = expression("Niche difference: " ~ 1 - rho),
    y = expression("Fitness difference: " ~ kappa[j] / kappa[i])
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.35, fill = NA),
    panel.spacing.x = grid::unit(0.35, "lines"),
    panel.spacing.y = grid::unit(0.35, "lines"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 9.5, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 8.5, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 8.5, colour = "black"),
    legend.text = element_text(size = 7.5, colour = "black"),
    legend.key.height = grid::unit(0.8, "lines"),
    legend.spacing.y = grid::unit(0.15, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

print(p_coex_contours_combined)

safe_save_plot(
  p_coex_contours_combined,
  fig_path(file_prefix, "cont_allcov"),
  width = 9,
  height = 4,
  dpi = save_png_dpi_contours
)


############################
## 20. Coexistence contour diagrams: Sun and Shade separated
############################

contour_plot_dat_sep <- contour_plot_dat %>%
  dplyr::mutate(
    Cov_plot = factor(as.character(Cov), levels = cover_plot_levels, labels = c("Sun", "Shade")),
    fac_label = factor(as.character(fac_label), levels = c("- Facilitator", "+ Facilitator"))
  ) %>%
  dplyr::group_by(Site, fac_label, Cov_plot) %>%
  dplyr::group_modify(~ {
    n_sample <- min(nrow(.x), max_draws_per_site_fac)
    dplyr::slice_sample(.x, n = n_sample)
  }) %>%
  dplyr::ungroup()

contour_group_check_sep <- contour_plot_dat_sep %>%
  dplyr::count(Site, fac_label, Cov_plot, name = "n_draws")

print(contour_group_check_sep, n = Inf)

if (any(contour_group_check_sep$n_draws < 50)) {
  warning(
    "Some Site x facilitator x cover groups have fewer than 50 posterior draws after filtering. ",
    "Contours may be unstable or missing."
  )
}

density_plot_dat_sep <- build_contour_density_data(
  contour_plot_dat_sep,
  group_vars = c("Site", "fac_label", "Cov_plot"),
  kde_n = kde_n
) %>%
  dplyr::mutate(
    fac_label = factor(as.character(fac_label), levels = c("- Facilitator", "+ Facilitator")),
    Cov_plot = factor(as.character(Cov_plot), levels = c("Sun", "Shade"))
  )

if (nrow(density_plot_dat_sep) == 0) {
  warning(
    "density_plot_dat_sep has zero rows for the cover-separated contour plot. ",
    "The background and boundary regions will still be plotted. ",
    "Check contour_group_check_sep and the x/y limits if you expected density shading."
  )
}

contour_facets_sep <- tidyr::expand_grid(
  fac_label = factor(c("- Facilitator", "+ Facilitator"), levels = c("- Facilitator", "+ Facilitator")),
  Cov_plot = factor(c("Sun", "Shade"), levels = c("Sun", "Shade")),
  Site = factor(site_order, levels = site_order)
)

background_coexist_sep <- expand_contour_layer(background_coexist, contour_facets_sep)
boundary_pos_sep      <- expand_contour_layer(boundary_pos, contour_facets_sep)
boundary_neg_sep      <- expand_contour_layer(boundary_neg, contour_facets_sep)
region_labels_sep    <- expand_contour_layer(region_labels, contour_facets_sep)

p_coex_contours_sep <- ggplot() +
  geom_blank(
    data = contour_facets_sep,
    aes(x = 0, y = 0),
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = background_coexist_sep,
    aes(x = niche_diff, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey85",
    colour = NA
  ) +
  geom_raster(
    data = density_plot_dat_sep %>% dplyr::filter(fac_label == "- Facilitator"),
    aes(x = niche_diff, y = fitness_ratio_log10, fill = density_level),
    interpolate = FALSE
  ) +
  scale_fill_manual(
    name = "Levels",
    values = blue_level_cols,
    breaks = as.character(6:1),
    labels = as.character(6:1),
    drop = FALSE,
    guide = guide_legend(title.position = "top", order = 1, override.aes = list(colour = NA))
  ) +
  ggnewscale::new_scale_fill() +
  geom_raster(
    data = density_plot_dat_sep %>% dplyr::filter(fac_label == "+ Facilitator"),
    aes(x = niche_diff, y = fitness_ratio_log10, fill = density_level),
    interpolate = FALSE
  ) +
  scale_fill_manual(
    name = "Levels",
    values = orange_level_cols,
    breaks = as.character(6:1),
    labels = as.character(6:1),
    drop = FALSE,
    guide = guide_legend(title.position = "top", order = 2, override.aes = list(colour = NA))
  ) +
  geom_line(data = boundary_pos_sep, aes(x = niche_diff, y = lower), inherit.aes = FALSE, colour = "red", linewidth = 0.5) +
  geom_line(data = boundary_pos_sep, aes(x = niche_diff, y = upper), inherit.aes = FALSE, colour = "red", linewidth = 0.5) +
  geom_line(data = boundary_neg_sep, aes(x = niche_diff, y = lower), inherit.aes = FALSE, colour = "grey40", linewidth = 0.5) +
  geom_line(data = boundary_neg_sep, aes(x = niche_diff, y = upper), inherit.aes = FALSE, colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey25") +
  geom_text(
    data = region_labels_sep,
    aes(x = niche_diff, y = fitness_ratio_log10, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    colour = "black",
    fontface = "bold",
    lineheight = 0.9
  ) +
  facet_grid(fac_label + Cov_plot ~ Site) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = c(-1.0, -0.50, 0, 0.5, 0.8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(y_min_log, y_max_log),
    breaks = log10(c(0.1, 0.5, 1, 3, 5)),
    labels = c("0.1", "0.5", "1.0", "3.0", "5.0"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = expression("Niche difference: " ~ 1 - rho),
    y = expression("Fitness difference: " ~ kappa[j] / kappa[i])
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.35, fill = NA),
    panel.spacing.x = grid::unit(0.35, "lines"),
    panel.spacing.y = grid::unit(0.35, "lines"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 9.5, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 8.5, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 8.5, colour = "black"),
    legend.text = element_text(size = 7.5, colour = "black"),
    legend.key.height = grid::unit(0.8, "lines"),
    legend.spacing.y = grid::unit(0.15, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

print(p_coex_contours_sep)

safe_save_plot(
  p_coex_contours_sep,
  fig_path(file_prefix, "cont_cov"),
  width = 9,
  height = 7.2,
  dpi = save_png_dpi_contours
)




  list(
    family_tag = family_tag,
    file_prefix = file_prefix,
    fit_file = fit_file,
    selected_fac_levels = selected_fac_levels,
    selected_fac_levels_wide = selected_fac_levels_wide,
    outcome_bar_summary_cover = outcome_bar_summary_cover,
    outcome_bar_long_cover = outcome_bar_long_cover,
    coex_summary = coex_summary,
    ldgr_summary = ldgr_summary,
    p_coexist_ldgr = p_coexist_ldgr,
    eta_validity_summary = eta_validity_summary,
    outcome_sum_check = outcome_sum_check,
    paths = list(
      tables_dir = tables_dir,
      figures_dir = figures_dir,
      summaries_dir = summaries_dir,
      draws_dir = draws_dir
    )
  )
}


############################
## 6. Compare BH versus Ricker outcomes
############################

compare_BH_Ricker <- function(results) {

  if (!all(c("BH", "Ricker") %in% names(results))) {
    warning("BH/Ricker comparison skipped because both families were not returned.")
    return(invisible(NULL))
  }

  outcome_cols_compare <- c(
    "TRCY wins" = "#83E8BA",
    "TROR wins" = "grey",
    "Coexistence" = "#72A1E5",
    "Priority effect" = "#9883E5"
  )

  make_outcome_long <- function(x, model_name) {
    x$outcome_bar_summary_cover %>%
      dplyr::select(
        Site, Cov, fac_present,
        p_TRCY_wins, p_TROR_wins, p_coexist, p_priority
      ) %>%
      tidyr::pivot_longer(
        cols = c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority),
        names_to = "outcome",
        values_to = model_name
      ) %>%
      dplyr::mutate(
        outcome = dplyr::case_when(
          outcome == "p_TRCY_wins" ~ "TRCY wins",
          outcome == "p_TROR_wins" ~ "TROR wins",
          outcome == "p_coexist" ~ "Coexistence",
          outcome == "p_priority" ~ "Priority effect",
          TRUE ~ outcome
        )
      )
  }

  ricker_long <- make_outcome_long(results$Ricker, "Ricker")
  bh_long     <- make_outcome_long(results$BH, "BH")

  outcome_delta <- ricker_long %>%
    dplyr::inner_join(
      bh_long,
      by = c("Site", "Cov", "fac_present", "outcome")
    ) %>%
    dplyr::mutate(
      delta_Ricker_minus_BH = Ricker - BH,
      abs_delta = abs(delta_Ricker_minus_BH),
      change_gt_5pct  = abs_delta >= 0.05,
      change_gt_10pct = abs_delta >= 0.10,
      Site = factor(as.character(Site), levels = levels(results$Ricker$outcome_bar_summary_cover$Site)),
      Cov = factor(as.character(Cov), levels = c("Sun", "Shade")),
      fac_present = factor(as.character(fac_present), levels = c("No", "Yes")),
      outcome = factor(outcome, levels = names(outcome_cols_compare))
    )

  dominant_from_summary <- function(x, model_name) {
    x$outcome_bar_summary_cover %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        dominant_outcome = c("TRCY wins", "TROR wins", "Coexistence", "Priority effect")[
          which.max(c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority))
        ],
        dominant_probability = max(c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority), na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(Site, Cov, fac_present, dominant_outcome, dominant_probability) %>%
      dplyr::rename(
        !!paste0("dominant_", model_name) := dominant_outcome,
        !!paste0("dominant_p_", model_name) := dominant_probability
      )
  }

  dominant_outcome_compare <- dominant_from_summary(results$Ricker, "Ricker") %>%
    dplyr::inner_join(
      dominant_from_summary(results$BH, "BH"),
      by = c("Site", "Cov", "fac_present")
    ) %>%
    dplyr::mutate(
      same_dominant_outcome = dominant_Ricker == dominant_BH,
      delta_dominant_probability = dominant_p_Ricker - dominant_p_BH
    ) %>%
    dplyr::arrange(Site, Cov, fac_present)

  delta_summary_by_outcome <- outcome_delta %>%
    dplyr::group_by(outcome) %>%
    dplyr::summarise(
      mean_delta_Ricker_minus_BH = mean(delta_Ricker_minus_BH, na.rm = TRUE),
      mean_abs_delta = mean(abs_delta, na.rm = TRUE),
      median_abs_delta = stats::median(abs_delta, na.rm = TRUE),
      max_abs_delta = max(abs_delta, na.rm = TRUE),
      n_cells_ge_5pct = sum(change_gt_5pct, na.rm = TRUE),
      n_cells_ge_10pct = sum(change_gt_10pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(max_abs_delta))

  overall_delta_summary <- outcome_delta %>%
    dplyr::summarise(
      n_probabilities_compared = dplyr::n(),
      mean_abs_delta = mean(abs_delta, na.rm = TRUE),
      median_abs_delta = stats::median(abs_delta, na.rm = TRUE),
      max_abs_delta = max(abs_delta, na.rm = TRUE),
      n_ge_5pct = sum(change_gt_5pct, na.rm = TRUE),
      n_ge_10pct = sum(change_gt_10pct, na.rm = TRUE)
    )

  dominant_summary <- dominant_outcome_compare %>%
    dplyr::summarise(
      n_site_cover_fac_cases = dplyr::n(),
      n_same_dominant = sum(same_dominant_outcome, na.rm = TRUE),
      n_different_dominant = sum(!same_dominant_outcome, na.rm = TRUE),
      prop_same_dominant = mean(same_dominant_outcome, na.rm = TRUE)
    )

  ldgr_delta <- results$Ricker$p_coexist_ldgr %>%
    dplyr::select(Site, panel_row, p_coexist_Ricker = p_coexist) %>%
    dplyr::inner_join(
      results$BH$p_coexist_ldgr %>%
        dplyr::select(Site, panel_row, p_coexist_BH = p_coexist),
      by = c("Site", "panel_row")
    ) %>%
    dplyr::mutate(
      delta_Ricker_minus_BH = p_coexist_Ricker - p_coexist_BH,
      abs_delta = abs(delta_Ricker_minus_BH)
    ) %>%
    dplyr::arrange(desc(abs_delta))

  ldgr_delta_summary <- ldgr_delta %>%
    dplyr::summarise(
      n_ldgr_cases = dplyr::n(),
      mean_abs_delta = mean(abs_delta, na.rm = TRUE),
      median_abs_delta = stats::median(abs_delta, na.rm = TRUE),
      max_abs_delta = max(abs_delta, na.rm = TRUE)
    )

  safe_write_csv(outcome_delta,             table_path("cmp", "out_delta"))
  safe_write_csv(dominant_outcome_compare,  table_path("cmp", "dominant"))
  safe_write_csv(delta_summary_by_outcome,  table_path("cmp", "out_sum"))
  safe_write_csv(overall_delta_summary,     table_path("cmp", "overall"))
  safe_write_csv(dominant_summary,          table_path("cmp", "dom_sum"))
  safe_write_csv(ldgr_delta,                table_path("cmp", "ldgr_delta"))
  safe_write_csv(ldgr_delta_summary,        table_path("cmp", "ldgr_sum"))

  p_delta <- ggplot(
    outcome_delta,
    aes(x = fac_present, y = delta_Ricker_minus_BH, fill = outcome)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      colour = "white",
      linewidth = 0.15
    ) +
    facet_grid(Cov ~ Site) +
    scale_fill_manual(values = outcome_cols_compare, drop = FALSE) +
    scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    labs(
      x = "Facilitators present",
      y = "Ricker - BH posterior probability",
      fill = "Outcome"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.position = "right"
    )

  print(p_delta)

  safe_save_plot(
    p_delta,
    fig_path("cmp", "out_delta"),
    width = 9,
    height = 5,
    dpi = save_png_dpi_main
  )

  p_coex_scatter_dat <- outcome_delta %>%
    dplyr::filter(outcome == "Coexistence")

  p_coex_scatter <- ggplot(
    p_coex_scatter_dat,
    aes(x = BH, y = Ricker, colour = fac_present, shape = Cov)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(size = 2.8, alpha = 0.9) +
    facet_wrap(~ Site, nrow = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(
      x = "BH posterior coexistence probability",
      y = "Ricker posterior coexistence probability",
      colour = "Facilitators present",
      shape = "Cover"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.position = "right"
    )

  print(p_coex_scatter)

  safe_save_plot(
    p_coex_scatter,
    fig_path("cmp", "coex_scatter"),
    width = 9,
    height = 3.5,
    dpi = save_png_dpi_main
  )

  p_ldgr_delta <- ggplot(
    ldgr_delta,
    aes(x = panel_row, y = delta_Ricker_minus_BH)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_col(width = 0.75) +
    facet_wrap(~ Site, nrow = 1) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(
      x = "Cover-facilitator condition",
      y = "Ricker - BH LDGR coexistence probability"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.title = element_text(colour = "black")
    )

  print(p_ldgr_delta)

  safe_save_plot(
    p_ldgr_delta,
    fig_path("cmp", "ldgr_delta"),
    width = 9,
    height = 3.5,
    dpi = save_png_dpi_main
  )

  capture.output(
    {
      cat("BH versus Ricker coexistence-outcome comparison\n")
      cat("================================================\n\n")

      cat("Overall posterior outcome probability differences, Ricker - BH:\n")
      print(overall_delta_summary)

      cat("\nDominant outcome agreement:\n")
      print(dominant_summary)

      cat("\nOutcome-specific differences:\n")
      print(delta_summary_by_outcome, n = Inf)

      cat("\nLDGR coexistence probability differences:\n")
      print(ldgr_delta_summary)

      cat("\nDominant outcome comparison by Site x Cover x Facilitator:\n")
      print(dominant_outcome_compare, n = Inf)

      cat("\nRows with posterior outcome probability difference >= 10 percentage points:\n")
      print(
        outcome_delta %>%
          dplyr::filter(change_gt_10pct) %>%
          dplyr::arrange(desc(abs_delta)),
        n = Inf
      )

      cat("\nLargest LDGR coexistence probability differences:\n")
      print(ldgr_delta, n = Inf)
    },
    file = summary_path("cmp", "summary", ext = "txt")
  )

  cat("\nBH versus Ricker comparison summary:\n")
  print(overall_delta_summary)
  print(dominant_summary)
  print(delta_summary_by_outcome, n = Inf)
  print(ldgr_delta_summary)

  invisible(
    list(
      outcome_delta = outcome_delta,
      dominant_outcome_compare = dominant_outcome_compare,
      delta_summary_by_outcome = delta_summary_by_outcome,
      overall_delta_summary = overall_delta_summary,
      dominant_summary = dominant_summary,
      ldgr_delta = ldgr_delta,
      ldgr_delta_summary = ldgr_delta_summary
    )
  )
}


############################
## 7. Run models and comparison
############################

results <- lapply(families_to_run, run_one_family)
names(results) <- vapply(results, function(x) x$family_tag, character(1))

comparison_results <- compare_BH_Ricker(results)


############################
## 8. Done
############################

cat("\nDone.\n", sep = "")
cat("Figures written to:   ", safe_path_report(figures_dir), "\n", sep = "")
cat("Tables written to:    ", safe_path_report(tables_dir), "\n", sep = "")
cat("Summaries written to: ", safe_path_report(summaries_dir), "\n", sep = "")
cat("Draws written to:     ", safe_path_report(draws_dir), "\n", sep = "")


############################
## 9. New figure that links the increase in mean seed production with the coexistence outcome
## The goal is to visually link the interpretation/main results that faciltators do increase fecundity
## BUT that does not translate into increased coexistence
############################
############################

## Uses Bayesian coexistence-model outputs only.
##
## x-axis:
##   facilitation asymmetry =
##   [log(lambda_TRCY_high) - log(lambda_TRCY_none)] -
##   [log(lambda_TROR_high) - log(lambda_TROR_none)]
##
## y-axis:
##   posterior probability of coexistence under high facilitator
##
## Interpretation:
##   x > 0 means facilitators increase TRCY lambda more than TROR lambda.
##   x < 0 means facilitators increase TROR lambda more than TRCY lambda.
##   x = 0 means symmetric facilitation effects.
############################


## Section 9 needs these objects in the global environment.
## In the main script they were created inside run_one_family(), so they disappear
## after each model finishes.

infer_factor_levels <- function(x, preferred = NULL) {
  vals <- unique(as.character(x))
  lv <- levels(x)
  
  out <- if (!is.null(lv) && length(lv) > 0) {
    lv[lv %in% vals]
  } else {
    vals
  }
  
  if (!is.null(preferred)) {
    out <- c(preferred[preferred %in% out], setdiff(out, preferred))
  }
  
  out
}

site_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Site,
  preferred = c("BEN", "GH", "NAM", "PJ", "CS")
)

cov_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Cov,
  preferred = c("Shade", "Sun")
)

cover_plot_levels <- c("Sun", "Shade")
cover_plot_levels <- cover_plot_levels[cover_plot_levels %in% cov_order]

if (length(cover_plot_levels) == 0) {
  cover_plot_levels <- cov_order
}

fac_join_digits <- 10

outcome_cols <- c(
  "TRCY wins" = "#83E8BA",
  "TROR wins" = "grey",
  "Coexistence" = "#72A1E5",
  "Priority effect" = "#9883E5"
)


make_fecundity_asymmetry_plot_data <- function(one_result) {
  
  model_name <- one_result$family_tag
  file_prefix_i <- one_result$file_prefix
  
  lambda_file <- draw_path(file_prefix_i, "lambda")
  
  if (!file.exists(lambda_file)) {
    stop("Cannot find lambda draw file for ", model_name, ": ", lambda_file)
  }
  
  lambda_i <- readRDS(lambda_file) %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      focsp = as.character(focsp),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    )
  
  selected_fac_i <- one_result$selected_fac_levels %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_level = as.character(fac_level),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    ) %>%
    dplyr::filter(fac_level %in% c("None", "High")) %>%
    dplyr::select(Site, Cov, fac_level, fac_sc_join) %>%
    dplyr::distinct()
  
  lambda_selected <- lambda_i %>%
    dplyr::inner_join(
      selected_fac_i,
      by = c("Site", "Cov", "fac_sc_join")
    ) %>%
    dplyr::select(draw, Site, Cov, focsp, fac_level, lambda)
  
  if (nrow(lambda_selected) == 0) {
    stop(
      "No lambda draws matched selected facilitator levels for ", model_name,
      ". Check fac_sc rounding/joining."
    )
  }
  
  asym_draws <- lambda_selected %>%
    dplyr::mutate(
      focsp = as.character(focsp),
      fac_level = as.character(fac_level),
      log_lambda = log(lambda)
    ) %>%
    dplyr::select(draw, Site, Cov, focsp, fac_level, log_lambda) %>%
    tidyr::pivot_wider(
      names_from = c(focsp, fac_level),
      values_from = log_lambda,
      names_glue = "loglam_{focsp}_{fac_level}"
    )
  
  required_asym_cols <- c(
    "loglam_TRCY_None", "loglam_TRCY_High",
    "loglam_TROR_None", "loglam_TROR_High"
  )
  
  missing_asym_cols <- setdiff(required_asym_cols, names(asym_draws))
  
  if (length(missing_asym_cols) > 0) {
    stop(
      "Missing columns needed for fecundity asymmetry in ", model_name, ": ",
      paste(missing_asym_cols, collapse = ", ")
    )
  }
  
  asym_summary <- asym_draws %>%
    dplyr::mutate(
      dlog_lambda_TRCY = loglam_TRCY_High - loglam_TRCY_None,
      dlog_lambda_TROR = loglam_TROR_High - loglam_TROR_None,
      fac_asym_TRCY_minus_TROR = dlog_lambda_TRCY - dlog_lambda_TROR
    ) %>%
    dplyr::group_by(Site, Cov) %>%
    dplyr::summarise(
      fac_asym_med = stats::median(fac_asym_TRCY_minus_TROR, na.rm = TRUE),
      fac_asym_lwr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.025, na.rm = TRUE),
      fac_asym_upr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.975, na.rm = TRUE),
      p_asym_gt0 = mean(fac_asym_TRCY_minus_TROR > 0, na.rm = TRUE),
      dlog_TRCY_med = stats::median(dlog_lambda_TRCY, na.rm = TRUE),
      dlog_TROR_med = stats::median(dlog_lambda_TROR, na.rm = TRUE),
      .groups = "drop"
    )
  
  outcome_high <- one_result$outcome_bar_summary_cover %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_present = as.character(fac_present)
    ) %>%
    dplyr::filter(fac_present == "Yes") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dominant_outcome = c(
        "TRCY wins",
        "TROR wins",
        "Coexistence",
        "Priority effect"
      )[
        which.max(c(
          p_TRCY_wins,
          p_TROR_wins,
          p_coexist,
          p_priority
        ))
      ],
      dominant_probability = max(
        c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority),
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Site, Cov,
      p_coexist,
      p_TRCY_wins,
      p_TROR_wins,
      p_priority,
      dominant_outcome,
      dominant_probability,
      p_valid_kappa
    )
  
  asym_summary %>%
    dplyr::left_join(outcome_high, by = c("Site", "Cov")) %>%
    dplyr::mutate(
      model_family = model_name,
      model_family = factor(model_family, levels = c("BH", "Ricker")),
      Cov = factor(as.character(Cov), levels = cover_plot_levels),
      dominant_outcome = factor(
        dominant_outcome,
        levels = c("TRCY wins", "TROR wins", "Coexistence", "Priority effect")
      ),
      point_label = paste(Site, Cov, sep = "-")
    )
}

fec_asym_plot_dat <- purrr::map_dfr(
  results,
  make_fecundity_asymmetry_plot_data
)

safe_write_csv(fec_asym_plot_dat, table_path("link", "fec_asym_pcoex"))

fec_asym_plot_dat <- fec_asym_plot_dat %>%
  dplyr::mutate(
    dominant_outcome = forcats::fct_drop(dominant_outcome)
  )

p_fec_asym_pcoex <- ggplot(
  fec_asym_plot_dat,
  aes(
    x = fac_asym_med,
    y = p_coexist,
    colour = dominant_outcome,
    shape = Cov
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    linewidth = 0.35,
    colour = "grey55"
  ) +
  geom_errorbar(
    aes(xmin = fac_asym_lwr, xmax = fac_asym_upr),
    orientation = "y",
    width = 0,
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    aes(label = Site),
    size = 3.1,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    show.legend = FALSE
  ) +
  facet_wrap(~ model_family, nrow = 1) +
  scale_colour_manual(
    values = outcome_cols,
    drop = TRUE,
    name = "Most-supported\noutcome"
  ) +
  scale_shape_manual(
    values = c("Sun" = 16, "Shade" = 17),
    drop = FALSE,
    name = "Cover"
  ) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = expression(
      "Facilitation effect asymmetry: " *
        Delta * " log(lambda)[" * TRCY * "] - " *
        Delta * " log(lambda)[" * TROR * "]"
    ),
    y = "Posterior probability of coexistence",
    caption = "Positive x-values indicate stronger facilitator effects on TRCY; negative values indicate stronger effects on TROR. Error bars show 95% posterior credible intervals for facilitation asymmetry."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10, colour = "black"),
    legend.text = element_text(size = 9, colour = "black"),
    plot.caption = element_text(size = 8.5, colour = "grey25", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

print(p_fec_asym_pcoex)

safe_save_plot(
  p_fec_asym_pcoex,
  fig_path("link", "fec_asym_pcoex"),
  width = 8.5,
  height = 4.8,
  dpi = save_png_dpi_main
)


###############################
## 10. try new version
##############################
infer_factor_levels <- function(x, preferred = NULL) {
  vals <- unique(as.character(x))
  lv <- levels(x)
  
  out <- if (!is.null(lv) && length(lv) > 0) {
    lv[lv %in% vals]
  } else {
    vals
  }
  
  if (!is.null(preferred)) {
    out <- c(preferred[preferred %in% out], setdiff(out, preferred))
  }
  
  out
}

site_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Site,
  preferred = c("BEN", "GH", "NAM", "PJ", "CS")
)

cov_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Cov,
  preferred = c("Shade", "Sun")
)

cover_plot_levels <- c("Sun", "Shade")
cover_plot_levels <- cover_plot_levels[cover_plot_levels %in% cov_order]

if (length(cover_plot_levels) == 0) {
  cover_plot_levels <- cov_order
}

fac_join_digits <- 10

outcome_cols <- c(
  "TRCY wins" = "#83E8BA",
  "TROR wins" = "grey",
  "Coexistence" = "#72A1E5",
  "Priority effect" = "#9883E5"
)


fec_asym_plot_dat <- fec_asym_plot_dat %>%
  dplyr::mutate(
    win_balance_TRCY_minus_TROR = p_TRCY_wins - p_TROR_wins,
    dominant_outcome = forcats::fct_drop(dominant_outcome)
  )


make_fecundity_asymmetry_plot_data <- function(one_result) {
  
  model_name <- one_result$family_tag
  file_prefix_i <- one_result$file_prefix
  
  lambda_file <- draw_path(file_prefix_i, "lambda")
  
  if (!file.exists(lambda_file)) {
    stop("Cannot find lambda draw file for ", model_name, ": ", lambda_file)
  }
  
  lambda_i <- readRDS(lambda_file) %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      focsp = as.character(focsp),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    )
  
  selected_fac_i <- one_result$selected_fac_levels %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_level = as.character(fac_level),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    ) %>%
    dplyr::filter(fac_level %in% c("None", "High")) %>%
    dplyr::select(Site, Cov, fac_level, fac_sc_join) %>%
    dplyr::distinct()
  
  lambda_selected <- lambda_i %>%
    dplyr::inner_join(
      selected_fac_i,
      by = c("Site", "Cov", "fac_sc_join")
    ) %>%
    dplyr::select(draw, Site, Cov, focsp, fac_level, lambda)
  
  if (nrow(lambda_selected) == 0) {
    stop(
      "No lambda draws matched selected facilitator levels for ", model_name,
      ". Check fac_sc rounding/joining."
    )
  }
  
  asym_draws <- lambda_selected %>%
    dplyr::mutate(
      focsp = as.character(focsp),
      fac_level = as.character(fac_level),
      log_lambda = log(lambda)
    ) %>%
    dplyr::select(draw, Site, Cov, focsp, fac_level, log_lambda) %>%
    tidyr::pivot_wider(
      names_from = c(focsp, fac_level),
      values_from = log_lambda,
      names_glue = "loglam_{focsp}_{fac_level}"
    )
  
  required_asym_cols <- c(
    "loglam_TRCY_None", "loglam_TRCY_High",
    "loglam_TROR_None", "loglam_TROR_High"
  )
  
  missing_asym_cols <- setdiff(required_asym_cols, names(asym_draws))
  
  if (length(missing_asym_cols) > 0) {
    stop(
      "Missing columns needed for fecundity asymmetry in ", model_name, ": ",
      paste(missing_asym_cols, collapse = ", ")
    )
  }
  
  asym_summary <- asym_draws %>%
    dplyr::mutate(
      dlog_lambda_TRCY = loglam_TRCY_High - loglam_TRCY_None,
      dlog_lambda_TROR = loglam_TROR_High - loglam_TROR_None,
      fac_asym_TRCY_minus_TROR = dlog_lambda_TRCY - dlog_lambda_TROR
    ) %>%
    dplyr::group_by(Site, Cov) %>%
    dplyr::summarise(
      fac_asym_med = stats::median(fac_asym_TRCY_minus_TROR, na.rm = TRUE),
      fac_asym_lwr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.025, na.rm = TRUE),
      fac_asym_upr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.975, na.rm = TRUE),
      p_asym_gt0 = mean(fac_asym_TRCY_minus_TROR > 0, na.rm = TRUE),
      dlog_TRCY_med = stats::median(dlog_lambda_TRCY, na.rm = TRUE),
      dlog_TROR_med = stats::median(dlog_lambda_TROR, na.rm = TRUE),
      .groups = "drop"
    )
  
  outcome_high <- one_result$outcome_bar_summary_cover %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_present = as.character(fac_present)
    ) %>%
    dplyr::filter(fac_present == "Yes") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dominant_outcome = c(
        "TRCY wins",
        "TROR wins",
        "Coexistence",
        "Priority effect"
      )[
        which.max(c(
          p_TRCY_wins,
          p_TROR_wins,
          p_coexist,
          p_priority
        ))
      ],
      dominant_probability = max(
        c(p_TRCY_wins, p_TROR_wins, p_coexist, p_priority),
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Site, Cov,
      p_coexist,
      p_TRCY_wins,
      p_TROR_wins,
      p_priority,
      dominant_outcome,
      dominant_probability,
      p_valid_kappa
    )
  
  asym_summary %>%
    dplyr::left_join(outcome_high, by = c("Site", "Cov")) %>%
    dplyr::mutate(
      model_family = model_name,
      model_family = factor(model_family, levels = c("BH", "Ricker")),
      Cov = factor(as.character(Cov), levels = cover_plot_levels),
      dominant_outcome = factor(
        dominant_outcome,
        levels = c("TRCY wins", "TROR wins", "Coexistence", "Priority effect")
      ),
      point_label = paste(Site, Cov, sep = "-")
    )
}




fec_asym_plot_dat <- purrr::map_dfr(
  results,
  make_fecundity_asymmetry_plot_data
) %>%
  dplyr::mutate(
    win_balance_TRCY_minus_TROR = p_TRCY_wins - p_TROR_wins,
    dominant_outcome = forcats::fct_drop(dominant_outcome)
  )


fec_asym_plot_dat <- fec_asym_plot_dat %>%
  dplyr::mutate(
    dominant_outcome = forcats::fct_drop(dominant_outcome)
  )


safe_write_csv(fec_asym_plot_dat, table_path("link", "fec_asym_winbalance"))


p_fec_asym_winbalance <- ggplot(
  fec_asym_plot_dat,
  aes(
    x = fac_asym_med,
    y = win_balance_TRCY_minus_TROR,
    colour = dominant_outcome,
    shape = Cov
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey35"
  ) +
  geom_errorbar(
    aes(xmin = fac_asym_lwr, xmax = fac_asym_upr),
    orientation = "y",
    width = 0,
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    aes(label = Site),
    size = 3.1,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    show.legend = FALSE
  ) +
  facet_wrap(~ model_family, nrow = 1) +
  scale_colour_manual(
    values = outcome_cols,
    drop = TRUE,
    name = "Most-supported\noutcome"
  ) +
  scale_shape_manual(
    values = c("Sun" = 16, "Shade" = 17),
    drop = FALSE,
    name = "Cover"
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    labels = function(x) paste0(round(100 * x), "%")
  ) +
  labs(
    x = expression(
      "Facilitation effect asymmetry: " *
        Delta * " log(lambda)[" * TRCY * "] - " *
        Delta * " log(lambda)[" * TROR * "]"
    ),
    y = "Posterior winner balance: P(TRCY wins) - P(TROR wins)",
    caption = "Positive x-values indicate stronger facilitator effects on TRCY; negative values indicate stronger effects on TROR.
    Positive y-values indicate stronger posterior support for TRCY winning; negative values indicate stronger support for TROR winning."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.3),
    strip.text = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10, colour = "black"),
    legend.text = element_text(size = 9, colour = "black"),
    plot.caption = element_text(size = 8.5, colour = "grey25", hjust = 0),
    plot.margin = margin(6, 5, 5, 5)
  )

print(p_fec_asym_winbalance)

safe_save_plot(
  p_fec_asym_winbalance,
  fig_path("link", "fec_asym_winbalance"),
  width = 8.5,
  height = 4.8,
  dpi = save_png_dpi_main
)

############################
## 11. Draw-level link between facilitation asymmetry and coexistence outcome
##
## This section links each posterior draw to the coexistence outcome calculated
## from that same posterior draw.
##
## x-axis:
##   [log(lambda_TRCY_high) - log(lambda_TRCY_none)] -
##   [log(lambda_TROR_high) - log(lambda_TROR_none)]
##
## Two figures are saved:
##
## 1. draw_fec_asym_outcome
##    x-axis = facilitation-effect asymmetry
##    y-axis = draw-level outcome category
##
## 2. draw_fec_asym_logfit
##    x-axis = facilitation-effect asymmetry
##    y-axis = log(fitness ratio TRCY/TROR)
##    This second figure only includes draws with positive fitness ratios,
##    because log(fitness ratio) is undefined for zero or negative values.
############################

## Recreate global plotting/order objects. These were created inside
## run_one_family(), so they are not available after the function finishes.
infer_factor_levels <- function(x, preferred = NULL) {
  vals <- unique(as.character(x))
  lv <- levels(x)
  
  out <- if (!is.null(lv) && length(lv) > 0) {
    lv[lv %in% vals]
  } else {
    vals
  }
  
  if (!is.null(preferred)) {
    out <- c(preferred[preferred %in% out], setdiff(out, preferred))
  }
  
  out
}

site_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Site,
  preferred = c("BEN", "GH", "NAM", "PJ", "CS")
)

cov_order <- infer_factor_levels(
  results[[1]]$selected_fac_levels$Cov,
  preferred = c("Shade", "Sun")
)

cover_plot_levels <- c("Sun", "Shade")
cover_plot_levels <- cover_plot_levels[cover_plot_levels %in% cov_order]

if (length(cover_plot_levels) == 0) {
  cover_plot_levels <- cov_order
}

fac_join_digits <- 10

outcome_cols_draw <- c(
  "TRCY wins" = "#0d7d87",
  "TROR wins" = "#4a2377",
  "Coexistence" = "#f55f74",
  "Priority effect" = "#8cc5e3",
  "Unclassified" = "black"
)

make_draw_level_link_data <- function(one_result) {
  
  model_name <- one_result$family_tag
  file_prefix_i <- one_result$file_prefix
  
  coex_all_file <- draw_path(file_prefix_i, "coex_all")
  
  if (!file.exists(coex_all_file)) {
    stop(
      "Cannot find coex_all draw file for ", model_name, ":\n",
      coex_all_file, "\n\n",
      "Add this inside run_one_family(), after coex_draws has outcome columns:\n",
      "safe_save_rds(coex_draws, draw_path(file_prefix, 'coex_all'))\n\n",
      "Then rerun the full outcome script before running Section 11."
    )
  }
  
  coex_all_i <- readRDS(coex_all_file) %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    )
  
  selected_fac_i <- one_result$selected_fac_levels %>%
    dplyr::mutate(
      Site = factor(as.character(Site), levels = site_order),
      Cov  = factor(as.character(Cov), levels = cov_order),
      fac_level = as.character(fac_level),
      fac_sc_join = round(as.numeric(fac_sc), fac_join_digits)
    ) %>%
    dplyr::filter(fac_level %in% c("None", "High")) %>%
    dplyr::select(Site, Cov, fac_level, fac_sc_join) %>%
    dplyr::distinct()
  
  coex_selected <- coex_all_i %>%
    dplyr::inner_join(
      selected_fac_i,
      by = c("Site", "Cov", "fac_sc_join")
    )
  
  if (nrow(coex_selected) == 0) {
    stop(
      "No coex_all draws matched selected facilitator levels for ", model_name,
      ". Check fac_sc rounding and whether selected_fac_levels$fac_sc values ",
      "were included in fac_values_all before lambda_grid was built."
    )
  }
  
  asym_draws <- coex_selected %>%
    dplyr::select(
      draw, Site, Cov, fac_level,
      lambda_TRCY, lambda_TROR
    ) %>%
    dplyr::mutate(
      log_lambda_TRCY = log(lambda_TRCY),
      log_lambda_TROR = log(lambda_TROR)
    ) %>%
    dplyr::select(
      draw, Site, Cov, fac_level,
      log_lambda_TRCY, log_lambda_TROR
    ) %>%
    tidyr::pivot_wider(
      names_from = fac_level,
      values_from = c(log_lambda_TRCY, log_lambda_TROR),
      names_glue = "{.value}_{fac_level}"
    )
  
  required_asym_cols <- c(
    "log_lambda_TRCY_None",
    "log_lambda_TRCY_High",
    "log_lambda_TROR_None",
    "log_lambda_TROR_High"
  )
  
  missing_asym_cols <- setdiff(required_asym_cols, names(asym_draws))
  
  if (length(missing_asym_cols) > 0) {
    stop(
      "Missing columns needed for draw-level facilitation asymmetry in ",
      model_name, ": ",
      paste(missing_asym_cols, collapse = ", ")
    )
  }
  
  asym_draws <- asym_draws %>%
    dplyr::mutate(
      dlog_lambda_TRCY = log_lambda_TRCY_High - log_lambda_TRCY_None,
      dlog_lambda_TROR = log_lambda_TROR_High - log_lambda_TROR_None,
      fac_asym_TRCY_minus_TROR = dlog_lambda_TRCY - dlog_lambda_TROR
    ) %>%
    dplyr::select(
      draw, Site, Cov,
      dlog_lambda_TRCY,
      dlog_lambda_TROR,
      fac_asym_TRCY_minus_TROR
    )
  
  outcome_high_draws <- coex_selected %>%
    dplyr::filter(fac_level == "High") %>%
    dplyr::mutate(
      fitness_ratio = fitness_ratio_TRCY_over_TROR,
      
      valid_draw =
        is.finite(rho) &
        rho > 0 &
        is.finite(fitness_ratio),
      
      inv_rho_draw = 1 / rho,
      lower_bound_draw = pmin(rho, inv_rho_draw),
      upper_bound_draw = pmax(rho, inv_rho_draw),
      
      outcome_draw = dplyr::case_when(
        valid_draw &
          rho < 1 &
          fitness_ratio > lower_bound_draw &
          fitness_ratio < upper_bound_draw ~ "Coexistence",
        
        valid_draw &
          rho > 1 &
          fitness_ratio > lower_bound_draw &
          fitness_ratio < upper_bound_draw ~ "Priority effect",
        
        valid_draw &
          fitness_ratio > upper_bound_draw ~ "TRCY wins",
        
        valid_draw &
          fitness_ratio < lower_bound_draw ~ "TROR wins",
        
        TRUE ~ "Unclassified"
      ),
      
      log_fitness_ratio_TRCY_over_TROR = {
        out <- rep(NA_real_, dplyr::n())
        ok <- valid_draw & fitness_ratio > 0
        out[ok] <- log(fitness_ratio[ok])
        out
      },
      
      fitness_ratio_positive = valid_draw & fitness_ratio > 0,
      fitness_ratio_zero_or_negative = valid_draw & fitness_ratio <= 0
    ) %>%
    dplyr::select(
      draw, Site, Cov,
      fitness_ratio_TRCY_over_TROR,
      log_fitness_ratio_TRCY_over_TROR,
      fitness_ratio_positive,
      fitness_ratio_zero_or_negative,
      niche_diff,
      rho,
      valid_draw,
      lower_bound_draw,
      upper_bound_draw,
      outcome_draw
    )
  
  asym_draws %>%
    dplyr::inner_join(
      outcome_high_draws,
      by = c("draw", "Site", "Cov")
    ) %>%
    dplyr::mutate(
      model_family = model_name,
      model_family = factor(model_family, levels = c("BH", "Ricker")),
      Site = factor(as.character(Site), levels = site_order),
      Cov = factor(as.character(Cov), levels = cover_plot_levels),
      outcome_draw = factor(
        outcome_draw,
        levels = c(
          "TRCY wins",
          "TROR wins",
          "Coexistence",
          "Priority effect",
          "Unclassified"
        )
      ),
      
      asym_direction = dplyr::case_when(
        fac_asym_TRCY_minus_TROR > 0 ~ "TRCY",
        fac_asym_TRCY_minus_TROR < 0 ~ "TROR",
        TRUE ~ NA_character_
      ),
      
      winner_direction = dplyr::case_when(
        outcome_draw == "TRCY wins" ~ "TRCY",
        outcome_draw == "TROR wins" ~ "TROR",
        TRUE ~ NA_character_
      ),
      
      asym_matches_winner = dplyr::case_when(
        !is.na(asym_direction) &
          !is.na(winner_direction) &
          asym_direction == winner_direction ~ TRUE,
        
        !is.na(asym_direction) &
          !is.na(winner_direction) &
          asym_direction != winner_direction ~ FALSE,
        
        TRUE ~ NA
      ),
      
      logfit_same_direction = dplyr::case_when(
        fac_asym_TRCY_minus_TROR > 0 &
          log_fitness_ratio_TRCY_over_TROR > 0 ~ TRUE,
        
        fac_asym_TRCY_minus_TROR < 0 &
          log_fitness_ratio_TRCY_over_TROR < 0 ~ TRUE,
        
        fac_asym_TRCY_minus_TROR == 0 |
          log_fitness_ratio_TRCY_over_TROR == 0 ~ NA,
        
        is.na(log_fitness_ratio_TRCY_over_TROR) ~ NA,
        
        TRUE ~ FALSE
      )
    )
}

draw_link_dat_all <- purrr::map_dfr(
  results,
  make_draw_level_link_data
) %>%
  dplyr::filter(
    is.finite(fac_asym_TRCY_minus_TROR),
    !is.na(Cov)
  )

if (nrow(draw_link_dat_all) == 0) {
  stop("draw_link_dat_all has zero finite rows. Check coex_all and lambda values.")
}

draw_link_outcome_check <- draw_link_dat_all %>%
  dplyr::count(model_family, outcome_draw, name = "n_draws") %>%
  dplyr::group_by(model_family) %>%
  dplyr::mutate(prop_draws = n_draws / sum(n_draws)) %>%
  dplyr::ungroup()

print(draw_link_outcome_check, n = Inf)

safe_write_csv(
  draw_link_outcome_check,
  table_path("link", "draw_fec_asym_outcome_check")
)

safe_save_rds(
  draw_link_dat_all,
  draw_path("link", "draw_fec_asym_all")
)

draw_link_summary <- draw_link_dat_all %>%
  dplyr::group_by(model_family, Site, Cov) %>%
  dplyr::summarise(
    n_draws = dplyr::n(),
    
    fac_asym_med = stats::median(fac_asym_TRCY_minus_TROR, na.rm = TRUE),
    fac_asym_lwr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.025, na.rm = TRUE),
    fac_asym_upr = stats::quantile(fac_asym_TRCY_minus_TROR, 0.975, na.rm = TRUE),
    
    dlog_TRCY_med = stats::median(dlog_lambda_TRCY, na.rm = TRUE),
    dlog_TROR_med = stats::median(dlog_lambda_TROR, na.rm = TRUE),
    
    p_asym_TRCY_gt_TROR = mean(fac_asym_TRCY_minus_TROR > 0, na.rm = TRUE),
    p_asym_TROR_gt_TRCY = mean(fac_asym_TRCY_minus_TROR < 0, na.rm = TRUE),
    
    p_valid_draw = mean(valid_draw, na.rm = TRUE),
    p_fitness_ratio_positive = mean(fitness_ratio_positive, na.rm = TRUE),
    p_fitness_ratio_zero_or_negative = mean(fitness_ratio_zero_or_negative, na.rm = TRUE),
    
    p_TRCY_wins = mean(outcome_draw == "TRCY wins", na.rm = TRUE),
    p_TROR_wins = mean(outcome_draw == "TROR wins", na.rm = TRUE),
    p_coexist = mean(outcome_draw == "Coexistence", na.rm = TRUE),
    p_priority = mean(outcome_draw == "Priority effect", na.rm = TRUE),
    p_unclassified = mean(outcome_draw == "Unclassified", na.rm = TRUE),
    
    p_asym_matches_winner = mean(asym_matches_winner, na.rm = TRUE),
    
    logfit_med = stats::median(log_fitness_ratio_TRCY_over_TROR, na.rm = TRUE),
    logfit_lwr = stats::quantile(log_fitness_ratio_TRCY_over_TROR, 0.025, na.rm = TRUE),
    logfit_upr = stats::quantile(log_fitness_ratio_TRCY_over_TROR, 0.975, na.rm = TRUE),
    p_logfit_TRCY_gt_TROR = mean(log_fitness_ratio_TRCY_over_TROR > 0, na.rm = TRUE),
    p_logfit_same_direction = mean(logfit_same_direction, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(model_family, Cov, Site)

safe_write_csv(
  draw_link_summary,
  table_path("link", "draw_fec_asym_sum")
)

print(draw_link_summary, n = Inf)

set.seed(123)

max_draws_per_panel_link <- 1000

sample_draws_by_panel <- function(dat, max_n = 1000) {
  dat %>%
    dplyr::group_by(model_family, Site, Cov) %>%
    dplyr::group_modify(~ {
      n_sample <- min(nrow(.x), max_n)
      dplyr::slice_sample(.x, n = n_sample)
    }) %>%
    dplyr::ungroup()
}

make_robust_limits <- function(x, probs = c(0.005, 0.995), pad_prop = 0.08) {
  lim <- stats::quantile(x, probs = probs, na.rm = TRUE)
  lim <- as.numeric(lim)
  
  lim <- c(
    min(lim[1], 0, na.rm = TRUE),
    max(lim[2], 0, na.rm = TRUE)
  )
  
  span <- diff(lim)
  
  if (!is.finite(span) || span == 0) {
    span <- 1
  }
  
  lim + c(-1, 1) * span * pad_prop
}

############################
## 11. Plotting: log-fitness-ratio only
############################

draw_link_logfit_dat <- draw_link_dat_all %>%
  dplyr::filter(
    is.finite(log_fitness_ratio_TRCY_over_TROR)
  )

if (nrow(draw_link_logfit_dat) == 0) {
  
  warning(
    "No draws have positive finite fitness ratios, so the log-fitness-ratio ",
    "plots were skipped."
  )
  
} else {
  
  draw_link_logfit_plot_dat <- draw_link_logfit_dat %>%
    sample_draws_by_panel(max_n = max_draws_per_panel_link) %>%
    dplyr::mutate(
      outcome_draw = forcats::fct_drop(outcome_draw)
    )
  
  x_lims_logfit <- make_robust_limits(
    draw_link_logfit_plot_dat$fac_asym_TRCY_minus_TROR
  )
  
  y_lims_logfit <- make_robust_limits(
    draw_link_logfit_plot_dat$log_fitness_ratio_TRCY_over_TROR
  )
  
 
  ############################
  ## 11b. Combined log-fitness-ratio plots
  ## Separate plot for each model family
  ##
  ## Sites = symbols
  ## Sun = open symbols
  ## Shade = closed symbols
  ############################
  
  site_shape_vals <- c(
    "BEN" = 21,
    "GH"  = 22,
    "NAM" = 23,
    "PJ"  = 24,
    "CS"  = 25
  )
  
  site_shape_vals <- site_shape_vals[
    names(site_shape_vals) %in% levels(draw_link_logfit_plot_dat$Site)
  ]
  
  ## Fallback in case site labels differ from BEN/GH/NAM/PJ/CS
  if (length(site_shape_vals) == 0) {
    site_shape_vals <- stats::setNames(
      c(21, 22, 23, 24, 25)[seq_along(levels(draw_link_logfit_plot_dat$Site))],
      levels(draw_link_logfit_plot_dat$Site)
    )
  }
  
  draw_link_logfit_plot_dat_combined <- draw_link_logfit_plot_dat %>%
    dplyr::mutate(
      cover_fill = dplyr::case_when(
        Cov == "Sun" ~ "white",
        Cov == "Shade" ~ as.character(outcome_cols_draw[as.character(outcome_draw)]),
        TRUE ~ "white"
      ),
      model_family_label = dplyr::case_when(
        model_family == "BH" ~ "Beverton-Holt",
        model_family == "Ricker" ~ "Ricker",
        TRUE ~ as.character(model_family)
      )
    )
  
  model_families_to_plot <- levels(draw_link_logfit_plot_dat_combined$model_family)
  model_families_to_plot <- model_families_to_plot[
    model_families_to_plot %in% unique(draw_link_logfit_plot_dat_combined$model_family)
  ]
  
  for (family_i in model_families_to_plot) {
    
    dat_i <- draw_link_logfit_plot_dat_combined %>%
      dplyr::filter(model_family == family_i)
    
    family_label_i <- unique(dat_i$model_family_label)[1]
    
    file_family_i <- dplyr::case_when(
      family_i == "BH" ~ "BH",
      family_i == "Ricker" ~ "Ricker",
      TRUE ~ gsub("[^A-Za-z0-9]+", "_", as.character(family_i))
    )
    
    p_draw_fec_asym_logfit_combined_i <- ggplot(
      dat_i,
      aes(
        x = fac_asym_TRCY_minus_TROR,
        y = log_fitness_ratio_TRCY_over_TROR,
        colour = outcome_draw,
        shape = Site
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.35,
        colour = "grey35"
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.35,
        colour = "grey35"
      ) +
      geom_point(
        aes(fill = cover_fill),
        size = 0.5,
        alpha = 0.2,
        stroke = 0.5
      ) +
      scale_colour_manual(
        values = outcome_cols_draw,
        drop = TRUE,
        name = "Draw-level\noutcome"
      ) +
      scale_shape_manual(
        values = site_shape_vals,
        drop = TRUE,
        name = "Site"
      ) +
      scale_fill_identity() +
      guides(
        colour = guide_legend(
          override.aes = list(
            size = 2,
            alpha = 1,
            shape = 16,
            stroke = 0
          )
        ),
        shape = guide_legend(
          override.aes = list(
            size = 2,
            alpha = 1,
            colour = "black",
            fill = "black",
            stroke = 0.5
          )
        )
      ) +
      coord_cartesian(
        #xlim = x_lims_logfit,
        xlim = c(-4,4),
        ylim = y_lims_logfit
      ) +
      labs(
        title = paste0(family_label_i, " model"),
        x = expression(
          "Facilitation effect asymmetry: " *
            Delta * " log(lambda)[" * TRCY * "] - " *
            Delta * " log(lambda)[" * TROR * "]"
        ),
        y = "Log fitness ratio: log(TRCY / TROR)",
        caption = "link fecundity with coexistence outcome posterior draws"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.18, colour = "grey92"),
        axis.title = element_text(size = 10.5, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        plot.title = element_text(size = 12, colour = "black", face = "bold"),
        plot.caption = element_text(size = 8.3, colour = "grey25", hjust = 0),
        plot.margin = margin(5, 5, 5, 5)
      )
    
    print(p_draw_fec_asym_logfit_combined_i)
    
    safe_save_plot(
      p_draw_fec_asym_logfit_combined_i,
      fig_path(
        "link",
        paste0("draw_fec_asym_logfit_combined_sites_cover_", file_family_i)
      ),
      width = 6,
      height = 4,
      dpi = 600
    )
  }
}

