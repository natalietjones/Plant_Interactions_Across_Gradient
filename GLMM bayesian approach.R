############################################################
## FULL MODEL COMPARISON WORKFLOW
## Context-dependent facilitation across site and cover
##
## Candidate models:
## 1) base: no Site, no Cov
## 2) cov: Cov only
## 3) site: Site only
## 4) site_cov: Site + Cov
## 5) covfac: Site + Cov + Cov:fac
## 6) sitefac: Site + Cov + Site:fac
## 7) bothfac: Site + Cov + Cov:fac + Site:fac
##
## Also:
## - germination adequacy checks
## - CS inclusion/exclusion sensitivity analysis
## - stable caching / load existing models
## - moment-matched LOO comparison
## - predicted GLMM plots from best model
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
  "tidybayes"
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
library(tidybayes)

## Backend
options(brms.backend = "rstan")
## If cmdstanr works on your machine:
## options(brms.backend = "cmdstanr")

options(mc.cores = max(1, parallel::detectCores() - 1))

############################
## 1. User settings
############################
analysis_id  <- "env_fac_model_comparison_v1"
force_refit  <- FALSE   # TRUE if you want to rerun models from scratch
run_with_CS  <- TRUE
run_no_CS    <- FALSE


############################
## 2. Read data
############################
dat0 <- read.csv("fullWAdata.csv")

############################################################
## CHECK AND MATCH LEVELS
############################################################

site_levels <- c("Ben", "GH", "Nam", "PJ", "CS")
cov_levels  <- c("Shade", "Sun")
sp_levels   <- c("TRCY", "TROR")

clean_site <- function(x) {
  x <- trimws(as.character(x))
  x_up <- toupper(x)
  
  dplyr::case_when(
    x_up == "BEN" ~ "Ben",
    x_up == "GH"  ~ "GH",
    x_up == "NAM" ~ "Nam",
    x_up == "PJ"  ~ "PJ",
    x_up == "CS"  ~ "CS",
    TRUE ~ NA_character_
  )
}

clean_cov <- function(x) {
  x <- trimws(as.character(x))
  x_low <- tolower(x)
  
  dplyr::case_when(
    x_low %in% c("shade", "sh") ~ "Shade",
    x_low %in% c("sun", "open") ~ "Sun",
    TRUE ~ NA_character_
  )
}

clean_focsp <- function(x) {
  x <- trimws(as.character(x))
  x_up <- toupper(x)
  
  dplyr::case_when(
    x_up == "TRCY" ~ "TRCY",
    x_up == "TROR" ~ "TROR",
    TRUE ~ NA_character_
  )
}

## first inspect raw values
cat("\nRAW Site values:\n")
print(sort(unique(trimws(as.character(dat0$Site)))))

cat("\nRAW Cov values:\n")
print(sort(unique(trimws(as.character(dat0$Cov)))))

cat("\nRAW focsp values:\n")
print(sort(unique(trimws(as.character(dat0$focsp)))))

## clean them
dat0 <- dat0 %>%
  mutate(
    Site  = clean_site(Site),
    Cov   = clean_cov(Cov),
    focsp = clean_focsp(focsp),
    BLK   = trimws(as.character(BLK))
  )

## inspect cleaned values before factoring
cat("\nCLEANED Site values:\n")
print(sort(unique(dat0$Site)))

cat("\nCLEANED Cov values:\n")
print(sort(unique(dat0$Cov)))

cat("\nCLEANED focsp values:\n")
print(sort(unique(dat0$focsp)))

## show anything that failed to match
cat("\nRows with unmatched Site/Cov/focsp after cleaning:\n")
print(
  dat0 %>%
    filter(is.na(Site) | is.na(Cov) | is.na(focsp)) %>%
    select(Site, Cov, focsp) %>%
    head(20)
)

## now set factor levels
dat0 <- dat0 %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    Cov   = factor(Cov, levels = cov_levels),
    focsp = factor(focsp, levels = sp_levels),
    BLK   = factor(BLK)
  )

## verify factor levels match expected levels
cat("\nEXPECTED Site levels:\n")
print(site_levels)
cat("ACTUAL Site levels:\n")
print(levels(dat0$Site))

cat("\nEXPECTED Cov levels:\n")
print(cov_levels)
cat("ACTUAL Cov levels:\n")
print(levels(dat0$Cov))

cat("\nEXPECTED focsp levels:\n")
print(sp_levels)
cat("ACTUAL focsp levels:\n")
print(levels(dat0$focsp))

## check whether any expected levels are absent from the data
cat("\nSite counts after cleaning:\n")
print(table(dat0$Site, useNA = "ifany"))

cat("\nCov counts after cleaning:\n")
print(table(dat0$Cov, useNA = "ifany"))

cat("\nSpecies counts after cleaning:\n")
print(table(dat0$focsp, useNA = "ifany"))

## strict checks: stop if levels are wrong
stopifnot(identical(levels(dat0$Site), site_levels))
stopifnot(identical(levels(dat0$Cov), cov_levels))
stopifnot(identical(levels(dat0$focsp), sp_levels))


############################
## 3. Prepare raw checks
############################
dat_check <- dat0 %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    Cov   = factor(Cov),
    BLK   = factor(BLK),
    focsp = factor(focsp),
    germ  = FOCAL_GERM %in% c("Y", "y", 1, "1", TRUE)
  )

############################
## 4. Germination adequacy checks
############################
germ_site_species <- dat_check %>%
  group_by(Site, focsp) %>%
  summarise(
    n_total   = n(),
    n_germ    = sum(germ, na.rm = TRUE),
    prop_germ = mean(germ, na.rm = TRUE),
    .groups   = "drop"
  )

germ_site_species_cov <- dat_check %>%
  group_by(Site, Cov, focsp) %>%
  summarise(
    n_total   = n(),
    n_germ    = sum(germ, na.rm = TRUE),
    prop_germ = mean(germ, na.rm = TRUE),
    .groups   = "drop"
  )

print(germ_site_species)
print(germ_site_species_cov)

############################
## 5. Build fecundity dataset
############################
dat_fec <- dat0 %>%
  mutate(
    Site     = factor(Site, levels = site_levels),
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
    intra_sc = intra,
    inter_sc = inter,
    fac_sc   = fac
  ) %>%
  droplevels()


############################
## 5B. Centre and scale facilitator density
############################
fac_center <- mean(dat_fec$fac_sc, na.rm = TRUE)
fac_scale  <- sd(dat_fec$fac_sc, na.rm = TRUE)

if (is.na(fac_scale) || fac_scale == 0) {
  stop("fac_sc has zero or undefined standard deviation; cannot fit linear/quadratic facilitation models.")
}

dat_fec <- dat_fec %>%
  mutate(
    fac_z  = (fac_sc - fac_center) / fac_scale,
    fac_z2 = fac_z^2
  )
############################
## 6. Design checks in fecundity data
############################
fec_site_species <- dat_fec %>%
  group_by(Site, focsp) %>%
  summarise(
    n_fec      = n(),
    mean_fit   = mean(fitness, na.rm = TRUE),
    mean_intra = mean(intra, na.rm = TRUE),
    mean_inter = mean(inter, na.rm = TRUE),
    mean_fac   = mean(fac, na.rm = TRUE),
    .groups    = "drop"
  )

fec_site_species_cov <- dat_fec %>%
  group_by(Site, Cov, focsp) %>%
  summarise(
    n_fec      = n(),
    mean_fit   = mean(fitness, na.rm = TRUE),
    mean_intra = mean(intra, na.rm = TRUE),
    mean_inter = mean(inter, na.rm = TRUE),
    mean_fac   = mean(fac, na.rm = TRUE),
    .groups    = "drop"
  )

print(fec_site_species)
print(fec_site_species_cov)

############################
## 7. Output folder: STABLE
############################
run_dir <- file.path(getwd(), analysis_id)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(germ_site_species,
          file.path(run_dir, "germ_site_species.csv"),
          row.names = FALSE)
write.csv(germ_site_species_cov,
          file.path(run_dir, "germ_site_species_cov.csv"),
          row.names = FALSE)
write.csv(fec_site_species,
          file.path(run_dir, "fec_site_species.csv"),
          row.names = FALSE)
write.csv(fec_site_species_cov,
          file.path(run_dir, "fec_site_species_cov.csv"),
          row.names = FALSE)
write.csv(dat_fec,
          file.path(run_dir, "fecundity_analysis_data_all_sites.csv"),
          row.names = FALSE)

############################
## 8. Exploratory plots
############################
p1 <- ggplot(dat_fec, aes(intra, fitness, colour = focsp)) +
  geom_point(alpha = 0.5) +
  facet_grid(Cov ~ Site) +
  theme_bw() +
  labs(title = "Fitness vs intraspecific density")

p2 <- ggplot(dat_fec, aes(inter, fitness, colour = focsp)) +
  geom_point(alpha = 0.5) +
  facet_grid(Cov ~ Site) +
  theme_bw() +
  labs(title = "Fitness vs interspecific density")

p3 <- ggplot(dat_fec, aes(fac, fitness, colour = focsp)) +
  geom_point(alpha = 0.5) +
  facet_grid(Cov ~ Site) +
  theme_bw() +
  labs(title = "Fitness vs facilitator density")

ggsave(file.path(run_dir, "explore_intra_all_sites.png"), p1, width = 10, height = 6, dpi = 300)
ggsave(file.path(run_dir, "explore_inter_all_sites.png"), p2, width = 10, height = 6, dpi = 300)
ggsave(file.path(run_dir, "explore_fac_all_sites.png"), p3, width = 10, height = 6, dpi = 300)

############################
## 9. Priors
############################
priors_loglin <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)

############################
## 10. Model formulas
############################
## base model: facilitator effect + densities only
form_base <- bf(
  fitness ~ 0 + focsp +
    focsp:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## add Cov only
form_cov <- bf(
  fitness ~ 0 + focsp +
    focsp:Cov +
    focsp:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## add Site only
form_site <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## add Site + Cov main effects
form_site_cov <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## does Cov modify facilitation?
form_covfac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Cov:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## does Site modify facilitation?
form_sitefac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

## do both Site and Cov modify facilitation?
form_bothfac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_sc +
    focsp:Site:fac_sc +
    focsp:Cov:fac_sc +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)


form_list <- list(
  base      = form_base,
  cov       = form_cov,
  site      = form_site,
  site_cov  = form_site_cov,
  covfac    = form_covfac,
  sitefac   = form_sitefac,
  bothfac   = form_bothfac
)

############################
## 10B. Focused context-dependent facilitation test
############################
## Compare:
## 1) no facilitation
## 2) linear facilitation
## 3) quadratic facilitation
## 4) site-varying quadratic facilitation
## 5) cover-varying quadratic facilitation
## 6) site + cover-varying quadratic facilitation
##
## fac_z is centred/scaled facilitator density
## fac_z2 is the quadratic term

form_site_cov_nofac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

form_site_cov_fac_linear <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_z +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

form_site_cov_fac_quad <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_z +
    focsp:fac_z2 +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

form_site_quadfac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_z +
    focsp:fac_z2 +
    focsp:Site:fac_z +
    focsp:Site:fac_z2 +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

form_cov_quadfac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_z +
    focsp:fac_z2 +
    focsp:Cov:fac_z +
    focsp:Cov:fac_z2 +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)

form_both_quadfac <- bf(
  fitness ~ 0 + focsp +
    focsp:Site +
    focsp:Cov +
    focsp:fac_z +
    focsp:fac_z2 +
    focsp:Site:fac_z +
    focsp:Site:fac_z2 +
    focsp:Cov:fac_z +
    focsp:Cov:fac_z2 +
    focsp:intra_sc +
    focsp:inter_sc +
    (1 | BLK)
)
############################
## 11. Fit settings
############################
n_cores <- min(4, max(1, parallel::detectCores() - 1))

fit_args <- list(
  family    = negbinomial(link = "log"),
  prior     = priors_loglin,
  chains    = 4,
  iter      = 4000,
  warmup    = 2000,
  cores     = n_cores,
  refresh   = 50,
  control   = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE)
)

############################
## 12. Helper functions
############################
fit_or_load <- function(formula, data, model_name, seed, out_dir, force_refit = FALSE) {
  
  rds_path <- file.path(out_dir, paste0(model_name, ".rds"))
  brms_file <- file.path(out_dir, model_name)
  
  if (file.exists(rds_path) && !force_refit) {
    cat("\nLoading existing model:", model_name, "\n")
    fit <- readRDS(rds_path)
  } else {
    cat("\nFitting model:", model_name, "\n")
    
    fit <- brm(
      formula   = formula,
      data      = data,
      family    = fit_args$family,
      prior     = fit_args$prior,
      chains    = fit_args$chains,
      iter      = fit_args$iter,
      warmup    = fit_args$warmup,
      cores     = fit_args$cores,
      refresh   = fit_args$refresh,
      control   = fit_args$control,
      save_pars = fit_args$save_pars,
      seed      = seed,
      file      = brms_file,
      file_refit = "on_change"
    )
    
    saveRDS(fit, rds_path)
  }
  
  capture.output(
    print(summary(fit), digits = 2),
    file = file.path(out_dir, paste0(model_name, "_summary.txt"))
  )
  
  png(file.path(out_dir, paste0(model_name, "_ppcheck.png")),
      width = 1400, height = 1000, res = 160)
  print(pp_check(fit))
  dev.off()
  
  fx <- as.data.frame(fixef(fit))
  fx$term <- rownames(fx)
  write.csv(fx, file.path(out_dir, paste0(model_name, "_fixef.csv")), row.names = FALSE)
  
  return(fit)
}

extract_key_terms <- function(fit, model_name, out_dir) {
  fx <- as.data.frame(fixef(fit))
  fx$term <- rownames(fx)
  
  key_fx <- fx %>%
    filter(
      grepl("fac_sc", term) |
        grepl("fac_z", term)  |
        grepl("Site", term)   |
        grepl("Cov", term)
    )
  
  write.csv(
    key_fx,
    file.path(out_dir, paste0(model_name, "_key_terms.csv")),
    row.names = FALSE
  )
  
  p <- ggplot(key_fx, aes(x = reorder(term, Estimate), y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0) +
    coord_flip() +
    theme_bw() +
    labs(title = paste(model_name, "key terms"), x = "", y = "Estimate")
  
  ggsave(
    file.path(out_dir, paste0(model_name, "_key_terms_plot.png")),
    p, width = 8, height = 6, dpi = 300
  )
}

compute_or_load_loo <- function(fit, loo_name, out_dir, moment_match = TRUE, force_refit = FALSE) {
  
  loo_path <- file.path(out_dir, paste0(loo_name, ".rds"))
  
  if (file.exists(loo_path) && !force_refit) {
    cat("Loading existing LOO:", loo_name, "\n")
    loo_obj <- readRDS(loo_path)
  } else {
    cat("Computing LOO:", loo_name, "\n")
    loo_obj <- loo(fit, moment_match = moment_match)
    saveRDS(loo_obj, loo_path)
  }
  
  return(loo_obj)
}

run_model_set <- function(data, label, out_dir_base, force_refit = FALSE) {
  
  out_dir <- file.path(out_dir_base, paste0("models_", label))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(data, file.path(out_dir, paste0("analysis_data_", label, ".csv")), row.names = FALSE)
  
  fits <- list()
  loos <- list()
  
  seed_base <- 1000 + nrow(data)
  
  for (i in seq_along(form_list)) {
    nm <- names(form_list)[i]
    
    fits[[nm]] <- fit_or_load(
      formula     = form_list[[nm]],
      data        = data,
      model_name  = paste0("fit_", nm, "_", label),
      seed        = seed_base + i,
      out_dir     = out_dir,
      force_refit = force_refit
    )
    
    extract_key_terms(
      fit        = fits[[nm]],
      model_name = paste0("fit_", nm, "_", label),
      out_dir    = out_dir
    )
  }
  
  cat("\nStarting LOO comparison for", label, "...\n")
  
  for (nm in names(fits)) {
    loos[[nm]] <- compute_or_load_loo(
      fit          = fits[[nm]],
      loo_name     = paste0("loo_", nm, "_", label),
      out_dir      = out_dir,
      moment_match = TRUE,
      force_refit  = force_refit
    )
  }
  
  loo_tab <- loo_compare(loos)
  print(loo_tab)
  
  capture.output(
    print(loo_tab),
    file = file.path(out_dir, paste0("loo_compare_", label, ".txt"))
  )
  
  best_model <- rownames(loo_tab)[1]
  writeLines(best_model, con = file.path(out_dir, paste0("best_model_", label, ".txt")))
  
  model_weights <- loo_model_weights(loos, method = "stacking")
  weights_df <- data.frame(
    model = names(model_weights),
    stacking_weight = as.numeric(model_weights)
  )
  
  weights_file <- file.path(out_dir, paste0("model_weights_", label, ".csv"))
  
  ok <- tryCatch({
    write.csv(weights_df, weights_file, row.names = FALSE)
    TRUE
  }, error = function(e) {
    message("Could not write model weights file: ", e$message)
    FALSE
  })
  
  return(list(
    fits       = fits,
    loos       = loos,
    loo_tab    = loo_tab,
    weights    = weights_df,
    best_model = best_model,
    out_dir    = out_dir
  ))
}
run_fac_context_test <- function(data, label, out_dir_base, force_refit = FALSE) {
  
  out_dir <- file.path(out_dir_base, paste0("fac_context_test_", label))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  fit_nofac <- fit_or_load(
    formula     = form_site_cov_nofac,
    data        = data,
    model_name  = paste0("fit_site_cov_nofac_", label),
    seed        = 7000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  fit_linear <- fit_or_load(
    formula     = form_site_cov_fac_linear,
    data        = data,
    model_name  = paste0("fit_site_cov_fac_linear_", label),
    seed        = 8000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  fit_quad <- fit_or_load(
    formula     = form_site_cov_fac_quad,
    data        = data,
    model_name  = paste0("fit_site_cov_fac_quad_", label),
    seed        = 9000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  fit_site_quad <- fit_or_load(
    formula     = form_site_quadfac,
    data        = data,
    model_name  = paste0("fit_site_quadfac_", label),
    seed        = 10000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  fit_cov_quad <- fit_or_load(
    formula     = form_cov_quadfac,
    data        = data,
    model_name  = paste0("fit_cov_quadfac_", label),
    seed        = 11000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  fit_both_quad <- fit_or_load(
    formula     = form_both_quadfac,
    data        = data,
    model_name  = paste0("fit_both_quadfac_", label),
    seed        = 12000 + nrow(data),
    out_dir     = out_dir,
    force_refit = force_refit
  )
  
  loo_nofac <- compute_or_load_loo(
    fit          = fit_nofac,
    loo_name     = paste0("loo_site_cov_nofac_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_linear <- compute_or_load_loo(
    fit          = fit_linear,
    loo_name     = paste0("loo_site_cov_fac_linear_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_quad <- compute_or_load_loo(
    fit          = fit_quad,
    loo_name     = paste0("loo_site_cov_fac_quad_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_site_quad <- compute_or_load_loo(
    fit          = fit_site_quad,
    loo_name     = paste0("loo_site_quadfac_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_cov_quad <- compute_or_load_loo(
    fit          = fit_cov_quad,
    loo_name     = paste0("loo_cov_quadfac_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_both_quad <- compute_or_load_loo(
    fit          = fit_both_quad,
    loo_name     = paste0("loo_both_quadfac_", label),
    out_dir      = out_dir,
    moment_match = TRUE,
    force_refit  = force_refit
  )
  
  loo_tab <- loo_compare(list(
    nofac     = loo_nofac,
    linear    = loo_linear,
    quad      = loo_quad,
    site_quad = loo_site_quad,
    cov_quad  = loo_cov_quad,
    both_quad = loo_both_quad
  ))
  
  print(loo_tab)
  
  capture.output(
    print(loo_tab),
    file = file.path(out_dir, paste0("loo_compare_fac_context_", label, ".txt"))
  )
  
  weights <- loo_model_weights(
    list(
      nofac     = loo_nofac,
      linear    = loo_linear,
      quad      = loo_quad,
      site_quad = loo_site_quad,
      cov_quad  = loo_cov_quad,
      both_quad = loo_both_quad
    ),
    method = "stacking"
  )
  
  weights_df <- data.frame(
    model = names(weights),
    stacking_weight = as.numeric(weights)
  )
  
  weights_file <- file.path(out_dir, paste0("model_weights_fac_context_", label, ".csv"))
  
  tryCatch({
    write.csv(weights_df, weights_file, row.names = FALSE)
  }, error = function(e) {
    message("Could not write fac-context model weights file: ", e$message)
  })
  
  fit_map <- list(
    linear    = fit_linear,
    quad      = fit_quad,
    site_quad = fit_site_quad,
    cov_quad  = fit_cov_quad,
    both_quad = fit_both_quad
  )
  
  for (nm in names(fit_map)) {
    fx <- as.data.frame(fixef(fit_map[[nm]]))
    fx$term <- rownames(fx)
    fx <- fx %>% filter(grepl("fac_z", term))
    
    write.csv(
      fx,
      file.path(out_dir, paste0(nm, "_fac_terms_", label, ".csv")),
      row.names = FALSE
    )
  }
  
  best_model <- rownames(loo_tab)[1]
  writeLines(best_model, con = file.path(out_dir, paste0("best_fac_context_model_", label, ".txt")))
  
  return(list(
    fits = list(
      nofac     = fit_nofac,
      linear    = fit_linear,
      quad      = fit_quad,
      site_quad = fit_site_quad,
      cov_quad  = fit_cov_quad,
      both_quad = fit_both_quad
    ),
    loos = list(
      nofac     = loo_nofac,
      linear    = loo_linear,
      quad      = loo_quad,
      site_quad = loo_site_quad,
      cov_quad  = loo_cov_quad,
      both_quad = loo_both_quad
    ),
    loo_tab = loo_tab,
    weights = weights_df,
    best_model = best_model,
    out_dir = out_dir
  ))
}

plot_best_fac_context_model <- function(res_obj, dat_use, run_dir, label, fac_center, fac_scale) {
  
  best_name <- res_obj$best_model
  best_fit  <- res_obj$fits[[best_name]]
  
  fac_seq <- seq(
    min(dat_use$fac_sc, na.rm = TRUE),
    max(dat_use$fac_sc, na.rm = TRUE),
    length.out = 80
  )
  
  newdat <- expand_grid(
    focsp  = levels(dat_use$focsp),
    Site   = levels(dat_use$Site),
    Cov    = levels(dat_use$Cov),
    fac_sc = fac_seq
  ) %>%
    mutate(
      fac_z    = (fac_sc - fac_center) / fac_scale,
      fac_z2   = fac_z^2,
      intra_sc = median(dat_use$intra_sc, na.rm = TRUE),
      inter_sc = median(dat_use$inter_sc, na.rm = TRUE),
      BLK      = levels(dat_use$BLK)[1]
    )
  
  pred <- fitted(
    best_fit,
    newdata = newdat,
    re_formula = NA,
    summary = TRUE
  )
  
  plot_dat <- bind_cols(newdat, as.data.frame(pred))
  
  p <- ggplot(plot_dat, aes(x = fac_sc, y = Estimate, colour = Cov, fill = Cov)) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    facet_grid(focsp ~ Site, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("Best facilitation-context model:", best_name, "-", label),
      x = "Facilitator density",
      y = "Predicted fitness"
    )
  
  print(p)
  
  ggsave(
    file.path(run_dir, paste0("best_fac_context_prediction_", label, ".png")),
    p, width = 12, height = 6, dpi = 300
  )
  
  write.csv(
    plot_dat,
    file.path(run_dir, paste0("best_fac_context_prediction_data_", label, ".csv")),
    row.names = FALSE
  )
  
  invisible(plot_dat)
}

plot_best_model <- function(res_obj, dat_use, run_dir, label) {
  
  best_name <- res_obj$best_model
  best_fit  <- res_obj$fits[[best_name]]
  
  newdat <- expand_grid(
    focsp  = levels(dat_use$focsp),
    Site   = levels(dat_use$Site),
    Cov    = levels(dat_use$Cov),
    fac_sc = seq(min(dat_use$fac_sc, na.rm = TRUE),
                 max(dat_use$fac_sc, na.rm = TRUE),
                 length.out = 60)
  ) %>%
    mutate(
      intra_sc = median(dat_use$intra_sc, na.rm = TRUE),
      inter_sc = median(dat_use$inter_sc, na.rm = TRUE),
      BLK      = levels(dat_use$BLK)[1]
    )
  
  pred <- fitted(
    best_fit,
    newdata = newdat,
    re_formula = NA,
    summary = TRUE
  )
  
  plot_dat <- bind_cols(newdat, as.data.frame(pred))
  
  p <- ggplot(plot_dat, aes(x = fac_sc, y = Estimate, colour = Cov, fill = Cov)) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    facet_grid(focsp ~ Site, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("Best model:", best_name, "-", label),
      x = "Facilitator density",
      y = "Predicted fitness"
    )
  
  print(p)
  
  ggsave(
    file.path(run_dir, paste0("best_model_prediction_", label, ".png")),
    p, width = 12, height = 6, dpi = 300
  )
  
  write.csv(plot_dat,
            file.path(run_dir, paste0("best_model_prediction_data_", label, ".csv")),
            row.names = FALSE)
  
  invisible(plot_dat)
}

extract_best_terms <- function(res_obj, run_dir, label) {
  best_fit <- res_obj$fits[[res_obj$best_model]]
  
  fx <- as.data.frame(fixef(best_fit))
  fx$term <- rownames(fx)
  
  write.csv(
    fx,
    file.path(run_dir, paste0("best_model_terms_", label, ".csv")),
    row.names = FALSE
  )
  
  fx
}




############################
## 13. Run WITH CS
############################
if (run_with_CS) {
  res_with_CS <- run_model_set(
    data        = dat_fec,
    label       = "with_CS",
    out_dir_base = run_dir,
    force_refit = force_refit
  )
  
  best_terms_with_CS <- extract_best_terms(
    res_obj = res_with_CS,
    run_dir = res_with_CS$out_dir,
    label   = "with_CS"
  )
  
  pred_with_CS <- plot_best_model(
    res_obj = res_with_CS,
    dat_use = dat_fec,
    run_dir = res_with_CS$out_dir,
    label   = "with_CS"
  )
}

############################
## 14. Run WITHOUT CS
############################
dat_no_CS <- dat_fec %>%
  filter(Site != "CS") %>%
  droplevels()

write.csv(dat_no_CS,
          file.path(run_dir, "fecundity_analysis_data_no_CS.csv"),
          row.names = FALSE)

if (run_no_CS) {
  res_no_CS <- run_model_set(
    data        = dat_no_CS,
    label       = "no_CS",
    out_dir_base = run_dir,
    force_refit = force_refit
  )
  
  best_terms_no_CS <- extract_best_terms(
    res_obj = res_no_CS,
    run_dir = res_no_CS$out_dir,
    label   = "no_CS"
  )
  
  pred_no_CS <- plot_best_model(
    res_obj = res_no_CS,
    dat_use = dat_no_CS,
    run_dir = res_no_CS$out_dir,
    label   = "no_CS"
  )
}






############################
## 14B. Focused test of context-dependent facilitation
############################
if (run_with_CS) {
  fac_context_with_CS <- run_fac_context_test(
    data         = dat_fec,
    label        = "with_CS",
    out_dir_base = run_dir,
    force_refit  = force_refit
  )
  
  pred_fac_context_with_CS <- plot_best_fac_context_model(
    res_obj     = fac_context_with_CS,
    dat_use     = dat_fec,
    run_dir     = fac_context_with_CS$out_dir,
    label       = "with_CS",
    fac_center  = fac_center,
    fac_scale   = fac_scale
  )
}

if (run_no_CS) {
  fac_context_no_CS <- run_fac_context_test(
    data         = dat_no_CS,
    label        = "no_CS",
    out_dir_base = run_dir,
    force_refit  = force_refit
  )
  
  pred_fac_context_no_CS <- plot_best_fac_context_model(
    res_obj     = fac_context_no_CS,
    dat_use     = dat_no_CS,
    run_dir     = fac_context_no_CS$out_dir,
    label       = "no_CS",
    fac_center  = fac_center,
    fac_scale   = fac_scale
  )
}
############################
## 15. Extract CS-specific terms
############################
extract_CS_terms <- function(fit, model_name, out_dir) {
  fx <- as.data.frame(fixef(fit))
  fx$term <- rownames(fx)
  
  cs_fx <- fx %>%
    filter(grepl("SiteCS", term))
  
  write.csv(
    cs_fx,
    file.path(out_dir, paste0(model_name, "_SiteCS_terms.csv")),
    row.names = FALSE
  )
  
  cs_fx
}

if (exists("res_with_CS")) {
  for (nm in names(res_with_CS$fits)) {
    extract_CS_terms(
      fit        = res_with_CS$fits[[nm]],
      model_name = paste0("fit_", nm, "_with_CS"),
      out_dir    = res_with_CS$out_dir
    )
  }
}

############################
## 16. Summaries
############################
summary_lines <- c(
  "Use these outputs to decide on the final ecological model:",
  "",
  "1. Germination adequacy:",
  "   - germ_site_species.csv",
  "   - germ_site_species_cov.csv",
  "",
  "2. Best candidate model within each dataset:",
  "   - loo_compare_with_CS.txt",
  "   - loo_compare_no_CS.txt",
  "   - model_weights_with_CS.csv",
  "   - model_weights_no_CS.csv",
  "",
  "3. Interpretation of the best model:",
  "   - best_model_terms_with_CS.csv or best_model_terms_no_CS.csv",
  "   - best_model_prediction_with_CS.png or best_model_prediction_no_CS.png",
  "",
  "4. CS sensitivity:",
  "   - compare best model type with and without CS",
  "   - inspect SiteCS terms",
  "   - inspect whether conclusions about facilitation are stable",
  "5. Facilitation-context test:",
  "   - loo_compare_fac_context_with_CS.txt (and/or no_CS)",
  "   - model_weights_fac_context_with_CS.csv",
  "   - linear_fac_terms_with_CS.csv",
  "   - quad_fac_terms_with_CS.csv",
  "   - site_quad_fac_terms_with_CS.csv",
  "   - cov_quad_fac_terms_with_CS.csv",
  "   - both_quad_fac_terms_with_CS.csv",
  "   - best_fac_context_prediction_with_CS.png",
  "   - compare nofac vs linear vs quad vs site-varying quad vs cover-varying quad vs both"
)

writeLines(summary_lines, con = file.path(run_dir, "analysis_notes.txt"))

cat("\nDone.\n")
cat("Outputs saved in:\n", normalizePath(run_dir), "\n")

# best model name
fac_context_with_CS$best_model

# fitted winning model object
best_fit <- fac_context_with_CS$fits[[fac_context_with_CS$best_model]]

# full model summary
summary(best_fit)

# fixed effects only
fixef(best_fit)

# posterior predictive check
pp_check(best_fit)


######################
# Plotting Raw fecundity data and fitted predictions from the best-supported Bayesian GLMM
######################
cov_colors <- c(
  "Sun"   = "#e9c716",
  "Shade" = "#0000a2"
)

site_levels <- c("Ben", "GH", "Nam", "PJ", "CS")

species_labs <- c(
  TRCY = "italic('T. cyanopetala')",
  TROR = "italic('T. ornata')"
)

best_fit <- fac_context_with_CS$fits[[fac_context_with_CS$best_model]]

raw_dat <- dat_fec %>%
  mutate(
    Site  = factor(Site, levels = site_levels),
    Cov   = factor(Cov, levels = c("Shade", "Sun")),
    focsp = factor(focsp, levels = c("TRCY", "TROR"))
  )

newdat <- expand_grid(
  focsp  = levels(raw_dat$focsp),
  Site   = factor(site_levels, levels = site_levels),
  Cov    = levels(raw_dat$Cov),
  fac_sc = seq(min(raw_dat$fac_sc, na.rm = TRUE),
               max(raw_dat$fac_sc, na.rm = TRUE),
               length.out = 200)
) %>%
  mutate(
    focsp    = factor(focsp, levels = c("TRCY", "TROR")),
    Cov      = factor(Cov, levels = c("Shade", "Sun")),
    fac_z    = (fac_sc - fac_center) / fac_scale,
    fac_z2   = fac_z^2,
    intra_sc = median(raw_dat$intra_sc, na.rm = TRUE),
    inter_sc = median(raw_dat$inter_sc, na.rm = TRUE),
    BLK      = levels(raw_dat$BLK)[1]
  )

pred <- fitted(
  best_fit,
  newdata = newdat,
  re_formula = NA,
  summary = TRUE
)

plot_dat <- bind_cols(newdat, as.data.frame(pred))

obs_ranges <- raw_dat %>%
  group_by(focsp, Site, Cov) %>%
  summarise(
    fac_min = min(fac_sc, na.rm = TRUE),
    fac_max = max(fac_sc, na.rm = TRUE),
    .groups = "drop"
  )

plot_dat_trim <- plot_dat %>%
  left_join(obs_ranges, by = c("focsp", "Site", "Cov")) %>%
  filter(fac_sc >= fac_min, fac_sc <= fac_max)

p_raw_fit <- ggplot() +
  geom_point(
    data = raw_dat,
    aes(x = fac_sc, y = fitness, colour = Cov),
    alpha = 0.28,
    size = 1.8,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_ribbon(
    data = plot_dat_trim,
    aes(x = fac_sc, ymin = Q2.5, ymax = Q97.5, fill = Cov),
    alpha = 0.08,
    colour = NA
  ) +
  geom_line(
    data = plot_dat_trim,
    aes(x = fac_sc, y = Estimate, colour = Cov),
    linewidth = 1.15
  ) +
  facet_grid(
    focsp ~ Site,
    scales = "free",
    labeller = labeller(focsp = as_labeller(species_labs, label_parsed))
  ) +
  scale_colour_manual(values = cov_colors) +
  scale_fill_manual(values = cov_colors) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    panel.grid.minor = element_line(linewidth = 0.2)
  ) +
  labs(
    x = "Facilitator density",
    y = "Fecundity",
    colour = "Canopy cover",
    fill   = "Canopy cover"
  )

ggsave(
  file.path(run_dir, "Figure_S2_raw_data_fitted_trends.png"),
  p_raw_fit,
  width = 12,
  height = 6,
  dpi = 300
)

print(p_raw_fit)

############################
## 17. Load and extract FINAL GLMM summary
############################

load_final_glmm <- function(run_dir,
                            label = "with_CS",
                            model_set = c("fac_context", "broad")) {
  
  model_set <- match.arg(model_set)
  
  if (model_set == "fac_context") {
    out_dir <- file.path(run_dir, paste0("fac_context_test_", label))
    best_file <- file.path(out_dir, paste0("best_fac_context_model_", label, ".txt"))
    
    fit_file_map <- c(
      nofac     = paste0("fit_site_cov_nofac_", label, ".rds"),
      linear    = paste0("fit_site_cov_fac_linear_", label, ".rds"),
      quad      = paste0("fit_site_cov_fac_quad_", label, ".rds"),
      site_quad = paste0("fit_site_quadfac_", label, ".rds"),
      cov_quad  = paste0("fit_cov_quadfac_", label, ".rds"),
      both_quad = paste0("fit_both_quadfac_", label, ".rds")
    )
    
  } else {
    out_dir <- file.path(run_dir, paste0("models_", label))
    best_file <- file.path(out_dir, paste0("best_model_", label, ".txt"))
    
    fit_file_map <- c(
      base     = paste0("fit_base_", label, ".rds"),
      cov      = paste0("fit_cov_", label, ".rds"),
      site     = paste0("fit_site_", label, ".rds"),
      site_cov = paste0("fit_site_cov_", label, ".rds"),
      covfac   = paste0("fit_covfac_", label, ".rds"),
      sitefac  = paste0("fit_sitefac_", label, ".rds"),
      bothfac  = paste0("fit_bothfac_", label, ".rds")
    )
  }
  
  if (!file.exists(best_file)) {
    stop("Cannot find best-model file: ", best_file)
  }
  
  best_name <- scan(best_file, what = "character", quiet = TRUE)
  
  if (!best_name %in% names(fit_file_map)) {
    stop("Best model name not recognised: ", best_name)
  }
  
  fit_path <- file.path(out_dir, fit_file_map[[best_name]])
  
  if (!file.exists(fit_path)) {
    stop("Cannot find saved model file: ", fit_path)
  }
  
  fit <- readRDS(fit_path)
  
  list(
    best_name = best_name,
    fit = fit,
    out_dir = out_dir,
    fit_path = fit_path,
    model_set = model_set,
    label = label
  )
}

extract_final_glmm_outputs <- function(final_obj) {
  
  fit <- final_obj$fit
  out_dir <- final_obj$out_dir
  stem <- paste0("final_", final_obj$model_set, "_", final_obj$label)
  
  # full printed summary
  summ_txt <- capture.output(print(summary(fit), digits = 3))
  writeLines(summ_txt, file.path(out_dir, paste0(stem, "_summary.txt")))
  
  # fixed effects table
  fx <- as.data.frame(fixef(fit))
  fx$term <- rownames(fx)
  rownames(fx) <- NULL
  fx <- fx %>%
    select(term, everything())
  
  write.csv(
    fx,
    file.path(out_dir, paste0(stem, "_fixef.csv")),
    row.names = FALSE
  )
  
  # key terms only
  fx_key <- fx %>%
    filter(
      grepl("fac_sc", term) |
        grepl("fac_z", term) |
        grepl("Site", term) |
        grepl("Cov", term) |
        grepl("intra_sc", term) |
        grepl("inter_sc", term)
    )
  
  write.csv(
    fx_key,
    file.path(out_dir, paste0(stem, "_key_terms.csv")),
    row.names = FALSE
  )
  
  # coefficient plot
  p_coef <- ggplot(fx_key, aes(x = reorder(term, Estimate), y = Estimate)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      x = NULL,
      y = "Posterior estimate",
      title = paste("Final", final_obj$model_set, "GLMM:", final_obj$best_name)
    )
  
  ggsave(
    file.path(out_dir, paste0(stem, "_coef_plot.png")),
    p_coef,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # posterior predictive check
  png(
    file.path(out_dir, paste0(stem, "_ppcheck.png")),
    width = 1400, height = 1000, res = 160
  )
  print(pp_check(fit))
  dev.off()
  
  # simple console print
  cat("\nLoaded final", final_obj$model_set, "model:", final_obj$best_name, "\n")
  print(fx_key)
  
  invisible(list(
    summary_text = summ_txt,
    fixef = fx,
    key_terms = fx_key,
    fit = fit
  ))
}

############################
## 18. Run extractor
############################

# final focused GLMM for interpretation
final_glmm <- load_final_glmm(
  run_dir = run_dir,
  label = "with_CS",
  model_set = "fac_context"
)

final_glmm_out <- extract_final_glmm_outputs(final_glmm)

# objects you can use directly
best_fit_final  <- final_glmm$fit
best_name_final <- final_glmm$best_name

# quick console checks
best_name_final
summary(best_fit_final)
fixef(best_fit_final)
