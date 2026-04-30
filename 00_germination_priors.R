############################################################
## 00_germination_priors.R
##
## Purpose:
## 1. Extract Trachymene cyanopetala (TRCY) and T. ornata (TROR)
##    germination data from collated_data_18_02_2020.csv
## 2. Build direct Beta posteriors for PJ and BEN
## 3. Build weakly informed Beta priors for GH, NAM, CS
## 4. Save outputs for downstream coexistence calculations
##
## Main-analysis germination definition:
## germination fraction = field_number_germinants / field_number_seeds
##
## Notes:
## - nonviable seeds are retained in the denominator
## - lost seeds are retained in the denominator
## - zero-germination rows are dropped
## - PJ  = Perenjori field rows
## - BEN = Kunjin field rows
## - Hao Ran rows excluded
############################################################

rm(list = ls())

############################
## 0. Install/load packages
############################
required_packages <- c("dplyr", "tidyr", "stringr", "tibble")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

############################
## 1. Paths and folders
############################
input_file <- "collated_data_18_02_2020.csv"
base_dir   <- "00_germination_priors"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "summaries"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "selection"), showWarnings = FALSE, recursive = TRUE)

############################
## 2. Read and clean data
############################
dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN", "NULL", "#N/A", "#VALUE!", "n/a")
)

trach <- dat %>%
  mutate(
    genus          = str_trim(str_to_lower(genus)),
    species        = str_trim(str_to_lower(species)),
    field_chambers = str_trim(str_to_lower(field_chambers)),
    location       = str_trim(location),
    researcher     = str_trim(researcher)
  ) %>%
  filter(
    genus == "trachymene",
    species %in% c("cyanopetala", "ornata")
  ) %>%
  mutate(
    sp_code = case_when(
      species == "cyanopetala" ~ "TRCY",
      species == "ornata"      ~ "TROR",
      TRUE ~ NA_character_
    )
  )

## numeric columns we may use
num_cols <- c(
  "field_number_seeds",
  "field_number_germinants",
  "field_seeds_scaled_viability",
  "field_seeds_lost"
)

num_cols <- intersect(num_cols, names(trach))
trach[num_cols] <- lapply(trach[num_cols], function(x) as.numeric(as.character(x)))

############################
## 3. Quick checks
############################
capture.output(
  {
    cat("Unique row-level locations in Trachymene subset:\n")
    print(sort(unique(trach$location)))
    cat("\nUnique field/chambers values:\n")
    print(sort(unique(trach$field_chambers)))
    cat("\nCounts by location, researcher, field_chambers, sp_code:\n")
    print(trach %>% count(location, researcher, field_chambers, sp_code, sort = TRUE))
  },
  file = file.path(base_dir, "summaries", "trach_data_overview.txt")
)

write.csv(
  trach %>% count(location, researcher, field_chambers, sp_code, sort = TRUE),
  file.path(base_dir, "summaries", "trach_counts_by_study.csv"),
  row.names = FALSE
)

############################
## 4. Keep only direct field evidence for PJ and BEN
############################
trach_direct <- trach %>%
  filter(
    field_chambers == "field",
    researcher != "Hao Ran"
  ) %>%
  mutate(
    evidence_site = case_when(
      location == "Perenjori" ~ "PJ",
      location == "Kunjin"    ~ "BEN",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(evidence_site))

write.csv(
  trach_direct %>%
    count(evidence_site, location, researcher, field_chambers, sp_code, sort = TRUE),
  file.path(base_dir, "summaries", "direct_evidence_breakdown.csv"),
  row.names = FALSE
)

############################
## 5. Build row-level germination data
############################
## Main-analysis definition:
## germination fraction = field_number_germinants / field_number_seeds
##
## Nonviable and lost seeds remain in denominator.
## Zero-germination rows are dropped.

drop_zero_germ_rows <- TRUE

germ_rows <- trach_direct %>%
  mutate(
    germ_trials  = as.numeric(as.character(field_number_seeds)),
    germ_success = as.numeric(as.character(field_number_germinants))
  ) %>%
  filter(
    !is.na(germ_trials),
    !is.na(germ_success),
    germ_trials > 0
  ) %>%
  mutate(
    germ_success = pmax(0, pmin(germ_success, germ_trials)),
    germ_fail    = germ_trials - germ_success,
    germ_frac_row = germ_success / germ_trials
  )

if (drop_zero_germ_rows) {
  germ_rows <- germ_rows %>%
    filter(germ_success > 0)
}

write.csv(
  germ_rows,
  file.path(base_dir, "data", "trach_germ_rows_used.csv"),
  row.names = FALSE
)

############################
## 6. Print row-level mean germination
##    for each species at each direct site
############################
germ_means_direct <- germ_rows %>%
  group_by(evidence_site, sp_code) %>%
  summarise(
    n_rows = dplyr::n(),
    mean_row_germination = mean(germ_frac_row, na.rm = TRUE),
    pooled_germination   = sum(germ_success, na.rm = TRUE) / sum(germ_trials, na.rm = TRUE),
    total_germ_success   = sum(germ_success, na.rm = TRUE),
    total_germ_trials    = sum(germ_trials, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nMean germination by direct site and species:\n")
print(germ_means_direct)

write.csv(
  germ_means_direct,
  file.path(base_dir, "selection", "germination_means_direct_sites.csv"),
  row.names = FALSE
)

############################
## 7. Direct posteriors for PJ and BEN
############################
## Beta posterior:
## g | data ~ Beta(success + 1, fail + 1)

germ_direct <- germ_rows %>%
  group_by(evidence_site, sp_code) %>%
  summarise(
    germ_success = sum(germ_success, na.rm = TRUE),
    germ_trials  = sum(germ_trials,  na.rm = TRUE),
    n_rows_used  = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    germ_fail = germ_trials - germ_success,
    g_alpha   = germ_success + 1,
    g_beta    = germ_fail + 1,
    g_mean    = g_alpha / (g_alpha + g_beta)
  )

cat("\nDirect Beta posterior means:\n")
print(germ_direct)

write.csv(
  germ_direct,
  file.path(base_dir, "selection", "germination_direct_posteriors.csv"),
  row.names = FALSE
)

############################
## 8. Weak priors for GH, NAM, CS
############################
## Centered on pooled PJ+BEN evidence by species
## with low effective strength.

weak_strength <- 4

germ_weak <- germ_direct %>%
  group_by(sp_code) %>%
  summarise(
    pooled_success = sum(germ_success, na.rm = TRUE),
    pooled_trials  = sum(germ_trials,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pooled_mean   = pooled_success / pooled_trials,
    g_alpha_weak  = pmax(pooled_mean * weak_strength, 0.5),
    g_beta_weak   = pmax((1 - pooled_mean) * weak_strength, 0.5)
  )

cat("\nWeak prior centres by species:\n")
print(germ_weak)

write.csv(
  germ_weak,
  file.path(base_dir, "selection", "germination_weak_priors.csv"),
  row.names = FALSE
)

############################
## 9. Final site x species germination table
############################
site_order <- c("BEN", "GH", "NAM", "PJ", "CS")

germ_site_priors <- expand.grid(
  Site  = site_order,
  focsp = c("TRCY", "TROR"),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(
    germ_direct %>%
      rename(Site = evidence_site, focsp = sp_code) %>%
      select(Site, focsp, g_alpha, g_beta, g_mean),
    by = c("Site", "focsp")
  ) %>%
  left_join(
    germ_weak %>%
      rename(focsp = sp_code) %>%
      select(focsp, g_alpha_weak, g_beta_weak, pooled_mean),
    by = "focsp"
  ) %>%
  mutate(
    g_source      = ifelse(is.na(g_alpha), "weak", "direct"),
    g_alpha_final = dplyr::coalesce(g_alpha, g_alpha_weak),
    g_beta_final  = dplyr::coalesce(g_beta,  g_beta_weak),
    g_mean_final  = dplyr::coalesce(g_mean,  pooled_mean)
  ) %>%
  select(
    Site, focsp, g_source,
    g_alpha_final, g_beta_final, g_mean_final
  )

cat("\nFinal site-level germination priors:\n")
print(germ_site_priors)

write.csv(
  germ_site_priors,
  file.path(base_dir, "selection", "germination_site_priors.csv"),
  row.names = FALSE
)

############################
## 10. Optional: print mean for each species for each final site
############################
germ_site_means <- germ_site_priors %>%
  select(Site, focsp, g_source, g_mean_final)

cat("\nMean germination used for each species at each final site:\n")
print(germ_site_means)

write.csv(
  germ_site_means,
  file.path(base_dir, "selection", "germination_site_means_final.csv"),
  row.names = FALSE
)

############################
## 11. Save final R object
############################
germination_priors <- list(
  input_file        = normalizePath(input_file),
  weak_strength     = weak_strength,
  site_order        = site_order,
  direct_rules      = list(
    PJ  = "location == 'Perenjori' & field_chambers == 'field' & researcher != 'Hao Ran'",
    BEN = "location == 'Kunjin' & field_chambers == 'field' & researcher != 'Hao Ran'"
  ),
  options = list(
    drop_zero_germ_rows = drop_zero_germ_rows,
    denominator_rule = "germ_trials = field_number_seeds",
    numerator_rule   = "germ_success = field_number_germinants",
    nonviable_rule   = "field nonviable retained in denominator",
    lost_seeds_rule  = "field_seeds_lost retained in denominator"
  ),
  row_level_data     = germ_rows,
  direct_site_means  = germ_means_direct,
  site_priors        = germ_site_priors,
  direct_posteriors  = germ_direct,
  weak_priors        = germ_weak
)

saveRDS(
  germination_priors,
  file.path(base_dir, "selection", "germination_priors.rds")
)

############################
## 12. Notes and session info
############################
note_lines <- c(
  "Germination priors created for TRCY and TROR.",
  "",
  "Interpretation:",
  "- PJ uses direct row-level Perenjori field evidence.",
  "- BEN uses direct row-level Kunjin field evidence as the Bendering analogue.",
  "- Hao Ran rows are excluded.",
  "- Zero-germination rows are dropped.",
  "- Germination fraction is field_number_germinants / field_number_seeds.",
  "- Nonviable and lost seeds are retained in the denominator.",
  "- GH, NAM, and CS use weakly informed priors centered on pooled PJ+BEN evidence."
)

writeLines(
  note_lines,
  con = file.path(base_dir, "selection", "germination_notes.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(base_dir, "selection", "sessionInfo.txt")
)

############################
## 13. Finish
############################
cat("\nDone.\n")
cat("Saved:\n")
cat(normalizePath(file.path(base_dir, "selection", "germination_priors.rds")), "\n")

############################
## 14. checking
############################