############################################################
## Environmental ordination of WA sites
## PCA of soil, canopy and climate variables
##
## Path update:
## - Script can live in scripts/raw_data_summaries_map/
## - Reads data from data/field/fullWAdata.csv
## - Writes figures/tables to figures/raw_data_summaries_map/
############################################################

rm(list = ls())

############################################################
## 0. Packages
############################################################
required_packages <- c(
  "tidyverse",
  "vegan",
  "ggrepel",
  "here",
  "readr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

############################################################
## 1. User settings
############################################################

## File is expected to be in data/field/
## Change this to "plot_location_data.csv" only if that file contains all env_vars below.
DATA_FILE_NAME <- "fullWAdata.csv"

SITE_LEVELS <- c("BEN", "GH", "NAM", "PJ", "CS")
COV_LEVELS  <- c("Sun", "Shade")

env_vars <- c(
  "CEC",
  "pH",
  "NO3.N",
  "mean_percent_open",
  "precip",
  "pet",
  "water_in_out"
)

## Site colours
site_colors <- c(
  "BEN" = "#89A6FB",
  "GH"  = "#3AB795",
  "NAM" = "#FFCF56",
  "PJ"  = "#FC7A57",
  "CS"  = "#FE4A49"
)

cover_shapes <- c(
  "Sun" = 1,      # open circle
  "Shade" = 19   # filled circle
)

BASE_FAMILY <- "Arial"
BASE_SIZE <- 11

############################################################
## 2. Robust project paths
############################################################

## Find project root by searching upward from common starting points.
## This is useful if the script is sourced from scripts/raw_data_summaries_map/
## or run from the project root.
find_project_root <- function(data_file_name = DATA_FILE_NAME) {
  start_dirs <- unique(c(
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    normalizePath(here::here(), winslash = "/", mustWork = TRUE)
  ))
  
  for (start in start_dirs) {
    current <- start
    
    for (i in seq_len(8)) {
      target_file <- file.path(current, "data", "field", data_file_name)
      field_dir   <- file.path(current, "data", "field")
      
      if (file.exists(target_file) || dir.exists(field_dir)) {
        return(normalizePath(current, winslash = "/", mustWork = TRUE))
      }
      
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  
  stop(
    "Could not find the project root containing data/field/", data_file_name, ".\n",
    "Current getwd(): ", getwd(), "\n",
    "Current here::here(): ", here::here()
  )
}

project_dir <- find_project_root(DATA_FILE_NAME)

field_data_dir <- file.path(project_dir, "data", "field")
figures_dir    <- file.path(project_dir, "figures", "raw_data_summaries_map")

DATA_FILE <- file.path(field_data_dir, DATA_FILE_NAME)
OUTDIR <- figures_dir

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(DATA_FILE)) {
  available_field_files <- list.files(field_data_dir, full.names = FALSE)
  stop(
    "Cannot find DATA_FILE: ", DATA_FILE, "\n",
    "Files currently available in data/field are:\n",
    paste(available_field_files, collapse = "\n")
  )
}

cat("\nProject dir: ", normalizePath(project_dir, winslash = "/"), "\n", sep = "")
cat("Data file:   ", normalizePath(DATA_FILE, winslash = "/"), "\n", sep = "")
cat("Output dir:  ", normalizePath(OUTDIR, winslash = "/"), "\n\n", sep = "")

############################################################
## 3. Load and clean environmental data
############################################################

dat_raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

missing_vars <- setdiff(c("Site", "Cov", "BLK", env_vars), names(dat_raw))

if (length(missing_vars) > 0) {
  stop(
    "These required columns are missing from the dataset: ",
    paste(missing_vars, collapse = ", "), "\n",
    "Data file used: ", DATA_FILE
  )
}

env_dat <- dat_raw %>%
  dplyr::mutate(
    Site = dplyr::case_when(
      Site %in% c("BEN", "Ben", "ben") ~ "BEN",
      Site %in% c("GH", "Gh", "gh") ~ "GH",
      Site %in% c("NAM", "Nam", "nam") ~ "NAM",
      Site %in% c("PJ", "Pj", "pj") ~ "PJ",
      Site %in% c("CS", "Cs", "cs") ~ "CS",
      TRUE ~ toupper(as.character(Site))
    ),
    Site = factor(Site, levels = SITE_LEVELS),
    Cov = dplyr::case_when(
      Cov %in% c("SUN", "Sun", "sun") ~ "Sun",
      Cov %in% c("SH", "SHADE", "Shade", "shade", "sh") ~ "Shade",
      TRUE ~ as.character(Cov)
    ),
    Cov = factor(Cov, levels = COV_LEVELS),
    BLK = factor(BLK)
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(env_vars),
      ~ if (is.numeric(.x)) .x else readr::parse_number(as.character(.x))
    )
  ) %>%
  dplyr::select(Site, Cov, BLK, dplyr::all_of(env_vars)) %>%
  dplyr::distinct()

if (any(is.na(env_dat$Site))) {
  bad_sites <- dat_raw %>%
    dplyr::mutate(Site_raw = as.character(Site)) %>%
    dplyr::filter(!toupper(Site_raw) %in% SITE_LEVELS) %>%
    dplyr::distinct(Site_raw) %>%
    dplyr::pull(Site_raw)
  stop("Some Site values could not be matched to SITE_LEVELS: ", paste(bad_sites, collapse = ", "))
}

if (any(is.na(env_dat$Cov))) {
  bad_cov <- dat_raw %>%
    dplyr::mutate(Cov_raw = as.character(Cov)) %>%
    dplyr::filter(!Cov_raw %in% c("SUN", "Sun", "sun", "SH", "SHADE", "Shade", "shade", "sh")) %>%
    dplyr::distinct(Cov_raw) %>%
    dplyr::pull(Cov_raw)
  stop("Some Cov values could not be matched to COV_LEVELS: ", paste(bad_cov, collapse = ", "))
}

## If the data are still repeated at plot level, average to Site x Cover x Block.
env_blk <- env_dat %>%
  dplyr::group_by(Site, Cov, BLK) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(env_vars), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(env_vars), ~ ifelse(is.nan(.x), NA_real_, .x))
  )

# ## Check missing data before ordination.
# missing_summary <- env_blk %>%
#   dplyr::summarise(
#     dplyr::across(dplyr::all_of(env_vars), ~ sum(is.na(.x)))
#   ) %>%
#   tidyr::pivot_longer(
#     everything(),
#     names_to = "variable",
#     values_to = "n_missing"
#   )
# 
# print(missing_summary)
# 
# write_csv(
#   missing_summary,
#   file.path(OUTDIR, "environmental_PCA_missing_summary.csv")
# )

## Drop incomplete rows rather than silently imputing.
env_ord <- env_blk %>%
  tidyr::drop_na(dplyr::all_of(env_vars))

if (nrow(env_ord) < nrow(env_blk)) {
  message(
    "Dropped ", nrow(env_blk) - nrow(env_ord),
    " Site x Cover x Block rows with missing environmental data."
  )
}

if (nrow(env_ord) < 3) {
  stop("Fewer than 3 complete rows remain for PCA. Check missing environmental data.")
}

zero_var <- env_ord %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(env_vars), ~ stats::sd(.x, na.rm = TRUE))
  ) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "sd") %>%
  dplyr::filter(is.na(sd) | sd == 0)

if (nrow(zero_var) > 0) {
  stop(
    "These environmental variables have zero or undefined variance and cannot be scaled in PCA: ",
    paste(zero_var$variable, collapse = ", ")
  )
}

############################################################
## 4. Run PCA using vegan::rda
############################################################

env_mat <- env_ord %>%
  dplyr::select(dplyr::all_of(env_vars)) %>%
  as.data.frame()

## Scale = TRUE is essential because variables have different units.
pca_env <- vegan::rda(env_mat, scale = TRUE)

pca_summary <- summary(pca_env)
print(pca_summary)

eig <- vegan::eigenvals(pca_env)
var_explained <- eig / sum(eig)

axis_lab <- paste0(
  c("PC1", "PC2"),
  " (",
  round(100 * var_explained[1:2], 1),
  "%)"
)

print(var_explained)

write_csv(
  tibble::tibble(
    axis = names(var_explained),
    eigenvalue = as.numeric(eig),
    proportion_explained = as.numeric(var_explained)
  ),
  file.path(OUTDIR, "environmental_PCA_variance_explained.csv")
)

############################################################
## 5. Extract site scores and variable loadings
############################################################

as_pc_score_df <- function(score_mat) {
  score_df <- as.data.frame(score_mat)
  if (ncol(score_df) < 2) stop("Scores object has fewer than two axes.")
  names(score_df)[1:2] <- c("PC1", "PC2")
  score_df[, c("PC1", "PC2"), drop = FALSE]
}

site_scores <- vegan::scores(
  pca_env,
  display = "sites",
  choices = 1:2,
  scaling = 2
) %>%
  as_pc_score_df() %>%
  dplyr::bind_cols(env_ord)

var_scores <- vegan::scores(
  pca_env,
  display = "species",
  choices = 1:2,
  scaling = 2
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable")

names(var_scores)[names(var_scores) == names(var_scores)[2]] <- "PC1"
names(var_scores)[names(var_scores) == names(var_scores)[3]] <- "PC2"

## Scale arrows to fit the ordination plot.
pc1_den <- diff(range(var_scores$PC1, na.rm = TRUE))
pc2_den <- diff(range(var_scores$PC2, na.rm = TRUE))

if (pc1_den == 0 || pc2_den == 0) {
  stop("Variable scores have zero range on PC1 or PC2, so arrows cannot be scaled.")
}

arrow_mult <- min(
  diff(range(site_scores$PC1, na.rm = TRUE)) / pc1_den,
  diff(range(site_scores$PC2, na.rm = TRUE)) / pc2_den
) * 0.75

## Create arrow endpoints, readable labels, and manual nudges in one block.
var_scores <- var_scores %>%
  dplyr::mutate(
    PC1_arrow = PC1 * arrow_mult,
    PC2_arrow = PC2 * arrow_mult,
    
    variable_label = dplyr::case_when(
      variable == "CEC" ~ "CEC",
      variable == "pH" ~ "pH",
      variable == "NO3.N" ~ "NO3-N",
      variable == "mean_percent_open" ~ "% open canopy",
      variable == "precip" ~ "Precipitation",
      variable == "pet" ~ "PET",
      variable == "water_in_out" ~ "Water balance",
      TRUE ~ variable
    ),
    
    ## Base label position, slightly beyond arrow tips.
    label_x = PC1_arrow * 1.12,
    label_y = PC2_arrow * 1.12,
    
    ## Manual label nudges to reduce overlap.
    label_x = dplyr::case_when(
      variable == "water_in_out" ~ label_x + 0.12,
      variable == "pet" ~ label_x - 0.08,
      variable == "precip" ~ label_x + 0.04,
      TRUE ~ label_x
    ),
    label_y = dplyr::case_when(
      variable == "water_in_out" ~ label_y + 0.24,
      variable == "pet" ~ label_y - 0.22,
      variable == "precip" ~ label_y - 0.10,
      variable == "CEC" ~ label_y - 0.04,
      TRUE ~ label_y
    )
  )

############################################################
## 6. Plot PCA
############################################################

## Symmetric y-axis based on site scores and arrow labels.
y_lim <- max(
  abs(site_scores$PC2),
  abs(var_scores$label_y),
  na.rm = TRUE
)
y_lim <- ceiling(y_lim * 2) / 2

fig_env_pca <- ggplot(site_scores, aes(PC1, PC2)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    colour = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.25,
    colour = "grey70"
  ) +
  
  ## Points: site = colour; cover = open vs filled.
  geom_point(
    aes(colour = Site, shape = Cov),
    size = 2.0,
    stroke = 1.2,
    alpha = 0.75
  ) +
  
  ## Environmental vectors.
  geom_segment(
    data = var_scores,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.18, "cm")),
    linewidth = 0.45,
    colour = "grey25"
  ) +
  geom_text(
    data = var_scores,
    aes(
      x = label_x,
      y = label_y,
      label = variable_label
    ),
    inherit.aes = FALSE,
    family = BASE_FAMILY,
    size = 11 / ggplot2::.pt,
    colour = "grey25",
    hjust = 0.5,
    vjust = 0.5
  ) +
  
  scale_colour_manual(
    values = site_colors,
    breaks = SITE_LEVELS,
    name = "Site"
  ) +
  scale_shape_manual(
    values = cover_shapes,
    breaks = COV_LEVELS,
    name = "Cover"
  ) +
  
  ## Symmetric x-axis as requested.
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = seq(-2, 2, by = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-y_lim, y_lim),
    breaks = seq(-ceiling(y_lim), ceiling(y_lim), by = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    x = axis_lab[1],
    y = axis_lab[2]
  ) +
  theme_bw(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
  theme(
    text = element_text(family = BASE_FAMILY),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = grid::unit(1.2, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.35
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(0.55, "lines"),
    legend.key.width = grid::unit(0.75, "lines"),
    plot.margin = margin(3, 3, 3, 3)
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        shape = 19,
        size = 2.0,
        alpha = 0.75
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        colour = "black",
        size = 2.0,
        alpha = 0.75
      )
    )
  )

print(fig_env_pca)

png_file <- file.path(OUTDIR, "environmental_PCA_site_cover.png")


ggsave(
  png_file,
  fig_env_pca,
  width = 4.5,
  height = 4.0,
  dpi = 600
)

############################################################
## 7. Loadings table for interpretation
############################################################

loadings_table <- var_scores %>%
  dplyr::mutate(
    abs_PC1 = abs(PC1),
    abs_PC2 = abs(PC2)
  ) %>%
  dplyr::arrange(desc(abs_PC1))

print(loadings_table)

# write_csv(
#   site_scores,
#   file.path(OUTDIR, "environmental_PCA_site_scores.csv")
# )
# 
# write_csv(
#   var_scores,
#   file.path(OUTDIR, "environmental_PCA_variable_loadings.csv")
# )
# 
# write_csv(
#   loadings_table,
#   file.path(OUTDIR, "environmental_PCA_variable_loadings_ranked.csv")
# )

cat("\nDone. Outputs written to: ", normalizePath(OUTDIR, winslash = "/"), "\n", sep = "")
cat("Figure written to: ", normalizePath(png_file, winslash = "/"), "\n", sep = "")

