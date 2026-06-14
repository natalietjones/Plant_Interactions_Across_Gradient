############################################################
## WA Wheatbelt map:
## 30-year April–November growing-season rainfall classes
##
## Path update:
## - This script can live in raw_data_summaries_map/.
## - It reads site coordinates from ../data/field/plot_location_data.csv
##   relative to the project root located by walking upward.
## - It downloads/writes rainfall files to ../data/rainfall/.
## - It writes plots to ../figures/.
############################################################

rm(list = ls())

############################################################
## 1. Packages
############################################################

packages <- c(
  "terra",
  "sf",
  "dplyr",
  "ggplot2",
  "ggspatial",
  "ozmaps",
  "readr",
  "purrr",
  "scales",
  "ggrepel",
  "cowplot",
  "tibble",
  "grid"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed], dependencies = TRUE)

invisible(lapply(packages, library, character.only = TRUE))

############################################################
## 1b. Robust project paths
############################################################

## The script is intended to live in raw_data_summaries_map/.
## The data and figures folders are expected to be one level above that folder:
##   ../data/field/plot_location_data.csv
##   ../data/rainfall/
##   ../figures/
##
## This root-finder is deliberately robust to whether you run the script from:
##   - raw_data_summaries_map/
##   - the project root
##   - another nearby folder

SITE_FILE_NAME <- "plot_location_data.csv"

normalise_existing <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

find_project_root <- function(site_file_name = SITE_FILE_NAME) {
  start_dirs <- unique(normalise_existing(c(
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )))
  start_dirs <- start_dirs[dir.exists(start_dirs)]

  for (start in start_dirs) {
    current <- start

    for (i in seq_len(10)) {
      field_file <- file.path(current, "data", "field", site_file_name)
      field_dir  <- file.path(current, "data", "field")
      figures_dir <- file.path(current, "figures")

      if (file.exists(field_file) || (dir.exists(field_dir) && dir.exists(figures_dir))) {
        return(normalizePath(current, winslash = "/", mustWork = TRUE))
      }

      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  stop(
    "Could not find the project root. Expected to find data/field/", site_file_name, "\n",
    "Starting from getwd(): ", getwd(), "\n",
    "Check that the project has data/field/ and figures/ folders."
  )
}

project_dir <- find_project_root(SITE_FILE_NAME)

field_data_dir <- file.path(project_dir, "data", "field")
rainfall_dir   <- file.path(project_dir, "data", "rainfall")
figures_dir    <- file.path(project_dir, "figures")

site_file <- file.path(field_data_dir, SITE_FILE_NAME)

## Create output folders only. Do not create data/field, because that should already exist.
dir.create(rainfall_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir,  showWarnings = FALSE, recursive = TRUE)

if (!file.exists(site_file)) {
  available_field_files <- if (dir.exists(field_data_dir)) {
    list.files(field_data_dir, full.names = FALSE)
  } else {
    character()
  }

  stop(
    "Cannot find site location file: ", site_file, "\n",
    "Files currently available in data/field are:\n",
    paste(available_field_files, collapse = "\n")
  )
}

cat("\nPath setup:\n")
cat("  project_dir:  ", normalizePath(project_dir, winslash = "/"), "\n", sep = "")
cat("  site_file:    ", normalizePath(site_file, winslash = "/"), "\n", sep = "")
cat("  rainfall_dir: ", normalizePath(rainfall_dir, winslash = "/"), "\n", sep = "")
cat("  figures_dir:  ", normalizePath(figures_dir, winslash = "/"), "\n\n", sep = "")

############################################################
## 2. Map extent and constants
############################################################

## Same extent as your old script
xmin <- 115.32739508354608
xmax <- 119.68146616967087
ymin <- -32.70626925914411
ymax <- -29.035739615127557

map_ext <- terra::ext(xmin, xmax, ymin, ymax)

## 30-year period
## Note: 1991:2020 is a 30-year period ending in 2020.
years_use <- 1991:2020

## Growing season: April to November
months_use <- 4:11

############################################################
## 3. Site locations
############################################################

site_raw <- readr::read_csv(site_file, show_col_types = FALSE)

## Try to find common column-name variants.
find_first_col <- function(dat, candidates, label) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) {
    stop(
      "Could not find ", label, " column. Tried: ", paste(candidates, collapse = ", "), "\n",
      "Columns available are:\n",
      paste(names(dat), collapse = ", ")
    )
  }
  hit[[1]]
}

site_name_col <- find_first_col(
  site_raw,
  candidates = c("site", "Site", "site_label", "Site_label", "site_name", "Site_name", "site_num", "Site_num"),
  label = "site label"
)

lat_col <- find_first_col(
  site_raw,
  candidates = c("lat", "Lat", "latitude", "Latitude", "LAT"),
  label = "latitude"
)

lon_col <- find_first_col(
  site_raw,
  candidates = c("long", "Long", "lon", "Lon", "longitude", "Longitude", "LONG", "LON"),
  label = "longitude"
)

sites_map <- site_raw %>%
  dplyr::transmute(
    site_label = as.character(.data[[site_name_col]]),
    lat        = as.numeric(.data[[lat_col]]),
    long       = as.numeric(.data[[lon_col]])
  ) %>%
  dplyr::filter(
    !is.na(site_label),
    !is.na(lat),
    !is.na(long)
  ) %>%
  dplyr::group_by(site_label) %>%
  dplyr::summarise(
    lat  = mean(lat, na.rm = TRUE),
    long = mean(long, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(sites_map) == 0) {
  stop("No valid site coordinate rows after reading: ", site_file)
}

sites_map <- sites_map %>%
  dplyr::mutate(
    label_long = long,
    label_lat  = lat,

    ## Shift Perenjori label down to reduce overlap.
    label_lat = dplyr::if_else(
      site_label %in% c("Perenjori", "PJ", "Perenjori "),
      lat - 0.12,
      label_lat
    )
  )

sites_sf <- sf::st_as_sf(
  sites_map,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)

## Perth point
perth <- tibble::tibble(
  city = "Perth",
  lon = 115.857048,
  lat = -31.953512
)

perth_sf <- sf::st_as_sf(
  perth,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

############################################################
## 4. Rainfall data: SILO monthly rainfall NetCDF
############################################################

download_silo_monthly_rain <- function(year, out_dir = rainfall_dir) {
  url <- sprintf(
    "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/monthly_rain/%s.monthly_rain.nc",
    year
  )

  dest <- file.path(out_dir, basename(url))

  if (!file.exists(dest)) {
    message("Downloading: ", year)
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  } else {
    message("Already downloaded: ", year)
  }

  dest
}

rain_files <- purrr::map_chr(years_use, download_silo_monthly_rain)

############################################################
## 5. Process rainfall
## For each year:
## - read monthly rainfall
## - crop to WA Wheatbelt extent
## - sum April–November
## Then:
## - average across 1991–2020
############################################################

calc_growing_season_total <- function(nc_file) {
  r <- terra::rast(nc_file)

  ## Usually this is 12 layers, one per month.
  if (terra::nlyr(r) < max(months_use)) {
    stop("Expected at least ", max(months_use), " monthly layers in: ", nc_file)
  }

  r <- r[[months_use]]
  r <- terra::crop(r, map_ext)

  ## Annual growing-season rainfall total
  terra::app(r, fun = sum, na.rm = TRUE)
}

rain_season_stack <- purrr::map(rain_files, calc_growing_season_total)
rain_season_stack <- terra::rast(rain_season_stack)

names(rain_season_stack) <- paste0("rain_", years_use)

## 30-year average growing-season rainfall in mm
rain_30yr_mean <- terra::app(rain_season_stack, mean, na.rm = TRUE)
names(rain_30yr_mean) <- "rain_mm"

rain_raster_file <- file.path(rainfall_dir, "rainfall_apr_nov_mean_1991_2020_mm.tif")
terra::writeRaster(
  rain_30yr_mean,
  rain_raster_file,
  overwrite = TRUE
)

############################################################
## 6. Convert rainfall raster to data frame
## Fixed rainfall classes to highlight site differences
############################################################

rain_df <- as.data.frame(rain_30yr_mean, xy = TRUE, na.rm = TRUE) %>%
  dplyr::rename(
    lon = x,
    lat = y,
    rain_mm = rain_mm
  )

print(range(rain_df$rain_mm, na.rm = TRUE))

## Custom rainfall classes:
## no <140 class; narrow bins where sites occur; >400 grouped
rain_breaks <- c(
  140,
  175,
  215,
  240,
  280,
  400,
  Inf
)

rain_labels <- c(
  "140–175",
  "176–215",
  "216–240",
  "241–280",
  "281–400",
  ">400"
)

rain_df <- rain_df %>%
  dplyr::mutate(
    rain_class = cut(
      rain_mm,
      breaks = rain_breaks,
      labels = rain_labels,
      include.lowest = TRUE,
      right = TRUE
    )
  )

## Cool rainfall colour scheme: light = drier, dark = wetter.
rain_cols <- c(
  "140–175" = "#F0F9E8",
  "176–215" = "#BAE4BC",
  "216–240" = "#7BCCC4",
  "241–280" = "#43A2CA",
  "281–400" = "#0868AC",
  ">400"    = "#084081"
)

print(table(rain_df$rain_class, useNA = "ifany"))

############################################################
## 7. Base map layers
############################################################

states <- ozmaps::ozmap_states

wa_state <- states %>%
  dplyr::filter(NAME == "Western Australia")

aus <- ozmaps::ozmap_country

inset_box <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    crs = sf::st_crs(4326)
  )
)

############################################################
## 8. Main map: 30-year mean rainfall in mm
############################################################

p_main <- ggplot2::ggplot() +
  ggplot2::geom_raster(
    data = rain_df,
    ggplot2::aes(x = lon, y = lat, fill = rain_class)
  ) +
  ggplot2::geom_sf(
    data = wa_state,
    fill = NA,
    colour = "black",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = sites_map,
    ggplot2::aes(x = long, y = lat),
    shape = 19,
    size = 2.5,
    stroke = 0.7,
    colour = "black",
    fill = "black"
  ) +
  ggplot2::geom_text(
    data = sites_map,
    ggplot2::aes(x = label_long, y = label_lat, label = site_label),
    hjust = -0.1,
    vjust = -0.1,
    size = 4,
    colour = "black",
    fontface = "bold"
  ) +
  ggplot2::geom_point(
    data = perth,
    ggplot2::aes(x = lon, y = lat),
    shape = 8,
    size = 3.2,
    colour = "white",
    stroke = 1
  ) +
  ggplot2::geom_text(
    data = perth,
    ggplot2::aes(x = lon, y = lat, label = city),
    hjust = -0.35,
    vjust = -0.35,
    size = 4,
    colour = "white",
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    name = "Mean Apr–Nov\nrainfall (mm)",
    values = rain_cols,
    drop = FALSE
  ) +
  ggplot2::coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    crs = 4326,
    expand = FALSE
  ) +
  ggspatial::annotation_north_arrow(
    which_north = "true",
    location = "tr",
    pad_x = grid::unit(0.5, "cm"),
    pad_y = grid::unit(0.2, "cm"),
    height = grid::unit(1.0, "cm"),
    width = grid::unit(0.75, "cm")
  ) +
  ggspatial::annotation_scale(
    location = "br",
    width_hint = 0.35,
    pad_x = grid::unit(0.2, "cm"),
    pad_y = grid::unit(0.2, "cm")
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 11, face = "bold", colour = "black"),
    axis.title = ggplot2::element_text(size = 12, face = "bold", colour = "black"),

    ## Legend outside plot panel, on right.
    legend.position = "right",
    legend.justification = "center",
    legend.box = "vertical",
    legend.background = ggplot2::element_rect(fill = "white", colour = NA),

    plot.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

print(p_main)

############################################################
## 9. Australia inset map - optional/not used
############################################################

# p_inset <- ggplot2::ggplot() +
#   ggplot2::geom_sf(
#     data = aus,
#     fill = "grey90",
#     colour = "black",
#     linewidth = 0.25
#   ) +
#   ggplot2::geom_sf(
#     data = inset_box,
#     fill = "black",
#     colour = "black",
#     alpha = 0.8
#   ) +
#   ggplot2::coord_sf(crs = 4326, expand = FALSE) +
#   ggplot2::theme_void() +
#   ggplot2::theme(
#     plot.background = ggplot2::element_rect(
#       fill = "white",
#       colour = "black",
#       linewidth = 0.3
#     )
#   )
#
# p_final <- cowplot::ggdraw() +
#   cowplot::draw_plot(p_main) +
#   cowplot::draw_plot(
#     p_inset,
#     x = 0.05,
#     y = 0.67,
#     width = 0.28,
#     height = 0.25
#   )

############################################################
## 10. Save figure
############################################################

main_png <- file.path(figures_dir, "map_wa_30yr_mean_apr_nov_rainfall_mm.png")

ggplot2::ggsave(
  filename = main_png,
  plot = p_main,
  width = 7,
  height = 8.5,
  dpi = 300
)

# main_pdf <- file.path(figures_dir, "map_wa_30yr_mean_apr_nov_rainfall_mm.pdf")
# ggplot2::ggsave(
#   filename = main_pdf,
#   plot = p_main,
#   width = 7,
#   height = 8.5
# )

cat("\nFiles written:\n")
cat("  rainfall raster: ", normalizePath(rain_raster_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  figure:          ", normalizePath(main_png, winslash = "/", mustWork = FALSE), "\n", sep = "")
