# Purpose:
# Regrid fine-resolution observed monthly precipitation data
# to each coarse GCM grid using area-averaging of fine pixels
# within coarse GCM cells.
#
# This script creates several observation products. Some are
# used directly in hindcast validation, while others are
# produced for different parts of the workflow or for checking.
#
# ------------------------------------------------------------
# OUTPUT PRODUCTS AND THEIR USES
# ------------------------------------------------------------
#
# 1. Shared 250 m observed climatology
#
#    File:
#      precip_obs_climatology_1986-2014_250m.tif
#
#    Used in this regridding script:
#      - Yes, internally, as the 12-month climatology that is
#        also aggregated to each GCM grid.
#
#    Used in hindcast validation:
#      - No.
#
#    Used in the current delta downscaling script:
#      - No. The current delta downscaling script uses the
#        original observed climatology NetCDF directly:
#        precip_3ds_normals_19862015.nc
#
#    Purpose:
#      - Provides a shared 250 m GeoTIFF version of the observed
#        precipitation climatology.
#      - Useful for checking, mapping, and optional workflows.
#
# ------------------------------------------------------------
#
# 2. Monthly observed precipitation on each GCM grid
#
#    File:
#      precip_obs_monthly_on_<MODEL>_grid.tif
#
#    Used in hindcast validation:
#      - Yes.
#
#    Used in the current delta downscaling script:
#      - No.
#
#    Purpose:
#      - This is the main observation input for hindcast
#        validation.
#      - It is compared directly with aligned historical GCM
#        monthly time series at native GCM resolution.
#      - Therefore, hindcast validation is performed at coarse
#        GCM resolution, not at 250 m.
#
# ------------------------------------------------------------
#
# 3. Observed precipitation climatology on each GCM grid
#
#    File:
#      precip_obs_climatology_1986-2014_on_<MODEL>_grid.tif
#
#    Used in hindcast validation:
#      - No.
#
#    Used in the current delta downscaling script:
#      - No.
#
#    Purpose:
#      - Produced for optional coarse-grid bias-correction,
#        climatology comparison, and diagnostic workflows.
#      - Useful when an observed monthly climatology is required
#        on the exact same grid as each GCM.
#
# ------------------------------------------------------------
# IMPORTANT METHOD NOTE
# ------------------------------------------------------------
#
# Fine-resolution observations are not directly resampled to the
# GCM grid using simple interpolation. Instead:
#
#   1. Observations are projected to the GCM CRS.
#   2. Fine-resolution pixels are area-averaged within each
#      coarse GCM cell.
#
# This preserves the spatial support of the fine-resolution
# observations before aggregation to the coarse GCM grid.
#
# Variable:
#   pr / precipitation
#
# ============================================================
#
# ============================================================

suppressPackageStartupMessages({
  library(terra)
})

terraOptions(progress = 1, memfrac = 0.7)

# ============================================================
# 1. SETTINGS
# ============================================================

project_root <- "~/Climate_downscaling/Downscaling"

obs_file <- file.path(
  dirname(project_root),
  "Observation",
  "precip_3ds_monthly_198505202008.nc"
)

gcm_dir <- file.path(project_root, "data", "gcm_monthly")

btn_file <- file.path(
  project_root,
  "data",
  "boundary",
  "Dzongkhag_projected.shp"
)

out_root <- file.path(
  project_root,
  "data",
  "observation_regridded_to_each_gcm_pr"
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

baseline_start <- as.Date("1986-01-01")
baseline_end   <- as.Date("2014-12-31")

target_crs <- "EPSG:5266"
cover_threshold <- 0.30

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

get_model_name <- function(f) {
  parts <- strsplit(basename(f), "_")[[1]]

  if (length(parts) >= 3) {
    return(parts[3])
  }

  NA_character_
}

choose_template_file <- function(model_name, files) {
  hits <- files[grepl(paste0("_", model_name, "_"), basename(files))]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  hist_hits <- hits[grepl("historical", basename(hits), ignore.case = TRUE)]

  if (length(hist_hits) > 0) {
    return(hist_hits[1])
  }

  hits[1]
}

aggregate_obs_to_gcm_cells <- function(obs_layer_proj, gcm_polys, gcm_window) {
  extracted <- terra::extract(
    obs_layer_proj,
    gcm_polys,
    fun = mean,
    na.rm = TRUE
  )

  out_poly <- gcm_polys
  out_poly$value <- extracted[, 2]

  rasterize(out_poly, gcm_window, field = "value")
}

# ============================================================
# 3. READ OBSERVATION AND BOUNDARY DATA
# ============================================================

obs <- rast(obs_file)
btn <- vect(btn_file)

# Crop and mask observed data to Bhutan
obs <- crop(obs, btn)
obs <- mask(obs, btn)

# Correct CRS metadata if missing or incorrect.
# This assigns the known CRS; it does not reproject the raster.
crs(obs) <- target_crs

cat("\n====================================================\n")
cat("OBSERVED PRECIPITATION STACK\n")
print(obs)

cat("\nBHUTAN BOUNDARY\n")
print(btn)

cat("\nObservation CRS:\n")
print(crs(obs))

# ============================================================
# 4. CHECK OBSERVATION TIME DIMENSION
# ============================================================

obs_time <- time(obs)

if (is.null(obs_time)) {
  stop("Observation stack has no time dimension. Please assign time(obs) first.")
}

cat("\nObservation time range:\n")
print(range(obs_time))

# ============================================================
# 5. CREATE SHARED 250 m OBSERVED CLIMATOLOGY
# ============================================================

base_idx <- which(obs_time >= baseline_start & obs_time <= baseline_end)

if (length(base_idx) == 0) {
  stop("No observation layers found within the baseline period.")
}

obs_base <- obs[[base_idx]]

obs_clim_250 <- tapp(
  obs_base,
  index = format(time(obs_base), "%m"),
  fun   = mean,
  na.rm = TRUE
)

names(obs_clim_250) <- month.abb

shared_dir <- file.path(out_root, "shared_250m")
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)

shared_file <- file.path(
  shared_dir,
  "precip_obs_climatology_1986-2014_250m.tif"
)

writeRaster(obs_clim_250, shared_file, overwrite = TRUE)

cat("\nShared 250 m precipitation climatology written to:\n")
cat(shared_file, "\n")

# ============================================================
# 6. LIST MONTHLY PRECIPITATION GCM FILES
# ============================================================

gcm_pr_files <- list.files(
  gcm_dir,
  pattern = "^pr_.*\\.nc$",
  full.names = TRUE
)

if (length(gcm_pr_files) == 0) {
  stop("No monthly precipitation GCM files found in: ", gcm_dir)
}

model_names <- sort(unique(na.omit(sapply(gcm_pr_files, get_model_name))))

cat("\n====================================================\n")
cat("Number of unique GCM models found:", length(model_names), "\n")
print(model_names)

# ============================================================
# 7. LOOP OVER GCM MODELS
# ============================================================

failed_models <- character(0)

for (mdl in model_names) {

  cat("\n====================================================\n")
  cat("Processing model:", mdl, "\n")

  tryCatch({

    template_file <- choose_template_file(mdl, gcm_pr_files)

    if (is.na(template_file)) {
      stop("No template GCM file found for model: ", mdl)
    }

    cat("Template file:", basename(template_file), "\n")

    # --------------------------------------------------------
    # Read one GCM layer as the model grid template
    # --------------------------------------------------------

    gcm <- rast(template_file)
    gcm_template <- gcm[[1]]

    cat("\nGCM template:\n")
    print(gcm_template)

    # --------------------------------------------------------
    # Project Bhutan boundary to GCM CRS
    # --------------------------------------------------------

    btn_gcm <- project(btn, crs(gcm_template))

    # --------------------------------------------------------
    # Create GCM window covering Bhutan using full coarse cells
    # --------------------------------------------------------

    gcm_window <- crop(gcm_template, ext(btn_gcm), snap = "out")

    if (ncell(gcm_window) == 0) {
      stop("GCM window has zero cells for model: ", mdl)
    }

    cat("\nGCM window:\n")
    print(gcm_window)

    # --------------------------------------------------------
    # Create polygons for each coarse GCM cell
    # --------------------------------------------------------

    gcm_polys_all <- as.polygons(gcm_window, aggregate = FALSE)
    gcm_polys_all$cell_id <- seq_len(nrow(gcm_polys_all))

    # Intersect coarse GCM cells with Bhutan boundary
    gcm_intersections_list <- intersect(gcm_polys_all, btn_gcm)
    gcm_intersections <- do.call(rbind, gcm_intersections_list)

    # Compute full cell area and Bhutan-covered area
    cell_area <- expanse(gcm_polys_all, unit = "km")

    int_area_df <- data.frame(
      cell_id  = gcm_intersections$cell_id,
      int_area = expanse(gcm_intersections, unit = "km")
    )

    int_area_sum <- aggregate(
      int_area ~ cell_id,
      data = int_area_df,
      sum
    )

    gcm_polys_all$cell_area <- cell_area
    gcm_polys_all$int_area  <- 0

    m <- match(gcm_polys_all$cell_id, int_area_sum$cell_id)

    gcm_polys_all$int_area[!is.na(m)] <-
      int_area_sum$int_area[m[!is.na(m)]]

    gcm_polys_all$cover_frac <-
      gcm_polys_all$int_area / gcm_polys_all$cell_area

    # Keep cells with sufficient Bhutan coverage
    gcm_polys <- gcm_polys_all[
      gcm_polys_all$cover_frac >= cover_threshold,
    ]

    cat("\nTotal coarse cells in window:", nrow(gcm_polys_all), "\n")
    cat("Cells kept after coverage filter:", nrow(gcm_polys), "\n")
    cat("Coverage threshold:", cover_threshold, "\n")

    if (nrow(gcm_polys) == 0) {
      stop("No GCM cells retained after coverage filtering for model: ", mdl)
    }

    # --------------------------------------------------------
    # Project observations to GCM CRS only
    #
    # Important:
    # The fine-resolution observation data are not directly
    # resampled to the GCM grid. They are first projected to the
    # GCM CRS, then area-averaged inside each GCM cell.
    # --------------------------------------------------------

    cat("\nProjecting monthly observations to GCM CRS...\n")
    obs_proj <- project(obs, crs(gcm_template), method = "bilinear")

    cat("Projecting 12-layer observed climatology to GCM CRS...\n")
    obs_clim_proj <- project(obs_clim_250, crs(gcm_template), method = "bilinear")

    # Crop projected observations to the GCM window for speed
    obs_proj <- crop(obs_proj, ext(gcm_window), snap = "out")
    obs_clim_proj <- crop(obs_clim_proj, ext(gcm_window), snap = "out")

    # --------------------------------------------------------
    # Create model-specific output folder
    # --------------------------------------------------------

    mdl_dir <- file.path(out_root, mdl)
    dir.create(mdl_dir, recursive = TRUE, showWarnings = FALSE)

    # --------------------------------------------------------
    # Output 1: monthly observed precipitation on GCM grid
    # --------------------------------------------------------

    cat("\nAggregating monthly observations to GCM cells...\n")

    obs_on_gcm_list <- vector("list", nlyr(obs_proj))

    for (i in seq_len(nlyr(obs_proj))) {
      cat("  Layer", i, "of", nlyr(obs_proj), "\n")

      obs_on_gcm_list[[i]] <- aggregate_obs_to_gcm_cells(
        obs_layer_proj = obs_proj[[i]],
        gcm_polys      = gcm_polys,
        gcm_window     = gcm_window
      )
    }

    obs_on_gcm <- rast(obs_on_gcm_list)

    time(obs_on_gcm) <- obs_time
    names(obs_on_gcm) <- paste0("precip_obs_", format(obs_time, "%Y%m"))

    file_monthly <- file.path(
      mdl_dir,
      paste0("precip_obs_monthly_on_", mdl, "_grid.tif")
    )

    writeRaster(obs_on_gcm, file_monthly, overwrite = TRUE)

    # --------------------------------------------------------
    # Output 2: observed precipitation climatology on GCM grid
    # --------------------------------------------------------

    cat("\nAggregating 12-layer climatology to GCM cells...\n")

    obs_clim_gcm_list <- vector("list", nlyr(obs_clim_proj))

    for (i in seq_len(nlyr(obs_clim_proj))) {
      cat("  Month", i, "of", nlyr(obs_clim_proj), "\n")

      obs_clim_gcm_list[[i]] <- aggregate_obs_to_gcm_cells(
        obs_layer_proj = obs_clim_proj[[i]],
        gcm_polys      = gcm_polys,
        gcm_window     = gcm_window
      )
    }

    obs_clim_gcm <- rast(obs_clim_gcm_list)
    names(obs_clim_gcm) <- month.abb

    file_clim <- file.path(
      mdl_dir,
      paste0("precip_obs_climatology_1986-2014_on_", mdl, "_grid.tif")
    )

    writeRaster(obs_clim_gcm, file_clim, overwrite = TRUE)

    # --------------------------------------------------------
    # Geometry and value checks
    # --------------------------------------------------------

    cat("\nGeometry check: monthly observation vs GCM window\n")
    print(compareGeom(obs_on_gcm[[1]], gcm_window, stopOnError = FALSE))

    cat("\nGeometry check: climatology vs GCM window\n")
    print(compareGeom(obs_clim_gcm[[1]], gcm_window, stopOnError = FALSE))

    rr1 <- global(obs_on_gcm[[1]], range, na.rm = TRUE)
    rr2 <- global(obs_clim_gcm[[1]], range, na.rm = TRUE)

    cat("\nSaved files:\n")
    cat("  ", file_monthly, "\n")
    cat("  ", file_clim, "\n")

    cat("\nFirst monthly layer range:",
        rr1[1, 1], "to", rr1[1, 2], "\n")

    cat("January climatology range:",
        rr2[1, 1], "to", rr2[1, 2], "\n")

    # --------------------------------------------------------
    # Clean memory
    # --------------------------------------------------------

    rm(
      gcm, gcm_template, btn_gcm, gcm_window,
      gcm_polys_all, gcm_intersections,
      int_area_df, int_area_sum, m, gcm_polys,
      obs_proj, obs_clim_proj,
      obs_on_gcm_list, obs_on_gcm,
      obs_clim_gcm_list, obs_clim_gcm
    )

    gc()

  }, error = function(e) {

    failed_models <<- c(failed_models, mdl)

    cat("\nERROR while processing model:", mdl, "\n")
    cat("Message:", conditionMessage(e), "\n")

  })
}

# ============================================================
# 8. FINAL SUMMARY
# ============================================================

cat("\n====================================================\n")
cat("Observed precipitation regridding completed.\n")
cat("Output directory:\n", out_root, "\n")
cat("Shared 250 m climatology:\n", shared_file, "\n")

if (length(failed_models) > 0) {
  cat("\nModels that failed:\n")
  print(failed_models)
} else {
  cat("\nAll models processed successfully.\n")
}

cat("====================================================\n")
