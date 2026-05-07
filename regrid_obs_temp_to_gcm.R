# Purpose:
# Regrid fine-resolution observed monthly temperature data
# to each coarse GCM grid using area-averaging of fine pixels
# within coarse GCM cells.
#
# This script can be used for:
#   - tasmax / maximum temperature
#   - tasmin / minimum temperature
#
# The variable is controlled by:
#   var_name <- "tasmax"
#   var_name <- "tasmin"
#
# This script creates several observation products. Some are
# used directly in hindcast validation, while others are
# produced for checking or optional workflows.
#
# ------------------------------------------------------------
# OUTPUT PRODUCTS AND THEIR USES
# ------------------------------------------------------------
#
# 1. Shared 250 m observed temperature climatology
#
#    Files:
#      tasmax_obs_climatology_1986-2014_250m.tif
#      tasmin_obs_climatology_1986-2014_250m.tif
#
#    Used in hindcast validation:
#      - No.
#
#    Used in the current delta downscaling script:
#      - No. The current delta downscaling script uses the
#        original observed climatology NetCDF files directly:
#        tmax_3ds_m_normals_19862015.nc
#        tmin_3ds_tm_normals_19862015.nc
#
#    Purpose:
#      - Provides shared 250 m GeoTIFF versions of the observed
#        temperature climatologies.
#      - Useful for checking, mapping, and optional workflows.
#
# ------------------------------------------------------------
#
# 2. Monthly observed temperature on each GCM grid
#
#    Files:
#      tasmax_obs_monthly_on_<MODEL>_grid.tif
#      tasmin_obs_monthly_on_<MODEL>_grid.tif
#
#    Used in hindcast validation:
#      - Yes.
#
#    Used in the current delta downscaling script:
#      - No.
#
#    Purpose:
#      - Main observation input for tasmax/tasmin hindcast
#        validation.
#      - Compared directly with aligned historical GCM monthly
#        time series at native GCM resolution.
#      - Therefore, hindcast validation is performed at coarse
#        GCM resolution, not at 250 m.
#
# ------------------------------------------------------------
#
# 3. Observed temperature climatology on each GCM grid
#
#    Files:
#      tasmax_obs_climatology_1986-2014_on_<MODEL>_grid.tif
#      tasmin_obs_climatology_1986-2014_on_<MODEL>_grid.tif
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
#      - Useful when observed monthly temperature climatology
#        is required on the exact same grid as each GCM.
#
# ------------------------------------------------------------
# IMPORTANT METHOD NOTE
# ------------------------------------------------------------
#
# Fine-resolution observations are not directly resampled to
# the GCM grid using simple interpolation. Instead:
#
#   1. Observations are projected to the GCM CRS.
#   2. Fine-resolution pixels are area-averaged within each
#      coarse GCM cell.
#
# This preserves the spatial support of the fine-resolution
# observations before aggregation to the coarse GCM grid.
#
# Variables:
#   tasmax / maximum temperature
#   tasmin / minimum temperature
#
# ============================================================
# 1. SETTINGS
# ============================================================

project_root <- "~/Climate_downscaling/Downscaling"

var_name <- "tasmax"   # use "tasmax" or "tasmin"

obs_file <- file.path(
  dirname(project_root),
  "Observation",
  "tmax_3ds_m_monthly_198505202008.nc"
)

gcm_dir <- file.path(
  project_root,
  "data",
  "gcm_monthly"
)

btn_file <- file.path(
  project_root,
  "data",
  "boundary",
  "Dzongkhag_projected.shp"
)

out_root <- file.path(
  project_root,
  "data",
  "observation_regridded_to_each_gcm_tmax"
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

baseline_start <- as.Date("1986-01-01")
baseline_end   <- as.Date("2014-12-31")

# Observation files are in Bhutan projected metres
obs_source_crs <- "EPSG:5266"


# ============================================================
# 2. READ INPUTS
# ============================================================

obs <- rast(obs_file)
btn <- vect(btn_file)

# Force correct source CRS for observation
crs(obs) <- obs_source_crs

cat("====================================================\n")
cat("VARIABLE:", var_name, "\n")
cat("OBSERVATION STACK\n")
print(obs)

cat("\nBOUNDARY\n")
print(btn)

cat("\nOBS CRS AFTER FORCING\n")
print(crs(obs))

obs_time <- time(obs)
if (is.null(obs_time)) {
  stop("Observation stack has no time dimension.")
}

cat("\nObservation time range:\n")
print(range(obs_time))

# ============================================================
# 3. BASELINE SUBSET + SHARED 250 m CLIMATOLOGY
# ============================================================

base_idx <- which(obs_time >= baseline_start & obs_time <= baseline_end)

if (length(base_idx) == 0) {
  stop("No observation layers found in baseline period.")
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
  paste0(var_name, "_obs_climatology_1986-2014_250m.tif")
)

writeRaster(obs_clim_250, shared_file, overwrite = TRUE)

cat("\nShared 250 m climatology written to:\n", shared_file, "\n")

# ============================================================
# 4. LIST GCM FILES
# ============================================================

gcm_files <- list.files(
  gcm_dir,
  pattern = paste0("^", var_name, "_.*\\.nc$"),
  full.names = TRUE
)

if (length(gcm_files) == 0) {
  stop("No matching GCM files found for variable: ", var_name)
}

get_model_name <- function(f) {
  parts <- strsplit(basename(f), "_")[[1]]
  if (length(parts) >= 3) parts[3] else NA_character_
}

model_names <- sort(unique(na.omit(sapply(gcm_files, get_model_name))))

cat("\n====================================================\n")
cat("NUMBER OF UNIQUE MODELS FOUND:", length(model_names), "\n")
print(model_names)

choose_template_file <- function(model_name, files) {
  hits <- files[grepl(paste0("_", model_name, "_"), basename(files))]
  if (length(hits) == 0) return(NA_character_)
  
  hist_hits <- hits[grepl("historical", basename(hits), ignore.case = TRUE)]
  if (length(hist_hits) > 0) return(hist_hits[1])
  
  hits[1]
}

# ============================================================
# 5. HELPER: AGGREGATE FINE OBS TO COARSE GCM CELLS
# ============================================================

aggregate_obs_to_gcm_cells <- function(obs_layer_proj, gcm_polys, gcm_window) {
  ex <- extract(obs_layer_proj, gcm_polys, fun = mean, na.rm = TRUE)
  
  out_poly <- gcm_polys
  out_poly$value <- ex[, 2]
  
  out_rast <- rasterize(out_poly, gcm_window, field = "value")
  out_rast
}

# ============================================================
# 6. LOOP OVER MODELS
# ============================================================

for (mdl in model_names) {
  
  cat("\n====================================================\n")
  cat("PROCESSING MODEL:", mdl, "\n")
  
  template_file <- choose_template_file(mdl, gcm_files)
  
  if (is.na(template_file)) {
    cat("  No template file found. Skipping.\n")
    next
  }
  
  cat("  Template file:", basename(template_file), "\n")
  
  # ----------------------------------------------------------
  # Read one GCM layer as template
  # ----------------------------------------------------------
  gcm <- rast(template_file)
  gcm_template <- gcm[[1]]
  
  cat("  GCM template summary:\n")
  print(gcm_template)
  
  # ----------------------------------------------------------
  # Project Bhutan boundary to GCM CRS
  # ----------------------------------------------------------
  btn_gcm <- project(btn, crs(gcm_template))
  
  # ----------------------------------------------------------
  # Build GCM window over Bhutan using full coarse cells
  # ----------------------------------------------------------
  gcm_window <- crop(gcm_template, ext(btn_gcm), snap = "out")
  
  if (ncell(gcm_window) == 0) {
    cat("  GCM window has zero cells. Skipping.\n")
    next
  }
  
  cat("  GCM window summary:\n")
  print(gcm_window)
  
  # ----------------------------------------------------------
  # One polygon per coarse GCM cell
  # ----------------------------------------------------------
  gcm_polys <- as.polygons(gcm_window, aggregate = FALSE)
  gcm_polys$cell_id <- seq_len(nrow(gcm_polys))
  
  # ----------------------------------------------------------
  # Project observations to GCM CRS only
  # ----------------------------------------------------------
  cat("  Projecting full monthly observation stack to GCM CRS only...\n")
  obs_proj <- project(obs, crs(gcm_template), method = "bilinear")
  
  cat("  Projecting 12-layer climatology to GCM CRS only...\n")
  obs_clim_proj <- project(obs_clim_250, crs(gcm_template), method = "bilinear")
  
  # Optional crop for speed
  obs_proj <- crop(obs_proj, ext(gcm_window), snap = "out")
  obs_clim_proj <- crop(obs_clim_proj, ext(gcm_window), snap = "out")
  
  cat("  Projected monthly obs summary:\n")
  print(obs_proj)
  
  # ----------------------------------------------------------
  # Output folder
  # ----------------------------------------------------------
  mdl_dir <- file.path(out_root, mdl)
  dir.create(mdl_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # FILE 1: monthly observed series on this GCM grid
  # ----------------------------------------------------------
  cat("  Aggregating full monthly observation stack to coarse GCM cells...\n")
  
  obs_on_gcm_list <- vector("list", nlyr(obs_proj))
  
  for (i in seq_len(nlyr(obs_proj))) {
    cat("    layer", i, "of", nlyr(obs_proj), "\n")
    obs_on_gcm_list[[i]] <- aggregate_obs_to_gcm_cells(
      obs_layer_proj = obs_proj[[i]],
      gcm_polys      = gcm_polys,
      gcm_window     = gcm_window
    )
  }
  
  obs_on_gcm <- rast(obs_on_gcm_list)
  time(obs_on_gcm) <- obs_time
  names(obs_on_gcm) <- paste0(var_name, "_obs_", seq_len(nlyr(obs_on_gcm)))
  
  file_monthly <- file.path(
    mdl_dir,
    paste0(var_name, "_obs_monthly_on_", mdl, "_grid.tif")
  )
  
  writeRaster(obs_on_gcm, file_monthly, overwrite = TRUE)
  
  # ----------------------------------------------------------
  # FILE 2: 12-layer climatology on this GCM grid
  # ----------------------------------------------------------
  cat("  Aggregating 12-layer climatology to coarse GCM cells...\n")
  
  obs_clim_gcm_list <- vector("list", nlyr(obs_clim_proj))
  
  for (i in seq_len(nlyr(obs_clim_proj))) {
    cat("    month", i, "of", nlyr(obs_clim_proj), "\n")
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
    paste0(var_name, "_obs_climatology_1986-2014_on_", mdl, "_grid.tif")
  )
  
  writeRaster(obs_clim_gcm, file_clim, overwrite = TRUE)
  
  # ----------------------------------------------------------
  # Checks
  # ----------------------------------------------------------
  cat("  Geometry check monthly vs template:\n")
  print(compareGeom(obs_on_gcm[[1]], gcm_window, stopOnError = FALSE))
  
  cat("  Geometry check climatology vs template:\n")
  print(compareGeom(obs_clim_gcm[[1]], gcm_window, stopOnError = FALSE))
  
  rr1 <- global(obs_on_gcm[[1]], range, na.rm = TRUE)
  rr2 <- global(obs_clim_gcm[[1]], range, na.rm = TRUE)
  
  cat("  Saved files:\n")
  cat("   ", file_monthly, "\n")
  cat("   ", file_clim, "\n")
  cat("  First monthly layer range :", rr1[1, 1], "to", rr1[1, 2], "\n")
  cat("  January climatology range :", rr2[1, 1], "to", rr2[1, 2], "\n")
  
  # ----------------------------------------------------------
  # Clean memory
  # ----------------------------------------------------------
  rm(gcm, gcm_template, btn_gcm, gcm_window, gcm_polys,
     obs_proj, obs_clim_proj,
     obs_on_gcm_list, obs_on_gcm,
     obs_clim_gcm_list, obs_clim_gcm)
  gc()
}

cat("\n====================================================\n")
cat("All observation products created for:", var_name, "\n")
cat("Shared 250 m climatology:\n", shared_file, "\n")
