# =============================================================================
# Modified delta change downscaling for Bhutan climate data
#
# This script downscales monthly GCM climate fields to a 250 m observational
# grid using a modified delta change approach.
#
# Method:
#   1. Derive GCM monthly climatology from historical simulations (1986-2014)
#   2. Compute monthly change factors at native GCM resolution
#      - precipitation: monthly value / monthly climatology
#      - temperature:   monthly value - monthly climatology
#   3. Interpolate change factors to the 250 m observation grid
#   4. Reconstruct downscaled monthly fields using observation climate normals
#      from 1986-2015
#
# Outputs:
#   Annual GeoTIFF stacks containing 12 monthly layers per year.
#
# -----------------------------------------------------------------------------
# DATA REQUIREMENTS
#
# Users must provide their own input datasets.
#
# The script assumes the following folder structure:
#
# CLIM_ROOT/
# ├── data/
# │   ├── observation/
# │   │   └── nc_files/
# │   │       ├── precip_3ds_normals_19862015.nc
# │   │       ├── tmax_3ds_m_normals_19862015.nc
# │   │       └── tmin_3ds_tm_normals_19862015.nc
# │   └── gcm_monthly/
# │       ├── pr_mon_ACCESS-CM2_historical_1980-2014.nc
# │       ├── pr_mon_ACCESS-CM2_ssp245_2015-2100.nc
# │       └── ...
#
# If your data are stored differently, modify:
#   - obs_root
#   - gcm_root
#
# -----------------------------------------------------------------------------
# HOW TO RUN
#
# Set your project root in R:
#
#   Sys.setenv(CLIM_ROOT = "path/to/your/data")
#
# Then run:
#
#   source("delta_downscaling.R")
#
# =============================================================================

rm(list = ls())

# =============================================================================
# 1. PACKAGE CHECK
# =============================================================================

packages <- c("terra", "stringr")
missing_pkgs <- packages[!packages %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install missing package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
})

terraOptions(progress = 1, memfrac = 0.7)

# =============================================================================
# 2. USER SETTINGS
# =============================================================================

# Variable to process
var_name <- "pr"     # options: "pr", "tasmax", "tasmin"

# Scenarios to process
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

# Run mode
run_mode <- "full"   # options: "full", "test"

# Historical GCM reference climatology period
clim_start_yr <- 1986
clim_end_yr   <- 2014

# Observation climate-normal period
obs_normal_period <- "1986-2015"

# GCM time periods
hist_start_ym <- "1980-01"
fut_start_ym  <- "2015-01"

hist_period_label <- "1980-2014"
fut_period_label  <- "2015-2100"

# Test periods, used only when run_mode = "test"
test_hist_start <- "1980-01"
test_hist_end   <- "1990-12"

test_fut_start  <- "2015-01"
test_fut_end    <- "2025-12"

# Numerical safety threshold for precipitation ratio calculation
pr_floor <- 0.1

# =============================================================================
# 3. PROJECT PATHS
# =============================================================================

# Recommended:
#   Sys.setenv(CLIM_ROOT = "C:/Users/dorji/Documents/Climate_downscaling/Downscaling")
#
# If CLIM_ROOT is not set, the current working directory is used.

root <- Sys.getenv("CLIM_ROOT", unset = getwd())

obs_root <- file.path(root, "data", "observation", "nc_files")
gcm_root <- file.path(root, "data", "gcm_monthly")

out_root <- file.path(root, "outputs", "delta_downscaled", run_mode, var_name)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# Observation climate normals, 250 m, 12 layers Jan-Dec
obs_clim_file <- switch(
  var_name,
  "pr"     = file.path(obs_root, "precip_3ds_normals_19862015.nc"),
  "tasmax" = file.path(obs_root, "tmax_3ds_m_normals_19862015.nc"),
  "tasmin" = file.path(obs_root, "tmin_3ds_tm_normals_19862015.nc"),
  stop("Unsupported var_name: ", var_name)
)

# =============================================================================
# 4. MODEL SETTINGS
# =============================================================================

known_models <- c(
  "ACCESS-CM2",
  "BCC-CSM2-MR",
  "CanESM5",
  "CMCC-ESM2",
  "CNRM-CM6-1",
  "CNRM-ESM2-1",
  "EC-Earth3-Veg-LR",
  "GFDL-ESM4",
  "INM-CM4-8",
  "INM-CM5-0",
  "IPSL-CM6A-LR",
  "MIROC6",
  "MPI-ESM1-2-HR",
  "MPI-ESM1-2-LR",
  "MRI-ESM2-0",
  "NorESM2-MM"
)

# Use NULL to process all available models
# model_keep <- NULL

model_keep <- c(
  "ACCESS-CM2",
  "MIROC6",
  "NorESM2-MM"
)

# =============================================================================
# 5. SAFETY CHECKS
# =============================================================================

if (!dir.exists(obs_root)) {
  stop("Observation folder not found: ", obs_root)
}

if (!dir.exists(gcm_root)) {
  stop("GCM folder not found: ", gcm_root)
}

if (!file.exists(obs_clim_file)) {
  stop("Observation climate-normal file not found: ", obs_clim_file)
}

if (!run_mode %in% c("full", "test")) {
  stop("run_mode must be either 'full' or 'test'")
}

if (!var_name %in% c("pr", "tasmax", "tasmin")) {
  stop("var_name must be one of: pr, tasmax, tasmin")
}

# =============================================================================
# 6. HELPER FUNCTIONS
# =============================================================================

make_monthly_dates <- function(start_ym, n_layers) {
  seq.Date(
    from = as.Date(paste0(start_ym, "-01")),
    by   = "month",
    length.out = n_layers
  )
}

get_model_name <- function(f, known_models) {
  m <- str_match(
    basename(f),
    str_c("(", str_c(known_models, collapse = "|"), ")")
  )
  
  if (!is.na(m[1, 2])) {
    return(m[1, 2])
  }
  
  tools::file_path_sans_ext(basename(f))
}

find_hist_gcm_files <- function(gcm_root, var_name, model_keep = NULL) {
  ff <- list.files(
    gcm_root,
    pattern = "\\.(nc|tif)$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  ff <- ff[str_detect(basename(ff), regex(var_name, ignore_case = TRUE))]
  ff <- ff[str_detect(basename(ff), regex("historical", ignore_case = TRUE))]
  
  if (!is.null(model_keep)) {
    ff <- ff[str_detect(basename(ff), str_c(model_keep, collapse = "|"))]
  }
  
  ff
}

find_future_gcm_files <- function(gcm_root, var_name, scenario_name,
                                  model_keep = NULL) {
  ff <- list.files(
    gcm_root,
    pattern = "\\.(nc|tif)$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  ff <- ff[str_detect(basename(ff), regex(var_name, ignore_case = TRUE))]
  ff <- ff[str_detect(basename(ff), regex(scenario_name, ignore_case = TRUE))]
  ff <- ff[!str_detect(basename(ff), regex("historical", ignore_case = TRUE))]
  
  if (!is.null(model_keep)) {
    ff <- ff[str_detect(basename(ff), str_c(model_keep, collapse = "|"))]
  }
  
  ff
}

compute_gcm_climatology <- function(r_gcm, dates, clim_start_yr, clim_end_yr) {
  yrs <- as.integer(format(dates, "%Y"))
  idx <- which(yrs >= clim_start_yr & yrs <= clim_end_yr)
  
  if (length(idx) == 0) {
    stop(
      sprintf(
        "No GCM layers fall within climatology period %d-%d",
        clim_start_yr,
        clim_end_yr
      )
    )
  }
  
  r_sub  <- r_gcm[[idx]]
  months <- format(dates[idx], "%m")
  
  clim_layers <- lapply(sprintf("%02d", 1:12), function(m) {
    mi <- which(months == m)
    
    if (length(mi) == 0) {
      stop("No layers for month ", m, " in climatology period")
    }
    
    mean(r_sub[[mi]], na.rm = TRUE)
  })
  
  clim <- rast(clim_layers)
  names(clim) <- month.abb
  
  clim
}

compute_gcm_anomaly <- function(r_gcm, gcm_clim, dates, var_name,
                                pr_floor = 0.1) {
  months <- as.integer(format(dates, "%m"))
  anomaly_list <- vector("list", nlyr(r_gcm))
  
  for (i in seq_len(nlyr(r_gcm))) {
    m <- months[i]
    clim_lyr <- gcm_clim[[m]]
    
    if (var_name == "pr") {
      clim_safe <- clamp(clim_lyr, lower = pr_floor, values = TRUE)
      gcm_safe  <- clamp(r_gcm[[i]], lower = pr_floor, values = TRUE)
      
      anom <- gcm_safe / clim_safe
      
    } else {
      anom <- r_gcm[[i]] - clim_lyr
    }
    
    anomaly_list[[i]] <- anom
  }
  
  anom_stack <- rast(anomaly_list)
  names(anom_stack) <- names(r_gcm)
  
  anom_stack
}

delta_reconstruct <- function(anom_stack, obs_clim_hr, dates, var_name) {
  months <- as.integer(format(dates, "%m"))
  downscaled_list <- vector("list", nlyr(anom_stack))
  
  for (i in seq_len(nlyr(anom_stack))) {
    m <- months[i]
    
    # Interpolate anomaly/change factor to the 250 m observation grid
    anom_hr <- project(
      anom_stack[[i]],
      obs_clim_hr[[m]],
      method = "bilinear"
    )
    
    # Reconstruct absolute downscaled field
    if (var_name == "pr") {
      ds <- anom_hr * obs_clim_hr[[m]]
      ds <- ifel(ds < 0, 0, ds)
    } else {
      ds <- anom_hr + obs_clim_hr[[m]]
    }
    
    downscaled_list[[i]] <- ds
  }
  
  ds_stack <- rast(downscaled_list)
  names(ds_stack) <- names(anom_stack)
  
  ds_stack
}

subset_test_period <- function(r, dates, start_ym, end_ym) {
  ym <- format(as.Date(dates), "%Y-%m")
  idx <- which(ym >= start_ym & ym <= end_ym)
  
  if (length(idx) == 0) {
    stop("No layers found in test period: ", start_ym, " to ", end_ym)
  }
  
  list(
    r = r[[idx]],
    dates = dates[idx]
  )
}

save_downscaled_by_year <- function(ds_stack, dates, out_dir, var_name,
                                    model, scenario, period_label) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  yrs <- as.integer(format(dates, "%Y"))
  
  for (yr in unique(yrs)) {
    idx <- which(yrs == yr)
    r_yr <- ds_stack[[idx]]
    
    mo_lbl <- format(dates[idx], "%m")
    names(r_yr) <- paste0(var_name, "_", yr, "_", mo_lbl)
    
    fname <- file.path(
      out_dir,
      sprintf(
        "%s_%s_%s_%s_%d.tif",
        var_name,
        model,
        scenario,
        period_label,
        yr
      )
    )
    
    writeRaster(
      r_yr,
      fname,
      overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6")
    )
    
    cat("    Saved:", basename(fname), "\n")
  }
}

write_processing_metadata <- function(out_root, var_name, scenarios, model_keep,
                                      clim_start_yr, clim_end_yr,
                                      obs_normal_period, run_mode,
                                      hist_period_label, fut_period_label,
                                      obs_clim_file, gcm_root) {
  
  metadata <- data.frame(
    item = c(
      "variable",
      "method",
      "historical_gcm_climatology_period",
      "observation_climate_normal_period",
      "historical_output_period",
      "future_output_period",
      "scenarios",
      "models",
      "run_mode",
      "observation_file",
      "gcm_root",
      "date_processed"
    ),
    value = c(
      var_name,
      "modified_delta_change",
      paste0(clim_start_yr, "-", clim_end_yr),
      obs_normal_period,
      hist_period_label,
      fut_period_label,
      paste(scenarios, collapse = ", "),
      ifelse(is.null(model_keep), "all_available_models", paste(model_keep, collapse = ", ")),
      run_mode,
      obs_clim_file,
      gcm_root,
      as.character(Sys.Date())
    )
  )
  
  write.csv(
    metadata,
    file.path(out_root, "processing_metadata.csv"),
    row.names = FALSE
  )
}

# =============================================================================
# 7. LOAD HIGH-RESOLUTION OBSERVATION CLIMATE NORMALS
# =============================================================================

cat("\n=== Loading high-resolution observation climate normals ===\n")

obs_clim_hr <- rast(obs_clim_file)

# Assign CRS only if missing
if (is.na(crs(obs_clim_hr)) || crs(obs_clim_hr) == "") {
  crs(obs_clim_hr) <- "EPSG:5266"
  cat("  CRS was missing. Assigned EPSG:5266.\n")
} else {
  cat("  Existing CRS detected:\n")
  cat("  ", crs(obs_clim_hr), "\n")
}

if (nlyr(obs_clim_hr) != 12) {
  stop("obs_clim_hr must have exactly 12 layers: Jan-Dec.")
}

names(obs_clim_hr) <- month.abb

cat("  Observation file:", obs_clim_file, "\n")
cat("  Resolution:", paste(res(obs_clim_hr), collapse = " x "), "\n")
cat("  Layers:", nlyr(obs_clim_hr), "\n")
cat("  Observation normal period:", obs_normal_period, "\n")

# =============================================================================
# 8. WRITE METADATA
# =============================================================================

write_processing_metadata(
  out_root = out_root,
  var_name = var_name,
  scenarios = scenarios,
  model_keep = model_keep,
  clim_start_yr = clim_start_yr,
  clim_end_yr = clim_end_yr,
  obs_normal_period = obs_normal_period,
  run_mode = run_mode,
  hist_period_label = hist_period_label,
  fut_period_label = fut_period_label,
  obs_clim_file = obs_clim_file,
  gcm_root = gcm_root
)

# =============================================================================
# 9. PRE-COMPUTE HISTORICAL GCM CLIMATOLOGIES
# =============================================================================

cat("\n=== Computing historical GCM climatologies ===\n")

hist_files <- find_hist_gcm_files(
  gcm_root = gcm_root,
  var_name = var_name,
  model_keep = model_keep
)

if (length(hist_files) == 0) {
  stop("No historical GCM files found for variable: ", var_name)
}

cat("  Number of historical files found:", length(hist_files), "\n")

model_clim_list <- list()

for (f in hist_files) {
  
  model <- get_model_name(f, known_models)
  
  cat("\n--- Historical climatology for model:", model, "---\n")
  cat("  File:", basename(f), "\n")
  
  r_hist <- rast(f)
  hist_dates <- make_monthly_dates(hist_start_ym, nlyr(r_hist))
  
  cat(
    sprintf(
      "  Historical file period assumed: %s to %s (%d layers)\n",
      min(hist_dates),
      max(hist_dates),
      length(hist_dates)
    )
  )
  
  gcm_clim <- compute_gcm_climatology(
    r_gcm = r_hist,
    dates = hist_dates,
    clim_start_yr = clim_start_yr,
    clim_end_yr = clim_end_yr
  )
  
  model_clim_list[[model]] <- gcm_clim
  
  rm(r_hist, hist_dates, gcm_clim)
  gc()
}

# =============================================================================
# 10. DOWNSCALE HISTORICAL MONTHLY SERIES
# =============================================================================

cat("\n=== Downscaling historical monthly series ===\n")

for (f in hist_files) {
  
  model <- get_model_name(f, known_models)
  
  cat("\n--- Model:", model, "| historical ---\n")
  cat("  File:", basename(f), "\n")
  
  r_gcm <- rast(f)
  dates <- make_monthly_dates(hist_start_ym, nlyr(r_gcm))
  
  if (run_mode == "test") {
    sub <- subset_test_period(
      r = r_gcm,
      dates = dates,
      start_ym = test_hist_start,
      end_ym = test_hist_end
    )
    
    r_gcm <- sub$r
    dates <- sub$dates
    
    cat(
      sprintf(
        "  TEST period: %s to %s (%d layers)\n",
        min(dates),
        max(dates),
        length(dates)
      )
    )
    
    period_label <- "test"
    
  } else {
    cat(
      sprintf(
        "  FULL period: %s to %s (%d layers)\n",
        min(dates),
        max(dates),
        length(dates)
      )
    )
    
    period_label <- hist_period_label
  }
  
  gcm_clim <- model_clim_list[[model]]
  
  cat("  Step 1: Computing monthly anomalies/change factors...\n")
  anom_stack <- compute_gcm_anomaly(
    r_gcm = r_gcm,
    gcm_clim = gcm_clim,
    dates = dates,
    var_name = var_name,
    pr_floor = pr_floor
  )
  
  cat("  Step 2: Interpolating anomaly and reconstructing downscaled fields...\n")
  ds_stack <- delta_reconstruct(
    anom_stack = anom_stack,
    obs_clim_hr = obs_clim_hr,
    dates = dates,
    var_name = var_name
  )
  
  out_dir_hist <- file.path(out_root, "historical", model)
  
  cat("  Step 3: Saving annual monthly stacks...\n")
  save_downscaled_by_year(
    ds_stack = ds_stack,
    dates = dates,
    out_dir = out_dir_hist,
    var_name = var_name,
    model = model,
    scenario = "historical",
    period_label = period_label
  )
  
  rm(r_gcm, dates, anom_stack, ds_stack, gcm_clim)
  gc()
}

# =============================================================================
# 11. DOWNSCALE FUTURE MONTHLY SERIES
# =============================================================================

for (scenario_name in scenarios) {
  
  cat(sprintf("\n=== Downscaling future monthly series: %s ===\n", scenario_name))
  
  fut_files <- find_future_gcm_files(
    gcm_root = gcm_root,
    var_name = var_name,
    scenario_name = scenario_name,
    model_keep = model_keep
  )
  
  if (length(fut_files) == 0) {
    warning("No future files found for scenario: ", scenario_name)
    next
  }
  
  cat("  Number of future files found:", length(fut_files), "\n")
  
  for (f in fut_files) {
    
    model <- get_model_name(f, known_models)
    
    cat(sprintf("\n--- Model: %s | Scenario: %s ---\n", model, scenario_name))
    cat("  File:", basename(f), "\n")
    
    if (is.null(model_clim_list[[model]])) {
      warning("No precomputed historical climatology for model ", model, " — skipping.")
      next
    }
    
    r_gcm <- rast(f)
    dates <- make_monthly_dates(fut_start_ym, nlyr(r_gcm))
    
    if (run_mode == "test") {
      sub <- subset_test_period(
        r = r_gcm,
        dates = dates,
        start_ym = test_fut_start,
        end_ym = test_fut_end
      )
      
      r_gcm <- sub$r
      dates <- sub$dates
      
      cat(
        sprintf(
          "  TEST period: %s to %s (%d layers)\n",
          min(dates),
          max(dates),
          length(dates)
        )
      )
      
      period_label <- "test"
      
    } else {
      cat(
        sprintf(
          "  FULL period: %s to %s (%d layers)\n",
          min(dates),
          max(dates),
          length(dates)
        )
      )
      
      period_label <- fut_period_label
    }
    
    gcm_clim <- model_clim_list[[model]]
    
    cat("  Step 1: Computing anomalies relative to historical climatology...\n")
    anom_stack <- compute_gcm_anomaly(
      r_gcm = r_gcm,
      gcm_clim = gcm_clim,
      dates = dates,
      var_name = var_name,
      pr_floor = pr_floor
    )
    
    cat("  Step 2: Interpolating anomaly and reconstructing downscaled fields...\n")
    ds_stack <- delta_reconstruct(
      anom_stack = anom_stack,
      obs_clim_hr = obs_clim_hr,
      dates = dates,
      var_name = var_name
    )
    
    out_dir_fut <- file.path(out_root, scenario_name, model)
    
    cat("  Step 3: Saving annual monthly stacks...\n")
    save_downscaled_by_year(
      ds_stack = ds_stack,
      dates = dates,
      out_dir = out_dir_fut,
      var_name = var_name,
      model = model,
      scenario = scenario_name,
      period_label = period_label
    )
    
    rm(r_gcm, dates, anom_stack, ds_stack, gcm_clim)
    gc()
  }
}

cat("\n=== Modified delta downscaling complete ===\n")
cat("Outputs saved in:\n")
cat(out_root, "\n")
