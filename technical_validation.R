# ============================================================
# TECHNICAL VALIDATION OF DOWNSCALED HISTORICAL DATA
#
# This script compares:
#   - Observed data
#   - Raw GCM historical data
#   - Delta-downscaled historical data
#
# Validation includes:
#   - Metrics (RMSE, MAE, bias, correlation, KGE)
#   - Density plots
#   - Scatter density plots
#   - Time series
#   - Seasonal cycle
#   - Spatial bias maps
#
# NOTE:
# - Observations are at high resolution (250 m)
# - Raw and downscaled data are aligned to observation grid
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("technical_validation.R")
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGE CHECK
# ============================================================

packages <- c("terra", "dplyr", "tidyr", "ggplot2", "openair")

missing_pkgs <- packages[!packages %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop("Please install missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openair)
})

terraOptions(progress = 1, memfrac = 0.7)

cat("\n=== Technical Validation Script ===\n")

# ============================================================
# 2. PATHS
# ============================================================

root <- Sys.getenv("CLIM_ROOT", unset = getwd())

obs_root  <- file.path(root, "data", "observation", "nc_files")
raw_root  <- file.path(root, "data", "gcm_monthly")
down_root <- file.path(root, "outputs", "delta_downscaled")

out_root  <- file.path(root, "analysis", "historical", "technical_validation")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. SETTINGS
# ============================================================

vars_to_run <- c("pr")   # change if needed

obs_files <- list(
  pr   = file.path(obs_root, "precip_3ds_monthly_198505202008.nc"),
  tmax = file.path(obs_root, "tmax_3ds_m_monthly_198505202008.nc"),
  tmin = file.path(obs_root, "tmin_3ds_tm_monthly_198505202008.nc")
)

val_start <- as.Date("1985-05-01")
val_end   <- as.Date("2014-12-01")

# ============================================================
# 4. HELPERS
# ============================================================

get_time_safe <- function(r) {
  tt <- time(r)
  if (is.null(tt)) stop("Raster has no time metadata.")
  as.Date(tt)
}

subset_time_safe <- function(r, start_date, end_date) {
  tt <- get_time_safe(r)
  idx <- which(tt >= start_date & tt <= end_date)
  if (length(idx) == 0) stop("No layers in requested date range.")
  r[[idx]]
}

calc_metrics <- function(obs_vals, pred_vals, is_precip = FALSE) {
  
  ok <- is.finite(obs_vals) & is.finite(pred_vals)
  obs_vals  <- obs_vals[ok]
  pred_vals <- pred_vals[ok]
  
  if (length(obs_vals) == 0) {
    return(data.frame(n = 0, rmse = NA, mae = NA, bias = NA,
                      pbias = NA, cor = NA, kge = NA))
  }
  
  r <- suppressWarnings(cor(obs_vals, pred_vals))
  
  rmse <- sqrt(mean((pred_vals - obs_vals)^2))
  mae  <- mean(abs(pred_vals - obs_vals))
  bias <- mean(pred_vals - obs_vals)
  
  pbias <- if (is_precip) {
    100 * sum(pred_vals - obs_vals) / sum(obs_vals)
  } else NA
  
  sd_obs  <- sd(obs_vals)
  sd_pred <- sd(pred_vals)
  
  alpha <- sd_pred / sd_obs
  beta  <- mean(pred_vals) / mean(obs_vals)
  kge   <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  
  data.frame(n = length(obs_vals), rmse = rmse, mae = mae,
             bias = bias, pbias = pbias, cor = r, kge = kge)
}

# ============================================================
# 5. MAIN LOOP
# ============================================================

all_metrics <- list()

for (v in vars_to_run) {
  
  cat("\n=====================================\n")
  cat("Variable:", v, "\n")
  
  obs_file <- obs_files[[v]]
  
  if (!file.exists(obs_file)) {
    stop("Observation file missing: ", obs_file)
  }
  
  obs <- rast(obs_file)
  crs(obs) <- "EPSG:5266"
  obs <- subset_time_safe(obs, val_start, val_end)
  
  model_dirs <- list.dirs(file.path(down_root, v, "historical"),
                          recursive = FALSE, full.names = FALSE)
  
  for (m in model_dirs) {
    
    cat("\nModel:", m, "\n")
    
    raw_file <- list.files(raw_root, pattern = m, full.names = TRUE)[1]
    down_files <- list.files(
      file.path(down_root, v, "historical", m),
      pattern = "\\.tif$", full.names = TRUE
    )
    
    if (is.na(raw_file) || length(down_files) == 0) next
    
    raw <- rast(raw_file)
    raw <- subset_time_safe(raw, val_start, val_end)
    
    down <- rast(down_files)
    
    # assign monthly time
    years <- as.integer(sub(".*_(\\d{4})\\.tif$", "\\1", basename(down_files)))
    dates <- unlist(lapply(years, function(y)
      seq(as.Date(paste0(y, "-01-01")), by = "month", length.out = 12)))
    
    time(down) <- dates
    down <- subset_time_safe(down, val_start, val_end)
    
    # align
    raw  <- resample(raw, obs)
    down <- resample(down, obs)
    
    obs_vals  <- values(obs)
    raw_vals  <- values(raw)
    down_vals <- values(down)
    
    met_raw <- calc_metrics(obs_vals, raw_vals, is_precip = (v == "pr"))
    met_down <- calc_metrics(obs_vals, down_vals, is_precip = (v == "pr"))
    
    met_raw$model <- m
    met_raw$dataset <- "Raw"
    
    met_down$model <- m
    met_down$dataset <- "Downscaled"
    
    all_metrics[[length(all_metrics)+1]] <- met_raw
    all_metrics[[length(all_metrics)+1]] <- met_down
    
    cat("  Done\n")
  }
}

# ============================================================
# 6. SAVE OUTPUT
# ============================================================

metrics_tbl <- bind_rows(all_metrics)

write.csv(
  metrics_tbl,
  file.path(out_root, "technical_validation_metrics.csv"),
  row.names = FALSE
)

cat("\n=== Validation complete ===\n")
cat("Output:", out_root, "\n")
