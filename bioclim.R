# ============================================================
# Generate bio1–bio19 from downscaled monthly climatologies
#
# This script generates the 19 standard bioclimatic variables
# from downscaled monthly climatology rasters for precipitation,
# maximum temperature, and minimum temperature. The method and variable naming convention
# follow the standard WorldClim bio1–bio19 definitions.
#
#   - Input rasters must contain 12 monthly layers.
#   - All input rasters must have matching geometry.
#   - Output rasters are saved with LZW compression.
#
# Input filename examples:
#   pr_historical_climatology_ACCESS-CM2.tif
#   pr_ssp126_2021_2040_climatology_ACCESS-CM2.tif
#
# Matching files:
#   tasmax_...
#   tasmin_...
#
# Output structure:
#   Historical:
#     outputs/bioclim/historical/MODEL/bio1.tif
#
#   Future:
#     outputs/bioclim/future/SCENARIO/PERIOD/MODEL/bio1.tif
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("bioclim.R")
# ============================================================

suppressPackageStartupMessages({
  library(terra)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

root <- Sys.getenv(
  "CLIM_ROOT",
  unset = "~/Climate_downscaling/Downscaling"
)

clim_root <- file.path(root, "outputs/climatology/climatology")
out_root  <- file.path(root, "outputs/bioclim")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- file.path(root, "tmp_terra_bioclim")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(
  progress = 1,
  memfrac  = 0.6,
  tempdir  = tmp_dir
)

# ------------------------------------------------------------
# 2. Helper function: extract pixel-wise quarter value
# ------------------------------------------------------------

extract_quarter <- function(stack, index) {
  app(c(stack, index), function(x) {
    idx <- x[length(x)]
    
    if (is.na(idx)) return(NA_real_)
    
    x[idx]
  })
}

# ------------------------------------------------------------
# 3. Function to calculate BIO1–BIO19
# ------------------------------------------------------------

make_bioclim <- function(pr, tmax, tmin) {
  
  names(pr)   <- month.abb
  names(tmax) <- month.abb
  names(tmin) <- month.abb
  
  tmean <- (tmax + tmin) / 2
  
  # Temperature variables
  bio1 <- mean(tmean)                 # Annual mean temperature
  bio2 <- mean(tmax - tmin)           # Mean diurnal range
  bio5 <- max(tmax)                   # Max temperature of warmest month
  bio6 <- min(tmin)                   # Min temperature of coldest month
  bio7 <- bio5 - bio6                 # Temperature annual range
  bio3 <- (bio2 / bio7) * 100         # Isothermality
  bio4 <- app(tmean, sd) * 100        # Temperature seasonality
  
  # Precipitation variables
  bio12 <- sum(pr)                    # Annual precipitation
  bio13 <- max(pr)                    # Precipitation of wettest month
  bio14 <- min(pr)                    # Precipitation of driest month
  bio15 <- (app(pr, sd) / mean(pr)) * 100  # Precipitation seasonality
  
  # Rolling three-month quarters
  q_index <- list(
    c(1, 2, 3),
    c(2, 3, 4),
    c(3, 4, 5),
    c(4, 5, 6),
    c(5, 6, 7),
    c(6, 7, 8),
    c(7, 8, 9),
    c(8, 9, 10),
    c(9, 10, 11),
    c(10, 11, 12),
    c(11, 12, 1),
    c(12, 1, 2)
  )
  
  q_temp <- rast(lapply(q_index, function(i) mean(tmean[[i]])))
  q_prec <- rast(lapply(q_index, function(i) sum(pr[[i]])))
  
  names(q_temp) <- paste0("q", 1:12)
  names(q_prec) <- paste0("q", 1:12)
  
  # Pixel-wise quarter indices
  wet_q  <- app(q_prec, which.max)
  dry_q  <- app(q_prec, which.min)
  warm_q <- app(q_temp, which.max)
  cold_q <- app(q_temp, which.min)
  
  # Quarter-based BIO variables
  bio8  <- extract_quarter(q_temp, wet_q)    # Mean temp of wettest quarter
  bio9  <- extract_quarter(q_temp, dry_q)    # Mean temp of driest quarter
  bio10 <- max(q_temp)                       # Mean temp of warmest quarter
  bio11 <- min(q_temp)                       # Mean temp of coldest quarter
  
  bio16 <- max(q_prec)                       # Precipitation of wettest quarter
  bio17 <- min(q_prec)                       # Precipitation of driest quarter
  bio18 <- extract_quarter(q_prec, warm_q)   # Precipitation of warmest quarter
  bio19 <- extract_quarter(q_prec, cold_q)   # Precipitation of coldest quarter
  
  bio <- c(
    bio1, bio2, bio3, bio4, bio5, bio6, bio7,
    bio8, bio9, bio10, bio11,
    bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19
  )
  
  names(bio) <- paste0("bio", 1:19)
  
  bio
}

# ------------------------------------------------------------
# 4. Find precipitation climatology files
# ------------------------------------------------------------

pr_files <- list.files(
  clim_root,
  pattern = "^pr_.*_climatology_.*\\.tif$",
  full.names = TRUE
)

if (length(pr_files) == 0) {
  stop("No precipitation climatology files found in: ", clim_root)
}

# ------------------------------------------------------------
# 5. Process all climatology files
# ------------------------------------------------------------

for (pr_file in pr_files) {
  
  fname <- basename(pr_file)
  
  message("\n====================================================")
  message("Processing: ", fname)
  message("====================================================")
  
  # Matching tasmax and tasmin files
  tasmax_file <- sub("^pr_", "tasmax_", fname)
  tasmin_file <- sub("^pr_", "tasmin_", fname)
  
  tasmax_path <- file.path(dirname(pr_file), tasmax_file)
  tasmin_path <- file.path(dirname(pr_file), tasmin_file)
  
  if (!file.exists(tasmax_path)) {
    warning("Missing tasmax file: ", tasmax_path)
    next
  }
  
  if (!file.exists(tasmin_path)) {
    warning("Missing tasmin file: ", tasmin_path)
    next
  }
  
  # Read rasters
  pr   <- rast(pr_file)
  tmax <- rast(tasmax_path)
  tmin <- rast(tasmin_path)
  
  # Check monthly layers
  if (nlyr(pr) != 12 || nlyr(tmax) != 12 || nlyr(tmin) != 12) {
    warning("One or more files do not have 12 monthly layers: ", fname)
    next
  }
  
  # Check geometry
  if (!compareGeom(pr, tmax, tmin, stopOnError = FALSE)) {
    warning("Geometry mismatch among pr, tasmax, and tasmin for: ", fname)
    next
  }
  
  # Generate BIO variables
  bio <- make_bioclim(pr, tmax, tmin)
  
  # ----------------------------------------------------------
  # Output folder structure
  # ----------------------------------------------------------
  
  base <- tools::file_path_sans_ext(fname)
  model <- sub(".*_climatology_", "", base)
  
  if (grepl("^pr_historical_climatology_", base)) {
    
    out_dir <- file.path(out_root, "historical", model)
    
  } else {
    
    scenario <- sub("^pr_(ssp[0-9]+)_.*", "\\1", base)
    period   <- sub(
      "^pr_ssp[0-9]+_([0-9]{4})_([0-9]{4})_.*",
      "\\1-\\2",
      base
    )
    
    out_dir <- file.path(out_root, "future", scenario, period, model)
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save BIO1–BIO19 as individual compressed GeoTIFFs
  for (i in seq_len(nlyr(bio))) {
    
    out_path <- file.path(out_dir, paste0("bio", i, ".tif"))
    
    writeRaster(
      bio[[i]],
      out_path,
      overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=LZW"))
    )
  }
  
  message("Saved BIO1–BIO19 to: ", out_dir)
}

message("\nAll BIOCLIM variables generated successfully.")
# To run:
#   source("generate_bioclim.R")
# ============================================================
