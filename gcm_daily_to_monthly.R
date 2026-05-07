# Purpose:
# Convert CMIP6 daily GCM NetCDF files to monthly NetCDF files.
#
# Main features:
#   - Handles pr, tasmax, tasmin
#   - Converts 0–360 longitude to -180–180 when needed
#   - Handles standard, noleap, 365_day, and 360_day calendars
#   - Converts precipitation from kg m-2 s-1 to mm/month
#   - Converts temperature from Kelvin to degree Celsius
#   - Corrects latitude orientation
#   - Handles extra singleton dimensions safely
#   - Saves compressed monthly NetCDF files
#
# Example output:
#   pr_mon_ACCESS-CM2_historical_1980-2014.nc
#
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
})

# ============================================================
# 1. SETTINGS
# ============================================================

project_root <- "~/Climate_downscaling/Downscaling"

in_dir  <- file.path(project_root, "data", "gcm")
out_dir <- file.path(project_root, "data", "gcm_monthly")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(progress = 1, memfrac = 0.7)

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

infer_varname <- function(fname) {
  b <- basename(fname)

  if (grepl("^pr_", b)) {
    return("pr")
  } else if (grepl("^tasmax_", b)) {
    return("tasmax")
  } else if (grepl("^tasmin_", b)) {
    return("tasmin")
  } else if (grepl("^tas_", b)) {
    return("tas")
  } else {
    return(NA_character_)
  }
}

get_nc_att_safe <- function(nc, varid, attname) {
  x <- tryCatch(ncatt_get(nc, varid, attname), error = function(e) NULL)
  if (is.null(x) || is.null(x$value)) return(NA)
  x$value
}

normalize_dimname <- function(x) {
  tolower(trimws(x))
}

find_dim_index <- function(dim_names, candidates) {
  idx <- which(dim_names %in% candidates)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

parse_nc_time <- function(tt, tunits, calendar = "standard") {
  if (is.na(tunits) || !nzchar(tunits)) {
    stop("Time units are missing.")
  }

  tunits_low <- tolower(tunits)
  calendar   <- tolower(ifelse(is.na(calendar), "standard", calendar))

  if (!grepl("since", tunits_low, fixed = TRUE)) {
    stop("Unsupported time units: ", tunits)
  }

  unit_part   <- trimws(sub("since.*$", "", tunits_low))
  origin_part <- trimws(sub("^.*since", "", tunits_low))
  origin_part <- gsub("t", " ", origin_part, fixed = TRUE)
  origin_part <- trimws(origin_part)

  origin_date_txt <- sub(" .*$", "", origin_part)

  if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", origin_date_txt)) {
    stop("Could not parse origin date from time units: ", tunits)
  }

  if (!unit_part %in% c("days", "day")) {
    stop("Only 'days since ...' time units are supported. Found: ", tunits)
  }

  origin_date <- as.Date(origin_date_txt)

  if (calendar %in% c("standard", "gregorian", "proleptic_gregorian", "julian")) {
    return(origin_date + tt)
  }

  if (calendar %in% c("noleap", "365_day")) {
    y0 <- as.integer(format(origin_date, "%Y"))
    m0 <- as.integer(format(origin_date, "%m"))
    d0 <- as.integer(format(origin_date, "%d"))

    mdays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

    doy0  <- sum(mdays[seq_len(m0 - 1)]) + d0
    total <- tt + doy0 - 1

    year <- y0 + floor(total / 365)
    doy  <- total %% 365 + 1

    cum_mdays <- cumsum(mdays)
    month <- integer(length(doy))
    day   <- integer(length(doy))

    for (i in seq_along(doy)) {
      month[i] <- which(doy[i] <= cum_mdays)[1]
      prev_cum <- if (month[i] == 1) 0 else cum_mdays[month[i] - 1]
      day[i]   <- doy[i] - prev_cum
    }

    return(as.Date(sprintf("%04d-%02d-%02d", year, month, day)))
  }

  if (calendar == "360_day") {
    y0 <- as.integer(format(origin_date, "%Y"))
    m0 <- as.integer(format(origin_date, "%m"))
    d0 <- as.integer(format(origin_date, "%d"))

    total_days <- tt + (d0 - 1) + (m0 - 1) * 30

    year  <- y0 + floor(total_days / 360)
    rem   <- total_days %% 360
    month <- floor(rem / 30) + 1
    day   <- rem %% 30 + 1

    day <- pmin(day, 28)

    return(as.Date(sprintf("%04d-%02d-%02d", year, month, day)))
  }

  stop("Unsupported calendar: ", calendar)
}

align_metadata_to_array <- function(dim_names_full, dim_lens_full, arr_dims) {
  if (length(dim_names_full) == length(arr_dims)) {
    return(list(dim_names = dim_names_full, dim_lens = dim_lens_full))
  }

  current_names <- dim_names_full
  current_lens  <- dim_lens_full

  while (length(current_names) > length(arr_dims)) {
    idx1 <- which(current_lens == 1)

    if (length(idx1) == 0) {
      stop("Cannot align metadata dimensions to array dimensions.")
    }

    drop_i <- idx1[1]

    current_names <- current_names[-drop_i]
    current_lens  <- current_lens[-drop_i]
  }

  list(dim_names = current_names, dim_lens = current_lens)
}

make_output_filename <- function(r_mon, f, varname, out_dir) {
  b <- basename(f)
  parts <- strsplit(gsub("\\.nc$", "", b), "_")[[1]]

  if (length(parts) < 5) {
    stop("Unexpected filename format: ", b)
  }

  model      <- parts[3]
  experiment <- parts[4]

  dates <- time(r_mon)

  start_year <- format(min(dates), "%Y")
  end_year   <- format(max(dates), "%Y")

  fname <- paste0(
    varname, "_mon_",
    model, "_",
    experiment, "_",
    start_year, "-", end_year, ".nc"
  )

  file.path(out_dir, fname)
}

# ============================================================
# 3. READ DAILY CMIP6 NETCDF
# ============================================================

read_cmip_daily_nc <- function(f) {
  nc <- nc_open(f)
  on.exit(nc_close(nc), add = TRUE)

  varname <- infer_varname(basename(f))

  if (is.na(varname)) {
    stop("Could not infer variable name from file: ", basename(f))
  }

  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  tt  <- ncvar_get(nc, "time")

  # ------------------------------------------------------------
  # Detect and standardise longitude convention
  # ------------------------------------------------------------
  if (max(lon, na.rm = TRUE) > 180) {
    cat("Detected 0–360 longitude. Converting to -180–180.\n")

    lon <- ifelse(lon > 180, lon - 360, lon)
    lon_order <- order(lon)
    lon <- lon[lon_order]

  } else {
    cat("Longitude already in -180–180 convention.\n")
    lon_order <- seq_along(lon)
  }

  tunits   <- get_nc_att_safe(nc, "time", "units")
  calendar <- get_nc_att_safe(nc, "time", "calendar")
  v_units  <- get_nc_att_safe(nc, varname, "units")
  fill1    <- get_nc_att_safe(nc, varname, "_FillValue")
  fill2    <- get_nc_att_safe(nc, varname, "missing_value")

  vobj <- nc$var[[varname]]

  if (is.null(vobj)) {
    stop("Variable '", varname, "' not found in file: ", basename(f))
  }

  dim_names_full <- vapply(vobj$dim, function(d) normalize_dimname(d$name), character(1))
  dim_lens_full  <- vapply(vobj$dim, function(d) d$len, numeric(1))

  cat("\n----------------------------------------------------\n")
  cat("Reading:", basename(f), "\n")
  cat("Variable:", varname, "\n")
  cat("Units:", v_units, "\n")
  cat("Calendar:", calendar, "\n")
  cat("Time units:", tunits, "\n")

  arr <- tryCatch(
    ncvar_get(nc, varname, collapse_degen = FALSE),
    error = function(e) ncvar_get(nc, varname)
  )

  if (!is.na(fill1)) arr[arr == fill1] <- NA
  if (!is.na(fill2)) arr[arr == fill2] <- NA

  aligned <- align_metadata_to_array(dim_names_full, dim_lens_full, dim(arr))

  dim_names <- aligned$dim_names

  lon_idx  <- find_dim_index(dim_names, c("lon", "longitude", "x"))
  lat_idx  <- find_dim_index(dim_names, c("lat", "latitude", "y"))
  time_idx <- find_dim_index(dim_names, c("time"))

  if (any(is.na(c(lon_idx, lat_idx, time_idx)))) {
    stop("Could not identify lon, lat, and time dimensions.")
  }

  extra_idx <- setdiff(seq_along(dim_names), c(lon_idx, lat_idx, time_idx))

  if (length(extra_idx) > 0) {
    idx_list <- lapply(seq_along(dim_names), function(i) {
      if (i %in% extra_idx) 1 else seq_len(dim(arr)[i])
    })

    arr <- do.call(`[`, c(list(arr), idx_list, list(drop = TRUE)))
    dim_names <- dim_names[-extra_idx]

    lon_idx  <- find_dim_index(dim_names, c("lon", "longitude", "x"))
    lat_idx  <- find_dim_index(dim_names, c("lat", "latitude", "y"))
    time_idx <- find_dim_index(dim_names, c("time"))
  }

  if (length(dim(arr)) != 3) {
    stop("Expected 3D array after dropping extra dimensions.")
  }

  # Reorder to lat x lon x time
  arr <- aperm(arr, c(lat_idx, lon_idx, time_idx))

  # Apply longitude reorder if 0–360 was converted
  arr <- arr[, lon_order, , drop = FALSE]

  dx <- abs(lon[2] - lon[1])
  dy <- abs(lat[2] - lat[1])

  r <- rast(arr)

  ext(r) <- c(
    min(lon) - dx / 2,
    max(lon) + dx / 2,
    min(lat) - dy / 2,
    max(lat) + dy / 2
  )

  crs(r) <- "EPSG:4326"

  if (lat[1] < lat[length(lat)]) {
    r <- flip(r, direction = "vertical")
  }

  dates <- parse_nc_time(tt, tunits, calendar)

  if (length(dates) != nlyr(r)) {
    stop("Time length does not match raster layers.")
  }

  time(r) <- dates

  list(
    rast     = r,
    varname  = varname,
    units_in = v_units,
    calendar = calendar,
    file     = f
  )
}

# ============================================================
# 4. CONVERT UNITS AND AGGREGATE TO MONTHLY
# ============================================================

convert_units_and_aggregate <- function(r, varname, units_in) {
  dates <- time(r)
  ym    <- format(dates, "%Y-%m")

  if (varname == "pr") {
    units_low <- tolower(trimws(units_in))

    if (
      grepl("kg", units_low) &&
      grepl("m-2", units_low) &&
      grepl("s-1", units_low)
    ) {
      cat("Converting precipitation from kg m-2 s-1 to mm/day.\n")
      r <- r * 86400

    } else if (
      grepl("mm/day", units_low, fixed = TRUE) ||
      grepl("mm d-1", units_low, fixed = TRUE) ||
      grepl("mm day-1", units_low, fixed = TRUE)
    ) {
      cat("Precipitation already in mm/day.\n")

    } else {
      stop("Unknown precipitation units: ", units_in)
    }

    r_mon <- tapp(r, index = ym, fun = function(x) sum(x, na.rm = TRUE))
    r_mon[r_mon < 0] <- 0

  } else if (varname %in% c("tas", "tasmax", "tasmin")) {
    units_low <- tolower(trimws(units_in))

    if (units_low %in% c("k", "kelvin")) {
      cat("Converting temperature from Kelvin to degree Celsius.\n")
      r <- r - 273.15

    } else if (units_low %in% c(
      "degc", "c", "celsius",
      "degrees_celsius", "degree_celsius"
    )) {
      cat("Temperature already in degree Celsius.\n")

    } else {
      stop("Unknown temperature units: ", units_in)
    }

    r_mon <- tapp(r, index = ym, fun = function(x) mean(x, na.rm = TRUE))

  } else {
    stop("Unsupported variable: ", varname)
  }

  mon_dates <- as.Date(paste0(sort(unique(ym)), "-01"))

  if (nlyr(r_mon) != length(mon_dates)) {
    stop("Monthly layer count does not match expected month sequence.")
  }

  time(r_mon) <- mon_dates
  names(r_mon) <- paste0(varname, "_", format(mon_dates, "%Y%m"))

  r_mon
}

# ============================================================
# 5. PROCESS ONE FILE
# ============================================================

aggregate_one_file <- function(f, out_dir) {
  cat("\n====================================================\n")
  cat("Processing:", basename(f), "\n")

  x <- read_cmip_daily_nc(f)

  r        <- x$rast
  varname  <- x$varname
  units_in <- x$units_in

  dates <- time(r)

  cat("Daily layers:", nlyr(r), "\n")
  cat("Date range:", as.character(min(dates)), "to", as.character(max(dates)), "\n")

  r_mon <- convert_units_and_aggregate(r, varname, units_in)

  cat("Monthly layers:", nlyr(r_mon), "\n")

  out_file <- make_output_filename(r_mon, f, varname, out_dir)

  cat("Writing:", out_file, "\n")

  writeCDF(
    r_mon,
    filename    = out_file,
    varname     = varname,
    overwrite   = TRUE,
    compression = 4
  )

  cat("Finished:", basename(out_file), "\n")

  invisible(out_file)
}

# ============================================================
# 6. RUN PREPROCESSING
# ============================================================

files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)

cat("\n====================================================\n")
cat("CMIP6 daily-to-monthly preprocessing\n")
cat("Input directory :", in_dir, "\n")
cat("Output directory:", out_dir, "\n")
cat("Files found     :", length(files), "\n")
cat("====================================================\n")

if (length(files) == 0) {
  stop("No NetCDF files found in: ", in_dir)
}

failed_files <- character(0)

for (f in files) {
  tryCatch(
    aggregate_one_file(f, out_dir),
    error = function(e) {
      failed_files <<- c(failed_files, basename(f))

      cat("\nERROR:", basename(f), "\n")
      cat("Message:", conditionMessage(e), "\n")
    }
  )
}

cat("\n====================================================\n")
cat("Preprocessing completed.\n")
cat("Successful output directory:", out_dir, "\n")

if (length(failed_files) > 0) {
  cat("\nFiles that failed:\n")
  print(failed_files)
} else {
  cat("All files processed successfully.\n")
}

cat("====================================================\n")
