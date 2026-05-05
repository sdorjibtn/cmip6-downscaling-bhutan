# ============================================================
# Scatter-density plots for historical validation
#
# This script compares observed, raw historical GCM, and
# downscaled historical climate data by:
#   1. Reading observed monthly climate data
#   2. Reading raw historical GCM monthly data
#   3. Reading downscaled historical yearly GeoTIFF stacks
#   4. Aligning observed, raw, and downscaled rasters
#   5. Extracting paired values
#   6. Producing scatter-density plots:
#        - observed vs raw GCM
#        - observed vs downscaled
#
# Notes:
#   - This script is for scatter-density plots only.
#   - It complements the metrics and Taylor diagram scripts.
#   - Observed data are read from original 250 m monthly NetCDF files.
#   - Raw GCM data are read from historical monthly NetCDF files.
#   - Downscaled data are read from yearly 250 m historical GeoTIFF stacks.
#   - Large paired datasets are randomly sampled before plotting.
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("scatter_density.R")
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

terraOptions(progress = 1, memfrac = 0.7)

# ============================================================
# 1. SETTINGS
# ============================================================

root <- Sys.getenv("CLIM_ROOT")

if (root == "") {
  root <- "~/Climate_downscaling/Downscaling"
}

vars_to_run <- c("pr", "tmax", "tmin")
# vars_to_run <- c("pr")   # for testing only

obs_files <- list(
  pr   = file.path(root, "data/observation/nc_files/precip_3ds_monthly_198505202008.nc"),
  tmax = file.path(root, "data/observation/nc_files/tmax_3ds_m_monthly_198505202008.nc"),
  tmin = file.path(root, "data/observation/nc_files/tmin_3ds_tm_monthly_198505202008.nc")
)

raw_root  <- file.path(root, "data/gcm_monthly")
down_root <- file.path(root, "outputs/delta_downscaling_pr_tmax_tmin")

out_root <- file.path(root, "analysis/historical/scatter_density")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

val_start <- as.Date("1985-05-01")
val_end   <- as.Date("2014-12-01")

max_points_per_model <- 200000

# ============================================================
# 2. HELPERS
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

get_down_var_name <- function(var_name) {
  switch(
    var_name,
    pr   = "pr",
    tmax = "tasmax",
    tmin = "tasmin",
    stop("Unknown variable: ", var_name)
  )
}

get_raw_var_token <- function(var_name) {
  switch(
    var_name,
    pr   = "pr",
    tmax = "tasmax",
    tmin = "tasmin",
    stop("Unknown variable: ", var_name)
  )
}

get_model_names <- function(root, var_name) {
  
  down_var <- get_down_var_name(var_name)
  hist_dir <- file.path(root, down_var, "historical")
  
  if (!dir.exists(hist_dir)) return(character(0))
  
  dirs <- list.dirs(hist_dir, recursive = FALSE, full.names = TRUE)
  model_names <- basename(dirs)
  model_names <- model_names[model_names != ""]
  
  model_names
}

find_raw_file <- function(var_name, model_name, raw_root) {
  
  raw_token <- get_raw_var_token(var_name)
  
  ff <- list.files(
    raw_root,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  ff <- ff[
    grepl(model_name, ff) &
      grepl("historical", ff, ignore.case = TRUE) &
      grepl(raw_token, ff, ignore.case = TRUE)
  ]
  
  if (length(ff) == 0) return(NA_character_)
  ff[1]
}

find_down_files <- function(var_name, model_name, down_root) {
  
  down_var <- get_down_var_name(var_name)
  hist_dir <- file.path(down_root, down_var, "historical", model_name)
  
  if (!dir.exists(hist_dir)) return(character(0))
  
  ff <- list.files(
    hist_dir,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(ff) == 0) return(character(0))
  
  years <- sub(".*_(\\d{4})\\.tif$", "\\1", basename(ff))
  ff[order(as.integer(years))]
}

assign_time_to_yearly_stack <- function(r, files) {
  
  years <- sub(".*_(\\d{4})\\.tif$", "\\1", basename(files))
  years <- as.integer(years)
  
  expected_layers <- length(years) * 12
  
  if (nlyr(r) != expected_layers) {
    stop(
      "Expected ", expected_layers,
      " layers from yearly TIFFs but found ", nlyr(r)
    )
  }
  
  dates <- unlist(lapply(years, function(y) {
    seq(as.Date(paste0(y, "-01-01")), by = "1 month", length.out = 12)
  }))
  
  time(r) <- dates
  r
}

align_three <- function(obs, raw, down) {
  
  if (!same.crs(obs, raw)) {
    raw <- project(raw, obs)
  }
  
  if (!same.crs(obs, down)) {
    down <- project(down, obs)
  }
  
  raw  <- resample(raw, obs, method = "bilinear")
  down <- resample(down, obs, method = "bilinear")
  
  raw  <- crop(raw, obs)
  down <- crop(down, obs)
  
  list(obs = obs, raw = raw, down = down)
}

# ============================================================
# 3. SCATTER-DENSITY PLOT FUNCTION
# ============================================================

make_scatter_density_plot <- function(df,
                                      xcol,
                                      ycol,
                                      title_txt,
                                      out_png,
                                      var_name = "pr",
                                      y_label_type = c("downscaled", "raw")) {
  
  y_label_type <- match.arg(y_label_type)
  
  if (var_name == "pr") {
    
    x_lab <- expression("Observed")
    
    y_lab <- if (y_label_type == "downscaled") {
      expression("Downscaled")
    } else {
      expression("Raw")
    }
    
    lims <- range(c(df[[xcol]], df[[ycol]]), na.rm = TRUE)
    lims <- c(0, ceiling(max(lims, na.rm = TRUE) / 100) * 100)
    
  } else {
    
    x_lab <- expression("Observed")
    
    y_lab <- if (y_label_type == "downscaled") {
      expression("Downscaled")
    } else {
      expression("Raw")
    }
    
    lims <- range(c(df[[xcol]], df[[ycol]]), na.rm = TRUE)
    lims <- c(
      floor(lims[1] / 10) * 10,
      ceiling(lims[2] / 10) * 10
    )
  }
  
  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    stat_bin2d(bins = 180) +
    geom_abline(
      slope = 1,
      intercept = 0,
      colour = "black",
      linewidth = 0.5
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      trans = "log10",
      guide = "none"
    ) +
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
    theme_bw(base_size = 20, base_family = "serif") +
    theme(
      panel.grid = element_line(colour = "grey90")
    ) +
    labs(
      # title = title_txt,
      x = x_lab,
      y = y_lab
    )
  
  ggsave(out_png, p, width = 5.5, height = 5, dpi = 300)
}

# ============================================================
# 4. MAIN LOOP
# ============================================================

fail_log <- list()

for (v in vars_to_run) {
  
  cat("\n====================================================\n")
  cat("VARIABLE:", v, "\n")
  cat("====================================================\n")
  
  obs_file <- obs_files[[v]]
  
  if (!file.exists(obs_file)) {
    fail_log[[length(fail_log) + 1]] <- data.frame(
      variable = v,
      model = NA_character_,
      reason = "Observed file not found"
    )
    next
  }
  
  obs_full <- rast(obs_file)
  crs(obs_full) <- "EPSG:5266"
  obs_full <- subset_time_safe(obs_full, val_start, val_end)
  
  model_names <- get_model_names(down_root, v)
  
  if (length(model_names) == 0) {
    fail_log[[length(fail_log) + 1]] <- data.frame(
      variable = v,
      model = NA_character_,
      reason = "No model folders found in downscaled directory"
    )
    next
  }
  
  for (m in model_names) {
    
    cat("\nMODEL:", m, "\n")
    
    raw_file   <- find_raw_file(v, m, raw_root)
    down_files <- find_down_files(v, m, down_root)
    
    if (is.na(raw_file) || length(down_files) == 0) {
      fail_log[[length(fail_log) + 1]] <- data.frame(
        variable = v,
        model = m,
        reason = paste(
          ifelse(is.na(raw_file), "raw missing", ""),
          ifelse(length(down_files) == 0, "downscaled missing", "")
        )
      )
      next
    }
    
    cat("  Raw file:", raw_file, "\n")
    cat("  Number of downscaled yearly files:", length(down_files), "\n")
    cat("  First downscaled file:", down_files[1], "\n")
    
    tryCatch({
      
      raw_full <- rast(raw_file)
      raw_full <- subset_time_safe(raw_full, val_start, val_end)
      
      down_full <- rast(down_files)
      down_full <- assign_time_to_yearly_stack(down_full, down_files)
      down_full <- subset_time_safe(down_full, val_start, val_end)
      
      n_common <- min(nlyr(obs_full), nlyr(raw_full), nlyr(down_full))
      
      if (n_common == 0) {
        stop("No overlapping layers after time subset.")
      }
      
      obs  <- obs_full[[1:n_common]]
      raw  <- raw_full[[1:n_common]]
      down <- down_full[[1:n_common]]
      
      aligned <- align_three(obs, raw, down)
      
      obs  <- aligned$obs
      raw  <- aligned$raw
      down <- aligned$down
      
      obs_vals  <- values(obs, mat = FALSE)
      raw_vals  <- values(raw, mat = FALSE)
      down_vals <- values(down, mat = FALSE)
      
      df_plot <- data.frame(
        obs  = obs_vals,
        raw  = raw_vals,
        down = down_vals
      ) |>
        filter(if_all(everything(), is.finite))
      
      if (nrow(df_plot) == 0) {
        stop("No finite values available after alignment.")
      }
      
      if (nrow(df_plot) > max_points_per_model) {
        set.seed(123)
        df_plot <- df_plot[sample(seq_len(nrow(df_plot)), max_points_per_model), ]
      }
      
      plot_dir <- file.path(out_root, "plots", v, m)
      dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
      
      make_scatter_density_plot(
        df = df_plot,
        xcol = "obs",
        ycol = "raw",
        title_txt = paste(v, "-", m, "- Raw GCM"),
        out_png = file.path(plot_dir, paste0(v, "_", m, "_scatter_raw.png")),
        var_name = v,
        y_label_type = "raw"
      )
      
      make_scatter_density_plot(
        df = df_plot,
        xcol = "obs",
        ycol = "down",
        title_txt = paste(v, "-", m, "- Downscaled"),
        out_png = file.path(plot_dir, paste0(v, "_", m, "_scatter_downscaled.png")),
        var_name = v,
        y_label_type = "downscaled"
      )
      
      cat("  Done.\n")
      
    }, error = function(e) {
      
      cat("  ERROR:", e$message, "\n")
      
      fail_log[[length(fail_log) + 1]] <<- data.frame(
        variable = v,
        model = m,
        reason = e$message
      )
    })
  }
}

# ============================================================
# 5. SAVE FAIL LOG
# ============================================================

if (length(fail_log) > 0) {
  
  fail_tbl <- bind_rows(fail_log)
  
  write.csv(
    fail_tbl,
    file.path(out_root, "scatter_density_fail_log.csv"),
    row.names = FALSE
  )
}

cat("\n====================================================\n")
cat("FINISHED SCATTER-DENSITY PLOTS\n")
cat("Output folder:", out_root, "\n")
cat("====================================================\n")
