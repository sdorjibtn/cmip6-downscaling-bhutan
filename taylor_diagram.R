# ============================================================
# Taylor diagram for historical validation of downscaled data
#
# This script evaluates historical downscaling performance by:
#   1. Reading observed monthly climate data
#   2. Reading raw historical GCM monthly data
#   3. Reading downscaled historical monthly rasters
#   4. Aligning observed, raw GCM, and downscaled data to a common grid
#   5. Extracting paired observed, raw, and downscaled values
#   6. Producing Taylor diagrams comparing raw and downscaled performance
#   7. Saving Taylor input CSV files and fail logs
#
# Notes:
#   - This script is for Taylor diagram generation only.
#   - It complements the hindcast validation and scatter-density scripts.
#   - Observed data are read from the original 250 m monthly NetCDF files.
#   - Raw GCM data are read from monthly historical GCM NetCDF files.
#   - Downscaled data are read from yearly 250 m historical GeoTIFF stacks.
#   - Taylor diagrams are generated for each variable separately.
#   - The default observed reference point from openair is hidden and replaced
#     manually with a black observed point and label.
#   - Raw GCM and downscaled outputs are shown using two colours only.
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("taylor_diagram.R")
# ============================================================
# ============================================================
# 1. USER SETTINGS
# ============================================================

root <- Sys.getenv("CLIM_ROOT")

if (root == "") {
  root <- "~/Climate_downscaling/Downscaling"
}

vars_to_run <- c("pr", "tmax", "tmin")
# vars_to_run <- c("pr")   # use this line for testing

obs_files <- list(
  pr   = file.path(root, "data/observation/nc_files/precip_3ds_monthly_198505202008.nc"),
  tmax = file.path(root, "data/observation/nc_files/tmax_3ds_m_monthly_198505202008.nc"),
  tmin = file.path(root, "data/observation/nc_files/tmin_3ds_tm_monthly_198505202008.nc")
)

raw_root  <- file.path(root, "data/gcm_monthly")
down_root <- file.path(root, "outputs/delta_downscaling_pr_tmax_tmin")

out_root <- file.path(root, "analysis/historical/taylor_diagram")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

val_start <- as.Date("1985-05-01")
val_end   <- as.Date("2014-12-01")

max_points_per_model <- 100000

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

get_time_safe <- function(r) {
  tt <- time(r)
  if (is.null(tt)) stop("Raster has no time metadata.")
  as.Date(tt)
}

subset_time_safe <- function(r, start_date, end_date) {
  tt <- get_time_safe(r)
  idx <- which(tt >= start_date & tt <= end_date)
  if (length(idx) == 0) stop("No layers found within requested date range.")
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

get_model_names <- function(root_dir, var_name) {
  down_var <- get_down_var_name(var_name)
  hist_dir <- file.path(root_dir, down_var, "historical")
  
  if (!dir.exists(hist_dir)) return(character(0))
  
  dirs <- list.dirs(hist_dir, recursive = FALSE, full.names = TRUE)
  basename(dirs)
}

find_raw_file <- function(var_name, model_name, raw_root) {
  raw_token <- get_raw_var_token(var_name)
  
  files <- list.files(
    raw_root,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  files <- files[
    grepl(model_name, files) &
      grepl("historical", files, ignore.case = TRUE) &
      grepl(raw_token, files, ignore.case = TRUE)
  ]
  
  if (length(files) == 0) return(NA_character_)
  files[1]
}

find_down_files <- function(var_name, model_name, down_root) {
  down_var <- get_down_var_name(var_name)
  hist_dir <- file.path(down_root, down_var, "historical", model_name)
  
  if (!dir.exists(hist_dir)) return(character(0))
  
  files <- list.files(hist_dir, pattern = "\\.tif$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  
  years <- sub(".*_(\\d{4})\\.tif$", "\\1", basename(files))
  files[order(as.integer(years))]
}

assign_time_to_yearly_stack <- function(r, files) {
  years <- sub(".*_(\\d{4})\\.tif$", "\\1", basename(files))
  years <- as.integer(years)
  
  dates <- unlist(lapply(years, function(y) {
    seq(as.Date(paste0(y, "-01-01")), by = "1 month", length.out = 12)
  }))
  
  if (nlyr(r) != length(dates)) {
    stop(
      "Layer-date mismatch: raster has ", nlyr(r),
      " layers but expected ", length(dates), "."
    )
  }
  
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

# =============================================================================
# 3. CUSTOM TAYLOR DIAGRAM FUNCTION
# =============================================================================

make_custom_taylor <- function(taylor_df, out_file) {
  
  group_levels <- unique(taylor_df$group)
  
  group_cols <- ifelse(
    grepl("_Raw$", group_levels),
    "#0072B2",   # Raw GCM
    "#D55E00"    # Downscaled
  )
  names(group_cols) <- group_levels
  
  png(out_file, width = 2000, height = 1600, res = 220)
  
  td <- openair::TaylorDiagram(
    mydata  = taylor_df,
    obs     = "obs",
    mod     = "mod",
    group   = "group",
    cols    = group_cols,
    key     = FALSE,
    rms.col = "black",
    
    # Hide default observed point/text
    obs.cex  = 0,
    text.obs = "",
    
    par.settings = list(
      fontsize      = list(fontfamily = "serif", text = 23),
      axis.text     = list(fontfamily = "serif", cex = 0.8),
      par.xlab.text = list(fontfamily = "serif", cex = 1.2, font = 1),
      par.ylab.text = list(fontfamily = "serif", cex = 1.2, font = 1),
      add.text      = list(fontfamily = "serif", font = 1, cex = 1.2)
    ),
    
    pch = 16,
    cex = 2.5
  )
  
  print(td)
  
  # ---------------------------------------------------------------------------
  # Manual edits to lattice Taylor diagram
  # ---------------------------------------------------------------------------
  
  lattice::trellis.focus("panel", 1, 1, highlight = FALSE)
  
  xlim_panel <- grid::current.viewport()$xscale
  ylim_panel <- grid::current.viewport()$yscale
  
  x_min <- xlim_panel[1]
  x_max <- xlim_panel[2]
  y_min <- ylim_panel[1]
  y_max <- ylim_panel[2]
  
  # Remove top and right box lines
  lattice::panel.segments(x_min, y_max, x_max, y_max, col = "white", lwd = 4)
  lattice::panel.segments(x_max, y_min, x_max, y_max, col = "white", lwd = 4)
  
  # Mask "centred RMS error" label
  lattice::panel.rect(
    xleft   = x_max * 0.68,
    ybottom = y_max * 0.84,
    xright  = x_max,
    ytop    = y_max,
    col     = "white",
    border  = "white"
  )
  
  # Add observed point manually in black
  obs_sd <- sd(taylor_df$obs, na.rm = TRUE)
  
  lattice::panel.points(
    x = obs_sd,
    y = 0,
    pch = 16,
    col = "black",
    cex = 1.4
  )
  
  lattice::panel.text(
    x = obs_sd,
    y = y_max * 0.035,
    labels = "observed",
    col = "black",
    cex = 1.1,
    fontfamily = "serif",
    adj = c(0.5, 0)
  )
  
  lattice::trellis.unfocus()
  
  dev.off()
}

# =============================================================================
# 4. MAIN PROCESSING LOOP
# =============================================================================

fail_log <- list()

for (v in vars_to_run) {
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("Processing variable:", v, "\n")
  cat(strrep("=", 70), "\n")
  
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
      reason = "No downscaled model folders found"
    )
    next
  }
  
  taylor_all <- list()
  
  for (m in model_names) {
    
    cat("\nModel:", m, "\n")
    
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
    
    tryCatch({
      
      raw_full <- rast(raw_file)
      raw_full <- subset_time_safe(raw_full, val_start, val_end)
      
      down_full <- rast(down_files)
      down_full <- assign_time_to_yearly_stack(down_full, down_files)
      down_full <- subset_time_safe(down_full, val_start, val_end)
      
      n_common <- min(nlyr(obs_full), nlyr(raw_full), nlyr(down_full))
      if (n_common == 0) stop("No overlapping layers.")
      
      obs  <- obs_full[[1:n_common]]
      raw  <- raw_full[[1:n_common]]
      down <- down_full[[1:n_common]]
      
      aligned <- align_three(obs, raw, down)
      
      obs_vals  <- values(aligned$obs, mat = FALSE)
      raw_vals  <- values(aligned$raw, mat = FALSE)
      down_vals <- values(aligned$down, mat = FALSE)
      
      valid_idx <- which(
        is.finite(obs_vals) &
          is.finite(raw_vals) &
          is.finite(down_vals)
      )
      
      if (length(valid_idx) == 0) {
        stop("No finite paired values after alignment.")
      }
      
      if (length(valid_idx) > max_points_per_model) {
        set.seed(123)
        valid_idx <- sample(valid_idx, max_points_per_model)
      }
      
      taylor_all[[m]] <- data.frame(
        obs = c(obs_vals[valid_idx], obs_vals[valid_idx]),
        mod = c(raw_vals[valid_idx], down_vals[valid_idx]),
        group = c(
          rep(paste0(m, "_Raw"), length(valid_idx)),
          rep(paste0(m, "_Downscaled"), length(valid_idx))
        )
      )
      
      cat("  Added successfully.\n")
      
    }, error = function(e) {
      
      cat("  ERROR:", e$message, "\n")
      
      fail_log[[length(fail_log) + 1]] <<- data.frame(
        variable = v,
        model = m,
        reason = e$message
      )
    })
  }
  
  if (length(taylor_all) == 0) {
    cat("No Taylor data available for", v, "- skipping.\n")
    next
  }
  
  taylor_df <- bind_rows(taylor_all)
  
  out_png <- file.path(
    out_root,
    paste0("taylor_diagram_", v, "_all_models.png")
  )
  
  make_custom_taylor(taylor_df, out_png)
  
  write.csv(
    taylor_df,
    file.path(out_root, paste0("taylor_input_values_", v, ".csv")),
    row.names = FALSE
  )
  
  cat("Taylor diagram saved:", out_png, "\n")
}

# =============================================================================
# 5. SAVE FAIL LOG
# =============================================================================

if (length(fail_log) > 0) {
  fail_tbl <- bind_rows(fail_log)
  
  write.csv(
    fail_tbl,
    file.path(out_root, "taylor_diagram_fail_log.csv"),
    row.names = FALSE
  )
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Finished Taylor diagram generation.\n")
cat("Output folder:", out_root, "\n")
cat(strrep("=", 70), "\n")
