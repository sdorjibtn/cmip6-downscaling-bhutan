# ============================================================
# Hindcast validation for modified delta downscaling
#
# This script evaluates historical downscaling performance by:
#   1. Splitting historical monthly data into calibration and validation periods
#   2. Applying monthly delta correction using calibration-period data
#   3. Comparing observed, raw GCM, and downscaled validation series
#   4. Computing correlation, RMSE, MAE, bias, and percent bias
#   5. Saving validation rasters, CSV summaries, and diagnostic plots
#
# Notes:
#   - Observations must already be regridded to each GCM grid.
#   - GCM files must already be aligned to the same grid and time axis.
#   - Metrics are computed using spatially averaged domain-mean monthly time series.
#   - This is a hindcast validation script, not the main 250 m downscaling script.
#
# Expected folder structure:
#
# CLIM_ROOT/
# ├── data/
# │   ├── gcm_monthly_aligned/
# │   ├── observation_regridded_to_each_gcm_pr/
# │   ├── observation_regridded_to_each_gcm_tmax/
# │   └── observation_regridded_to_each_gcm_tmin/
# └── analysis/
#     └── historical/
#         └── hindcast_validation/
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("hindcast_validation.R")
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGE CHECK
# ============================================================

packages <- c("terra", "dplyr", "tidyr", "ggplot2", "readr")

missing_pkgs <- packages[!packages %in% rownames(installed.packages())]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install missing package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

terraOptions(progress = 1, memfrac = 0.7)

cat("\n=== Hindcast Validation for Modified Delta Downscaling ===\n")

# ============================================================
# 2. MAIN SETTINGS
# ============================================================

root <- Sys.getenv("CLIM_ROOT", unset = getwd())

vars_to_run <- c("pr", "tasmax", "tasmin")

gcm_aligned_dir <- file.path(root, "data", "gcm_monthly_aligned")

out_base <- file.path(root, "analysis", "historical", "hindcast_validation")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# Calibration period: used to estimate monthly delta correction
cal_start <- as.Date("1986-01-01")
cal_end   <- as.Date("2005-12-31")

# Validation period: independent period not used for calibration
val_start <- as.Date("2006-01-01")
val_end   <- as.Date("2014-12-31")

# ============================================================
# 3. SAFETY CHECKS
# ============================================================

if (!dir.exists(gcm_aligned_dir)) {
  stop("Aligned GCM folder not found: ", gcm_aligned_dir)
}

if (!dir.exists(out_base)) {
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
}

if (!all(vars_to_run %in% c("pr", "tasmax", "tasmin"))) {
  stop("vars_to_run must contain only: pr, tasmax, tasmin")
}

cat("Project root:", root, "\n")
cat("Aligned GCM directory:", gcm_aligned_dir, "\n")
cat("Output directory:", out_base, "\n")

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

get_settings <- function(var_name, root) {
  
  obs_root <- switch(
    var_name,
    "pr"     = file.path(root, "data", "observation_regridded_to_each_gcm_pr"),
    "tasmax" = file.path(root, "data", "observation_regridded_to_each_gcm_tmax"),
    "tasmin" = file.path(root, "data", "observation_regridded_to_each_gcm_tmin"),
    stop("Unsupported variable: ", var_name)
  )
  
  obs_varname <- switch(
    var_name,
    "pr"     = "precip",
    "tasmax" = "tasmax",
    "tasmin" = "tasmin"
  )
  
  y_label <- switch(
    var_name,
    "pr"     = "Precipitation",
    "tasmax" = "Maximum temperature",
    "tasmin" = "Minimum temperature"
  )
  
  list(
    obs_root = obs_root,
    obs_varname = obs_varname,
    y_label = y_label
  )
}

find_gcm_file <- function(model_name, var_name, gcm_aligned_dir) {
  
  patt <- paste0("^", var_name, "_mon_", model_name, "_historical_.*\\.nc$")
  
  hits <- list.files(
    gcm_aligned_dir,
    pattern = patt,
    full.names = TRUE
  )
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  hits[1]
}

calc_metrics <- function(obs_vals, sim_vals) {
  
  ok <- is.finite(obs_vals) & is.finite(sim_vals)
  obs_vals <- obs_vals[ok]
  sim_vals <- sim_vals[ok]
  
  if (length(obs_vals) < 2) {
    return(data.frame(
      n = length(obs_vals),
      cor = NA_real_,
      rmse = NA_real_,
      mae = NA_real_,
      bias = NA_real_,
      pbias = NA_real_
    ))
  }
  
  data.frame(
    n     = length(obs_vals),
    cor   = suppressWarnings(cor(obs_vals, sim_vals)),
    rmse  = sqrt(mean((sim_vals - obs_vals)^2)),
    mae   = mean(abs(sim_vals - obs_vals)),
    bias  = mean(sim_vals - obs_vals),
    pbias = ifelse(
      sum(obs_vals, na.rm = TRUE) == 0,
      NA_real_,
      100 * sum(sim_vals - obs_vals, na.rm = TRUE) / sum(obs_vals, na.rm = TRUE)
    )
  )
}

monthly_delta_correct <- function(sim_cal, obs_cal, sim_val, var_name = "pr") {
  
  out <- rep(NA_real_, length(sim_val))
  
  ok_cal <- is.finite(sim_cal) & is.finite(obs_cal)
  ok_val <- is.finite(sim_val)
  
  if (sum(ok_cal) < 3 || sum(ok_val) == 0) {
    return(out)
  }
  
  sim_cal_mean <- mean(sim_cal[ok_cal], na.rm = TRUE)
  obs_cal_mean <- mean(obs_cal[ok_cal], na.rm = TRUE)
  
  if (!is.finite(sim_cal_mean) || !is.finite(obs_cal_mean)) {
    return(out)
  }
  
  if (var_name == "pr") {
    
    if (sim_cal_mean <= 0) {
      return(out)
    }
    
    ratio <- sim_val[ok_val] / sim_cal_mean
    corrected <- obs_cal_mean * ratio
    corrected[corrected < 0] <- 0
    
  } else {
    
    delta <- sim_val[ok_val] - sim_cal_mean
    corrected <- obs_cal_mean + delta
  }
  
  out[ok_val] <- corrected
  
  out
}

get_season <- function(month_num) {
  
  dplyr::case_when(
    month_num %in% c(12, 1, 2)  ~ "DJF",
    month_num %in% c(3, 4, 5)   ~ "MAM",
    month_num %in% c(6, 7, 8)   ~ "JJA",
    month_num %in% c(9, 10, 11) ~ "SON"
  )
}

# ============================================================
# 5. SUMMARY PLOTS
# ============================================================

make_summary_plots <- function(summary_df, all_ts_df, var_name, out_root, y_label) {
  
  out_dir <- file.path(out_root, "summary_plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  summary_df$dataset <- factor(summary_df$dataset, levels = c("raw", "downscaled"))
  
  wide <- summary_df %>%
    select(model, dataset, cor, rmse, mae, bias, pbias) %>%
    pivot_wider(
      names_from = dataset,
      values_from = c(cor, rmse, mae, bias, pbias)
    )
  
  write_csv(
    wide,
    file.path(out_dir, paste0("hindcast_raw_vs_downscaled_wide_", var_name, ".csv"))
  )
  
  improve_tbl <- wide %>%
    mutate(
      d_cor   = cor_downscaled - cor_raw,
      d_rmse  = rmse_downscaled - rmse_raw,
      d_mae   = mae_downscaled - mae_raw,
      d_bias  = abs(bias_downscaled) - abs(bias_raw),
      d_pbias = abs(pbias_downscaled) - abs(pbias_raw)
    )
  
  write_csv(
    improve_tbl,
    file.path(out_dir, paste0("hindcast_improvement_table_", var_name, ".csv"))
  )
  
  rank_tbl <- improve_tbl %>%
    mutate(
      rank_cor   = rank(-cor_downscaled, ties.method = "average", na.last = "keep"),
      rank_rmse  = rank(rmse_downscaled, ties.method = "average", na.last = "keep"),
      rank_mae   = rank(mae_downscaled, ties.method = "average", na.last = "keep"),
      rank_bias  = rank(abs(bias_downscaled), ties.method = "average", na.last = "keep"),
      rank_pbias = rank(abs(pbias_downscaled), ties.method = "average", na.last = "keep")
    ) %>%
    mutate(
      mean_rank = rowMeans(
        select(., rank_cor, rank_rmse, rank_mae, rank_bias, rank_pbias),
        na.rm = TRUE
      )
    ) %>%
    arrange(mean_rank)
  
  write_csv(
    rank_tbl,
    file.path(out_dir, paste0("hindcast_model_ranking_", var_name, ".csv"))
  )
  
  make_dumbbell <- function(data, raw_col, down_col, xlab, title, out_file) {
    
    plot_df <- data %>%
      select(model, raw = all_of(raw_col), downscaled = all_of(down_col)) %>%
      mutate(model = reorder(model, downscaled))
    
    p <- ggplot(plot_df, aes(y = model)) +
      geom_segment(
        aes(x = raw, xend = downscaled, yend = model),
        linewidth = 0.7
      ) +
      geom_point(aes(x = raw), size = 2) +
      geom_point(aes(x = downscaled), size = 2) +
      labs(title = title, x = xlab, y = "Model") +
      theme_bw()
    
    ggsave(out_file, p, width = 8, height = 6, dpi = 300)
  }
  
  make_dumbbell(
    wide, "cor_raw", "cor_downscaled",
    "Correlation",
    paste("Hindcast correlation:", var_name),
    file.path(out_dir, "dumbbell_correlation.png")
  )
  
  make_dumbbell(
    wide, "rmse_raw", "rmse_downscaled",
    "RMSE",
    paste("Hindcast RMSE:", var_name),
    file.path(out_dir, "dumbbell_rmse.png")
  )
  
  make_dumbbell(
    wide, "mae_raw", "mae_downscaled",
    "MAE",
    paste("Hindcast MAE:", var_name),
    file.path(out_dir, "dumbbell_mae.png")
  )
  
  make_dumbbell(
    wide, "bias_raw", "bias_downscaled",
    "Bias",
    paste("Hindcast bias:", var_name),
    file.path(out_dir, "dumbbell_bias.png")
  )
  
  make_dumbbell(
    wide, "pbias_raw", "pbias_downscaled",
    "Percent bias (%)",
    paste("Hindcast percent bias:", var_name),
    file.path(out_dir, "dumbbell_pbias.png")
  )
  
  heat_df <- summary_df %>%
    select(model, dataset, cor, rmse, mae, bias, pbias) %>%
    pivot_longer(
      cols = c(cor, rmse, mae, bias, pbias),
      names_to = "metric",
      values_to = "value"
    )
  
  p_heat <- ggplot(heat_df, aes(x = metric, y = model, fill = value)) +
    geom_tile() +
    facet_wrap(~dataset) +
    labs(
      title = paste("Hindcast metrics heatmap:", var_name),
      x = "Metric",
      y = "Model"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_dir, "heatmap_metrics.png"),
    p_heat,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  p_rank <- ggplot(rank_tbl, aes(x = reorder(model, mean_rank), y = mean_rank)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("Overall hindcast ranking:", var_name),
      x = "Model",
      y = "Mean rank; lower is better"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_dir, "hindcast_model_ranking.png"),
    p_rank,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ensemble_ts <- all_ts_df %>%
    group_by(date) %>%
    summarise(
      observed = mean(observed, na.rm = TRUE),
      raw_mean = mean(raw, na.rm = TRUE),
      raw_sd = sd(raw, na.rm = TRUE),
      downscaled_mean = mean(downscaled, na.rm = TRUE),
      downscaled_sd = sd(downscaled, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    ensemble_ts,
    file.path(out_dir, paste0("ensemble_hindcast_timeseries_", var_name, ".csv"))
  )
  
  p_ens <- ggplot(ensemble_ts, aes(x = date)) +
    geom_line(aes(y = observed), linewidth = 0.8) +
    geom_ribbon(
      aes(ymin = raw_mean - raw_sd, ymax = raw_mean + raw_sd),
      alpha = 0.20
    ) +
    geom_line(aes(y = raw_mean), linewidth = 0.7, linetype = "dashed") +
    geom_ribbon(
      aes(
        ymin = downscaled_mean - downscaled_sd,
        ymax = downscaled_mean + downscaled_sd
      ),
      alpha = 0.20
    ) +
    geom_line(aes(y = downscaled_mean), linewidth = 0.7) +
    labs(
      title = paste("Ensemble hindcast validation:", var_name),
      x = NULL,
      y = y_label
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_dir, "ensemble_hindcast_timeseries.png"),
    p_ens,
    width = 9,
    height = 4.5,
    dpi = 300
  )
  
  ensemble_clim <- all_ts_df %>%
    group_by(month) %>%
    summarise(
      observed = mean(observed, na.rm = TRUE),
      raw_mean = mean(raw, na.rm = TRUE),
      raw_sd = sd(raw, na.rm = TRUE),
      downscaled_mean = mean(downscaled, na.rm = TRUE),
      downscaled_sd = sd(downscaled, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    ensemble_clim,
    file.path(out_dir, paste0("ensemble_hindcast_monthly_climatology_", var_name, ".csv"))
  )
  
  p_clim <- ggplot(ensemble_clim, aes(x = month)) +
    geom_line(aes(y = observed), linewidth = 0.9) +
    geom_ribbon(
      aes(ymin = raw_mean - raw_sd, ymax = raw_mean + raw_sd),
      alpha = 0.20
    ) +
    geom_line(aes(y = raw_mean), linewidth = 0.8, linetype = "dashed") +
    geom_ribbon(
      aes(
        ymin = downscaled_mean - downscaled_sd,
        ymax = downscaled_mean + downscaled_sd
      ),
      alpha = 0.20
    ) +
    geom_line(aes(y = downscaled_mean), linewidth = 0.8) +
    scale_x_continuous(breaks = 1:12) +
    labs(
      title = paste("Ensemble monthly climatology:", var_name),
      x = "Month",
      y = y_label
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_dir, "ensemble_monthly_climatology.png"),
    p_clim,
    width = 7,
    height = 4.5,
    dpi = 300
  )
  
  ensemble_scatter_df <- all_ts_df %>%
    group_by(date) %>%
    summarise(
      observed = mean(observed, na.rm = TRUE),
      raw = mean(raw, na.rm = TRUE),
      downscaled = mean(downscaled, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    ensemble_scatter_df,
    file.path(out_dir, paste0("ensemble_scatter_data_", var_name, ".csv"))
  )
  
  calc_stats <- function(obs, sim) {
    
    ok <- is.finite(obs) & is.finite(sim)
    obs <- obs[ok]
    sim <- sim[ok]
    
    data.frame(
      r    = suppressWarnings(cor(obs, sim)),
      rmse = sqrt(mean((sim - obs)^2)),
      bias = mean(sim - obs)
    )
  }
  
  stats_raw <- calc_stats(
    ensemble_scatter_df$observed,
    ensemble_scatter_df$raw
  )
  
  stats_down <- calc_stats(
    ensemble_scatter_df$observed,
    ensemble_scatter_df$downscaled
  )
  
  label_raw <- paste0(
    "Raw\n",
    "R = ", round(stats_raw$r, 2), "\n",
    "RMSE = ", round(stats_raw$rmse, 2), "\n",
    "Bias = ", round(stats_raw$bias, 2)
  )
  
  label_down <- paste0(
    "Downscaled\n",
    "R = ", round(stats_down$r, 2), "\n",
    "RMSE = ", round(stats_down$rmse, 2), "\n",
    "Bias = ", round(stats_down$bias, 2)
  )
  
  ensemble_scatter_long <- ensemble_scatter_df %>%
    pivot_longer(
      cols = c(raw, downscaled),
      names_to = "dataset",
      values_to = "simulated"
    )
  
  x_rng <- range(ensemble_scatter_long$observed, na.rm = TRUE)
  y_rng <- range(ensemble_scatter_long$simulated, na.rm = TRUE)
  
  x_pos <- x_rng[1] + 0.97 * diff(x_rng)
  y_down <- y_rng[1] + 0.70 * diff(y_rng)
  y_raw  <- y_rng[1] + 0.30 * diff(y_rng)
  
  p_ensemble_scatter <- ggplot(
    ensemble_scatter_long,
    aes(x = observed, y = simulated, color = dataset)
  ) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    annotate(
      "text",
      x = x_pos,
      y = y_down,
      label = label_down,
      hjust = 1,
      vjust = 1,
      size = 4
    ) +
    annotate(
      "text",
      x = x_pos,
      y = y_raw,
      label = label_raw,
      hjust = 1,
      vjust = 1,
      size = 4
    ) +
    labs(
      title = paste("Ensemble observed vs simulated:", var_name),
      x = paste("Observed", y_label),
      y = paste("Simulated", y_label),
      color = "Dataset"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_dir, paste0("ensemble_scatter_observed_vs_simulated_", var_name, ".png")),
    p_ensemble_scatter,
    width = 6.5,
    height = 5,
    dpi = 300
  )
}

# ============================================================
# 6. RUN ONE VARIABLE
# ============================================================

run_one_variable <- function(var_name) {
  
  cat("\n====================================================\n")
  cat("RUNNING VARIABLE:", var_name, "\n")
  cat("====================================================\n")
  
  st <- get_settings(var_name, root)
  
  obs_root <- st$obs_root
  obs_varname <- st$obs_varname
  y_label <- st$y_label
  
  if (!dir.exists(obs_root)) {
    warning("Observation directory not found for ", var_name, ": ", obs_root)
    return(NULL)
  }
  
  out_root <- file.path(out_base, var_name)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  model_names <- list.dirs(obs_root, recursive = FALSE, full.names = FALSE)
  model_names <- sort(model_names[nzchar(model_names)])
  model_names <- setdiff(model_names, c("shared_250m", "_summary_plots", "summary_plots"))
  
  if (length(model_names) == 0) {
    warning("No model folders found in: ", obs_root)
    return(NULL)
  }
  
  cat("Models found:\n")
  print(model_names)
  
  summary_tbl <- list()
  all_ts_tbl <- list()
  
  for (model_name in model_names) {
    
    cat("\n----------------------------------------------------\n")
    cat("PROCESSING MODEL:", model_name, "\n")
    
    obs_file <- file.path(
      obs_root,
      model_name,
      paste0(obs_varname, "_obs_monthly_on_", model_name, "_grid.tif")
    )
    
    gcm_file <- find_gcm_file(model_name, var_name, gcm_aligned_dir)
    
    if (!file.exists(obs_file)) {
      cat("Skipping: observation file not found:\n", obs_file, "\n")
      next
    }
    
    if (is.na(gcm_file) || !file.exists(gcm_file)) {
      cat("Skipping: aligned GCM file not found for model:", model_name, "\n")
      next
    }
    
    model_out <- file.path(out_root, model_name)
    dir.create(model_out, recursive = TRUE, showWarnings = FALSE)
    
    obs_on_gcm  <- rast(obs_file)
    gcm_aligned <- rast(gcm_file)
    
    if (is.null(time(obs_on_gcm))) {
      cat("Skipping: observation raster has no time metadata:", obs_file, "\n")
      next
    }
    
    if (is.null(time(gcm_aligned))) {
      cat("Skipping: GCM raster has no time metadata:", gcm_file, "\n")
      next
    }
    
    names(obs_on_gcm) <- paste0(var_name, "_obs_", seq_len(nlyr(obs_on_gcm)))
    names(gcm_aligned) <- paste0(var_name, "_gcm_", seq_len(nlyr(gcm_aligned)))
    
    if (!compareGeom(obs_on_gcm[[1]], gcm_aligned[[1]], stopOnError = FALSE)) {
      cat("Skipping: geometry mismatch for model:", model_name, "\n")
      next
    }
    
    obs_cal <- obs_on_gcm[[time(obs_on_gcm) >= cal_start & time(obs_on_gcm) <= cal_end]]
    gcm_cal <- gcm_aligned[[time(gcm_aligned) >= cal_start & time(gcm_aligned) <= cal_end]]
    
    obs_val <- obs_on_gcm[[time(obs_on_gcm) >= val_start & time(obs_on_gcm) <= val_end]]
    gcm_val <- gcm_aligned[[time(gcm_aligned) >= val_start & time(gcm_aligned) <= val_end]]
    
    if (nlyr(obs_cal) == 0 || nlyr(gcm_cal) == 0 || nlyr(obs_val) == 0 || nlyr(gcm_val) == 0) {
      cat("Skipping: empty calibration or validation period for model:", model_name, "\n")
      next
    }
    
    if (nlyr(obs_cal) != nlyr(gcm_cal) || nlyr(obs_val) != nlyr(gcm_val)) {
      cat("Skipping: layer count mismatch for model:", model_name, "\n")
      next
    }
    
    obs_cal_ym <- format(time(obs_cal), "%Y-%m")
    gcm_cal_ym <- format(time(gcm_cal), "%Y-%m")
    obs_val_ym <- format(time(obs_val), "%Y-%m")
    gcm_val_ym <- format(time(gcm_val), "%Y-%m")
    
    if (!all(obs_cal_ym == gcm_cal_ym)) {
      cat("Skipping: calibration year-month mismatch for model:", model_name, "\n")
      next
    }
    
    if (!all(obs_val_ym == gcm_val_ym)) {
      cat("Skipping: validation year-month mismatch for model:", model_name, "\n")
      next
    }
    
    dates_val <- as.Date(paste0(obs_val_ym, "-01"))
    
    mon_cal <- format(time(gcm_cal), "%m")
    mon_val <- format(time(gcm_val), "%m")
    
    gcm_cal_mat <- values(gcm_cal)
    obs_cal_mat <- values(obs_cal)
    gcm_val_mat <- values(gcm_val)
    
    n_cells <- nrow(gcm_val_mat)
    n_time  <- ncol(gcm_val_mat)
    
    hindcast_mat <- matrix(NA_real_, nrow = n_cells, ncol = n_time)
    
    for (m in sprintf("%02d", 1:12)) {
      
      idx_cal <- which(mon_cal == m)
      idx_val <- which(mon_val == m)
      
      cat("  Month:", m, "\n")
      
      if (length(idx_cal) == 0 || length(idx_val) == 0) {
        next
      }
      
      for (i in seq_len(n_cells)) {
        hindcast_mat[i, idx_val] <- monthly_delta_correct(
          sim_cal = gcm_cal_mat[i, idx_cal],
          obs_cal = obs_cal_mat[i, idx_cal],
          sim_val = gcm_val_mat[i, idx_val],
          var_name = var_name
        )
      }
    }
    
    hindcast_downscaled <- rast(gcm_val)
    values(hindcast_downscaled) <- hindcast_mat
    time(hindcast_downscaled) <- time(gcm_val)
    
    out_raster <- file.path(
      model_out,
      paste0("hindcast_monthly_delta_downscaled_", model_name, "_", var_name, ".tif")
    )
    
    writeRaster(hindcast_downscaled, out_raster, overwrite = TRUE)
    
    obs_val_mean  <- global(obs_val, mean, na.rm = TRUE)[, 1]
    raw_val_mean  <- global(gcm_val, mean, na.rm = TRUE)[, 1]
    down_val_mean <- global(hindcast_downscaled, mean, na.rm = TRUE)[, 1]
    
    metrics_raw  <- calc_metrics(obs_val_mean, raw_val_mean)
    metrics_down <- calc_metrics(obs_val_mean, down_val_mean)
    
    metrics_raw$model <- model_name
    metrics_raw$variable <- var_name
    metrics_raw$dataset <- "raw"
    
    metrics_down$model <- model_name
    metrics_down$variable <- var_name
    metrics_down$dataset <- "downscaled"
    
    summary_tbl[[paste0(model_name, "_raw")]] <- metrics_raw
    summary_tbl[[paste0(model_name, "_downscaled")]] <- metrics_down
    
    write_csv(metrics_raw,  file.path(model_out, "metrics_raw_domain.csv"))
    write_csv(metrics_down, file.path(model_out, "metrics_downscaled_domain.csv"))
    
    ts_wide <- data.frame(
      model = model_name,
      variable = var_name,
      date = dates_val,
      year = as.integer(format(dates_val, "%Y")),
      month = as.integer(format(dates_val, "%m")),
      season = get_season(as.integer(format(dates_val, "%m"))),
      observed = obs_val_mean,
      raw = raw_val_mean,
      downscaled = down_val_mean
    )
    
    write_csv(ts_wide, file.path(model_out, "hindcast_domain_mean_timeseries.csv"))
    
    all_ts_tbl[[model_name]] <- ts_wide
    
    clim_df <- ts_wide %>%
      group_by(month) %>%
      summarise(
        observed = mean(observed, na.rm = TRUE),
        raw = mean(raw, na.rm = TRUE),
        downscaled = mean(downscaled, na.rm = TRUE),
        .groups = "drop"
      )
    
    write_csv(clim_df, file.path(model_out, "hindcast_monthly_climatology.csv"))
    
    ts_long <- ts_wide %>%
      select(date, observed, raw, downscaled) %>%
      pivot_longer(-date, names_to = "series", values_to = "value")
    
    p_ts <- ggplot(ts_long, aes(date, value, color = series)) +
      geom_line(linewidth = 0.7) +
      labs(
        title = paste("Hindcast validation:", model_name, var_name),
        x = NULL,
        y = y_label,
        color = "Series"
      ) +
      theme_bw()
    
    ggsave(
      file.path(model_out, "hindcast_timeseries.png"),
      p_ts,
      width = 9,
      height = 4,
      dpi = 300
    )
    
    clim_long <- clim_df %>%
      pivot_longer(
        cols = c(observed, raw, downscaled),
        names_to = "series",
        values_to = "value"
      )
    
    p_clim <- ggplot(clim_long, aes(month, value, color = series, group = series)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2) +
      scale_x_continuous(breaks = 1:12) +
      labs(
        title = paste("Monthly climatology:", model_name, var_name),
        x = "Month",
        y = y_label,
        color = "Series"
      ) +
      theme_bw()
    
    ggsave(
      file.path(model_out, "hindcast_monthly_climatology.png"),
      p_clim,
      width = 7,
      height = 4,
      dpi = 300
    )
    
    scatter_long <- ts_wide %>%
      select(observed, raw, downscaled) %>%
      pivot_longer(
        cols = c(raw, downscaled),
        names_to = "dataset",
        values_to = "simulated"
      )
    
    p_scatter <- ggplot(scatter_long, aes(x = observed, y = simulated, color = dataset)) +
      geom_point(alpha = 0.65, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
      labs(
        title = paste("Observed vs simulated:", model_name, var_name),
        x = paste("Observed", y_label),
        y = paste("Simulated", y_label),
        color = "Dataset"
      ) +
      theme_bw()
    
    ggsave(
      file.path(model_out, "scatter_observed_vs_simulated.png"),
      p_scatter,
      width = 6,
      height = 5,
      dpi = 300
    )
    
    cat("Finished:", model_name, "\n")
    
    rm(
      obs_on_gcm, gcm_aligned,
      obs_cal, gcm_cal, obs_val, gcm_val,
      gcm_cal_mat, obs_cal_mat, gcm_val_mat,
      hindcast_mat, hindcast_downscaled
    )
    gc()
  }
  
  if (length(summary_tbl) == 0 || length(all_ts_tbl) == 0) {
    warning("No valid results generated for variable: ", var_name)
    return(NULL)
  }
  
  summary_df <- bind_rows(summary_tbl) %>%
    select(model, variable, dataset, everything())
  
  write_csv(
    summary_df,
    file.path(out_root, paste0("hindcast_summary_metrics_", var_name, ".csv"))
  )
  
  all_ts_df <- bind_rows(all_ts_tbl)
  
  write_csv(
    all_ts_df,
    file.path(out_root, paste0("hindcast_domain_mean_timeseries_all_models_", var_name, ".csv"))
  )
  
  make_summary_plots(summary_df, all_ts_df, var_name, out_root, y_label)
  
  cat("\nFinished variable:", var_name, "\n")
}

# ============================================================
# 7. RUN ALL VARIABLES
# ============================================================

for (v in vars_to_run) {
  run_one_variable(v)
}

cat("\nAll hindcast validation completed.\n")
