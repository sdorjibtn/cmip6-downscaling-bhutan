# ============================================================
# PDF-style historical validation plots
#
# This script compares observed, raw historical GCM, and
# downscaled historical climate data by:
#   1. Reading observed monthly climate data
#   2. Reading raw historical GCM monthly data
#   3. Reading downscaled historical yearly GeoTIFF stacks
#   4. Aligning observed, raw, and downscaled rasters
#   5. Extracting Bhutan-domain monthly mean values
#   6. Aggregating values to seasonal and annual scales
#   7. Plotting probability density functions with ensemble mean ± 1 SD
#
# Notes:
#   - This script is for PDF-style historical validation plots only.
#   - It complements the metrics, scatter-density, and Taylor diagram scripts.
#   - Precipitation is aggregated as seasonal/annual totals.
#   - Temperature is aggregated as seasonal/annual means.
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("pdf.R")
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================

root <- Sys.getenv("CLIM_ROOT")

if (root == "") {
  root <- "~/Climate_downscaling/Downscaling"
}

var <- "pr"   # options: "pr", "tasmax", "tasmin"

obs_root  <- file.path(root, "data/observation/nc_files")
raw_root  <- file.path(root, "data/gcm_monthly")
down_root <- file.path(root, "outputs/delta_downscaling_pr_tmax_tmin")

out_root   <- file.path(root, "analysis/historical/pdf", var)
plot_root  <- file.path(out_root, "plots")
value_root <- file.path(out_root, "values")

dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
dir.create(value_root, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- file.path(root, "tmp_terra_pdf_validation", var)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

terraOptions(
  progress = 1,
  memfrac = 0.6,
  tempdir = tmp_dir,
  threads = 2
)

message("Running variable: ", var)
message("Project root: ", root)
message("Temporary directory: ", tmp_dir)

val_start <- "1985-05"
val_end   <- "2014-12"

# Use NULL to run all available models
test_model <- NULL
# test_model <- c("ACCESS-CM2", "MIROC6", "BCC-CSM2-MR")

obs_file <- switch(
  var,
  pr     = file.path(obs_root, "precip_3ds_monthly_198505202008.nc"),
  tasmax = file.path(obs_root, "tmax_3ds_m_monthly_198505202008.nc"),
  tasmin = file.path(obs_root, "tmin_3ds_tm_monthly_198505202008.nc"),
  stop("Unknown variable: ", var)
)

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

get_time_safe <- function(r) {
  tt <- time(r)
  if (is.null(tt)) stop("Raster has no time metadata.")
  as.Date(tt)
}

subset_year_month_by_index <- function(r, start_ym, end_ym, data_start_ym) {
  
  start_date <- as.Date(paste0(start_ym, "-01"))
  end_date   <- as.Date(paste0(end_ym, "-01"))
  data_start <- as.Date(paste0(data_start_ym, "-01"))
  
  all_months <- seq(data_start, by = "month", length.out = terra::nlyr(r))
  
  idx <- which(all_months >= start_date & all_months <= end_date)
  
  if (length(idx) == 0) {
    stop("No layers in requested date range using index-based subset.")
  }
  
  r2 <- r[[idx]]
  terra::time(r2) <- all_months[idx]
  names(r2) <- format(all_months[idx], "%Y_%m")
  
  r2
}

get_models <- function() {
  p <- file.path(down_root, var, "historical")
  if (!dir.exists(p)) return(character(0))
  list.dirs(p, full.names = FALSE, recursive = FALSE)
}

load_obs <- function() {
  if (!file.exists(obs_file)) {
    stop("Observation file not found: ", obs_file)
  }
  
  r <- rast(obs_file)
  crs(r) <- "EPSG:5266"
  r
}

load_raw_list <- function(models) {
  
  raw_list <- list()
  
  ff_all <- list.files(
    raw_root,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  for (m in models) {
    
    ff <- ff_all[
      grepl(m, ff_all) &
        grepl("historical", ff_all, ignore.case = TRUE) &
        grepl(var, ff_all, ignore.case = TRUE)
    ]
    
    if (length(ff) == 0) {
      cat("Missing RAW:", m, "\n")
      next
    }
    
    raw_list[[m]] <- rast(ff[1])
    cat("Loaded RAW:", m, "\n")
  }
  
  raw_list
}

load_downscaled_list <- function(models) {
  
  down_list <- list()
  
  for (m in models) {
    
    folder <- file.path(down_root, var, "historical", m)
    
    if (!dir.exists(folder)) {
      cat("Missing DOWN folder:", folder, "\n")
      next
    }
    
    files <- list.files(
      folder,
      pattern = "\\.tif$",
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      cat("Missing DOWN tif files:", folder, "\n")
      next
    }
    
    yrs <- sub(".*_(\\d{4})\\.tif$", "\\1", basename(files))
    keep <- !is.na(suppressWarnings(as.integer(yrs)))
    
    files <- files[keep]
    yrs   <- yrs[keep]
    
    ord <- order(as.integer(yrs))
    files <- files[ord]
    yrs   <- yrs[ord]
    
    r <- rast(files)
    
    dates <- seq(
      as.Date(paste0(min(as.integer(yrs)), "-01-01")),
      as.Date(paste0(max(as.integer(yrs)), "-12-01")),
      by = "month"
    )
    
    if (length(dates) != nlyr(r)) {
      stop("Date count does not match layer count for ", m)
    }
    
    time(r) <- dates
    names(r) <- format(dates, "%Y_%m")
    
    down_list[[m]] <- r
    cat("Loaded DOWN:", m, " -> ", length(files), " yearly tif files\n")
  }
  
  down_list
}

match_by_yearmonth <- function(obs_full, raw_full, down_full) {
  
  t_obs  <- get_time_safe(obs_full)
  t_raw  <- get_time_safe(raw_full)
  t_down <- get_time_safe(down_full)
  
  ym_obs  <- format(t_obs, "%Y-%m")
  ym_raw  <- format(t_raw, "%Y-%m")
  ym_down <- format(t_down, "%Y-%m")
  
  common_ym <- Reduce(intersect, list(ym_obs, ym_raw, ym_down))
  
  if (length(common_ym) == 0) {
    stop("No common year-months across observation, raw and downscaled.")
  }
  
  list(
    obs  = obs_full[[match(common_ym, ym_obs)]],
    raw  = raw_full[[match(common_ym, ym_raw)]],
    down = down_full[[match(common_ym, ym_down)]]
  )
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

get_season <- function(dates) {
  
  m <- format(dates, "%m")
  
  out <- ifelse(
    m %in% c("12", "01", "02"), "DJF",
    ifelse(
      m %in% c("03", "04", "05"), "MAM",
      ifelse(m %in% c("06", "07", "08"), "JJA", "SON")
    )
  )
  
  factor(out, levels = c("DJF", "MAM", "JJA", "SON"))
}

get_season_year <- function(dates) {
  
  y <- as.integer(format(dates, "%Y"))
  m <- format(dates, "%m")
  
  ifelse(m == "12", y + 1, y)
}

extract_monthly_bhutan_mean <- function(r) {
  
  vals <- global(r, mean, na.rm = TRUE)[, 1]
  dates <- get_time_safe(r)
  
  data.frame(
    date = dates,
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    season = get_season(dates),
    season_year = get_season_year(dates),
    value = vals
  )
}

monthly_to_seasonal <- function(df) {
  
  if (var == "pr") {
    df |>
      group_by(season, season_year) |>
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
  } else {
    df |>
      group_by(season, season_year) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
}

monthly_to_annual <- function(df) {
  
  if (var == "pr") {
    df |>
      group_by(year) |>
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
  } else {
    df |>
      group_by(year) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
}

make_value_tables <- function(obs_r, raw_list, down_list) {
  
  obs_monthly <- extract_monthly_bhutan_mean(obs_r)
  
  obs_seasonal <- monthly_to_seasonal(obs_monthly) |>
    mutate(
      dataset = "Observed",
      model = NA_character_,
      variable = var
    ) |>
    rename(period_year = season_year, group = season)
  
  obs_annual <- monthly_to_annual(obs_monthly) |>
    mutate(
      dataset = "Observed",
      model = NA_character_,
      variable = var,
      group = "Annual",
      period_year = year
    ) |>
    select(group, period_year, value, dataset, model, variable)
  
  raw_seasonal_all <- lapply(names(raw_list), function(m) {
    extract_monthly_bhutan_mean(raw_list[[m]]) |>
      monthly_to_seasonal() |>
      mutate(
        dataset = "Raw historical",
        model = m,
        variable = var
      ) |>
      rename(period_year = season_year, group = season)
  }) |> bind_rows()
  
  raw_annual_all <- lapply(names(raw_list), function(m) {
    extract_monthly_bhutan_mean(raw_list[[m]]) |>
      monthly_to_annual() |>
      mutate(
        dataset = "Raw historical",
        model = m,
        variable = var,
        group = "Annual",
        period_year = year
      ) |>
      select(group, period_year, value, dataset, model, variable)
  }) |> bind_rows()
  
  down_seasonal_all <- lapply(names(down_list), function(m) {
    extract_monthly_bhutan_mean(down_list[[m]]) |>
      monthly_to_seasonal() |>
      mutate(
        dataset = "Downscaled historical",
        model = m,
        variable = var
      ) |>
      rename(period_year = season_year, group = season)
  }) |> bind_rows()
  
  down_annual_all <- lapply(names(down_list), function(m) {
    extract_monthly_bhutan_mean(down_list[[m]]) |>
      monthly_to_annual() |>
      mutate(
        dataset = "Downscaled historical",
        model = m,
        variable = var,
        group = "Annual",
        period_year = year
      ) |>
      select(group, period_year, value, dataset, model, variable)
  }) |> bind_rows()
  
  seasonal_all <- bind_rows(
    obs_seasonal,
    raw_seasonal_all,
    down_seasonal_all
  ) |>
    mutate(group = factor(group, levels = c("DJF", "MAM", "JJA", "SON")))
  
  annual_all <- bind_rows(
    obs_annual,
    raw_annual_all,
    down_annual_all
  )
  
  list(seasonal = seasonal_all, annual = annual_all)
}

make_padded_grid <- function(x_all, n = 512) {
  
  x_all <- x_all[is.finite(x_all)]
  
  if (length(x_all) < 2) {
    return(seq(0, 1, length.out = n))
  }
  
  xr <- range(x_all, na.rm = TRUE)
  dx <- diff(xr)
  
  if (!is.finite(dx) || dx == 0) dx <- 1
  
  left_pad  <- 0.25 * dx
  right_pad <- 0.30 * dx
  
  seq(
    xr[1] - left_pad,
    xr[2] + right_pad,
    length.out = n
  )
}

density_on_grid <- function(x, grid, adjust = 1.2) {
  
  x <- x[is.finite(x)]
  
  if (length(unique(x)) < 2) {
    return(rep(NA_real_, length(grid)))
  }
  
  d <- density(
    x,
    from = min(grid),
    to = max(grid),
    n = length(grid),
    adjust = adjust,
    na.rm = TRUE
  )
  
  d$y
}

build_density_band <- function(df_vals, obs_vals, other_vals, groups) {
  
  out <- list()
  
  for (g in groups) {
    
    df_g    <- df_vals[df_vals$group == g, ]
    obs_g   <- obs_vals[obs_vals$group == g, ]
    other_g <- other_vals[other_vals$group == g, ]
    
    if (nrow(df_g) == 0) next
    
    x_all <- c(
      df_g$value[is.finite(df_g$value)],
      obs_g$value[is.finite(obs_g$value)],
      other_g$value[is.finite(other_g$value)]
    )
    
    if (length(unique(x_all)) < 2) next
    
    grid <- make_padded_grid(x_all, n = 512)
    
    model_ids <- unique(df_g$model)
    
    dens_mat <- sapply(model_ids, function(mm) {
      density_on_grid(
        df_g$value[df_g$model == mm],
        grid,
        adjust = 1.2
      )
    })
    
    if (is.vector(dens_mat)) {
      dens_mat <- matrix(dens_mat, ncol = 1)
    }
    
    out[[as.character(g)]] <- data.frame(
      group = g,
      x = grid,
      y_mean = rowMeans(dens_mat, na.rm = TRUE),
      y_sd = apply(dens_mat, 1, sd, na.rm = TRUE)
    )
  }
  
  bind_rows(out)
}

build_obs_density <- function(obs_vals, groups) {
  
  out <- list()
  
  for (g in groups) {
    
    x <- obs_vals$value[obs_vals$group == g]
    x <- x[is.finite(x)]
    
    if (length(unique(x)) < 2) next
    
    xr <- range(x, na.rm = TRUE)
    dx <- diff(xr)
    
    if (!is.finite(dx) || dx == 0) dx <- 1
    
    left_pad  <- 0.25 * dx
    right_pad <- 0.30 * dx
    
    d <- density(
      x,
      from = xr[1] - left_pad,
      to = xr[2] + right_pad,
      n = 512,
      adjust = 1.3,
      na.rm = TRUE
    )
    
    out[[as.character(g)]] <- data.frame(
      group = g,
      x = d$x,
      y = d$y
    )
  }
  
  bind_rows(out)
}

make_pdf <- function(
    obs_r,
    raw_list,
    down_list,
    out_file,
    time_scale = c("seasonal", "annual")
) {
  
  time_scale <- match.arg(time_scale)
  
  obs_monthly <- extract_monthly_bhutan_mean(obs_r)
  
  if (time_scale == "seasonal") {
    
    obs_vals <- monthly_to_seasonal(obs_monthly) |>
      mutate(group = factor(season, levels = c("DJF", "MAM", "JJA", "SON")))
    
    group_levels <- c("DJF", "MAM", "JJA", "SON")
    
    x_lab <- if (var == "pr") {
      "Seasonal precipitation (mm)"
    } else {
      "Seasonal temperature (°C)"
    }
    
  } else {
    
    obs_vals <- monthly_to_annual(obs_monthly) |>
      mutate(group = factor("Annual", levels = "Annual"))
    
    group_levels <- "Annual"
    
    x_lab <- if (var == "pr") {
      "Annual precipitation (mm)"
    } else {
      "Annual temperature (°C)"
    }
  }
  
  raw_vals_all <- lapply(names(raw_list), function(m) {
    
    x <- extract_monthly_bhutan_mean(raw_list[[m]])
    
    if (time_scale == "seasonal") {
      monthly_to_seasonal(x) |>
        mutate(group = factor(season, levels = group_levels), model = m)
    } else {
      monthly_to_annual(x) |>
        mutate(group = factor("Annual", levels = group_levels), model = m)
    }
    
  }) |> bind_rows()
  
  down_vals_all <- lapply(names(down_list), function(m) {
    
    x <- extract_monthly_bhutan_mean(down_list[[m]])
    
    if (time_scale == "seasonal") {
      monthly_to_seasonal(x) |>
        mutate(group = factor(season, levels = group_levels), model = m)
    } else {
      monthly_to_annual(x) |>
        mutate(group = factor("Annual", levels = group_levels), model = m)
    }
    
  }) |> bind_rows()
  
  raw_band    <- build_density_band(raw_vals_all, obs_vals, down_vals_all, group_levels)
  down_band   <- build_density_band(down_vals_all, obs_vals, raw_vals_all, group_levels)
  obs_density <- build_obs_density(obs_vals, group_levels)
  
  p <- ggplot() +
    
    geom_ribbon(
      data = raw_band,
      aes(
        x = x,
        ymin = pmax(0, y_mean - y_sd),
        ymax = y_mean + y_sd,
        fill = "Raw GCM"
      ),
      alpha = 0.15
    ) +
    
    geom_line(
      data = raw_band,
      aes(x = x, y = y_mean, colour = "Raw GCM"),
      linewidth = 1
    ) +
    
    geom_ribbon(
      data = down_band,
      aes(
        x = x,
        ymin = pmax(0, y_mean - y_sd),
        ymax = y_mean + y_sd,
        fill = "Downscaled GCM"
      ),
      alpha = 0.15
    ) +
    
    geom_line(
      data = down_band,
      aes(x = x, y = y_mean, colour = "Downscaled GCM"),
      linewidth = 1
    ) +
    
    geom_line(
      data = obs_density,
      aes(x = x, y = y, colour = "Observed"),
      linewidth = 1,
      linetype = "dashed"
    ) +
    
    facet_wrap(
      ~group,
      scales = "free",
      labeller = labeller(
        group = c(
          DJF = "DJF (Winter)",
          MAM = "MAM (Spring)",
          JJA = "JJA (Summer)",
          SON = "SON (Autumn)",
          Annual = "Annual"
        )
      )
    ) +
    
    # scale_colour_manual(
    #   name = NULL,
    #   values = c(
    #     "Raw GCM" = "red",
    #     "Downscaled GCM" = "blue",
    #     "Observed" = "black"
    #   )
    # ) +
    # 
    # scale_fill_manual(
    #   name = NULL,
    #   values = c(
    #     "Raw GCM" = "red",
    #     "Downscaled GCM" = "blue"
    #   )
    # ) +
    # Using colour blind friendly
    scale_colour_manual(
      name = NULL,
      values = c(
        "Raw GCM" = "#D55E00",          # vermillion
        "Downscaled GCM" = "#0072B2",   # blue
        "Observed" = "black"
      )
    ) +
    
    scale_fill_manual(
      name = NULL,
      values = c(
        "Raw GCM" = "#D55E00",
        "Downscaled GCM" = "#0072B2"
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::extended_breaks(n = 8),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(family = "serif", size = 20),
      legend.key.size = unit(1.5, "lines"),
      legend.key.width = unit(2, "lines"),
      plot.title = element_blank(),   # ensures no title even if added later
      axis.text = element_text(family = "serif", size = 20),
      axis.title = element_text(family = "serif", size = 20),
      strip.background = element_blank(),
      strip.text = element_text(family = "serif", size = 20),
      panel.grid.minor = element_blank()
    ) +
    
    labs(
      x = x_lab,
      y = "PDF"
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 12,
    height = 9,
    units = "in",
    dpi = 600
  )
}

# ============================================================
# 3. RUN ANALYSIS
# ============================================================

obs_full <- load_obs()

obs_full <- subset_year_month_by_index(
  r = obs_full,
  start_ym = val_start,
  end_ym = val_end,
  data_start_ym = "1985-05"
)

models <- get_models()

if (length(models) == 0) {
  stop(paste("No model folders found for", var))
}

if (!is.null(test_model)) {
  models <- intersect(models, test_model)
  
  if (length(models) == 0) {
    stop(paste("Requested test_model not found ->", paste(test_model, collapse = ", ")))
  }
}

message("Models to process:")
print(models)

raw_loaded  <- load_raw_list(models)
down_loaded <- load_downscaled_list(models)

common_models <- intersect(names(raw_loaded), names(down_loaded))

if (length(common_models) == 0) {
  stop(paste("No common models loaded for", var))
}

message("Common models loaded:")
print(common_models)

raw_list  <- list()
down_list <- list()
obs_aligned <- NULL

for (m in common_models) {
  
  cat("Matching and aligning model:", m, "\n")
  
  raw_full <- subset_year_month_by_index(
    r = raw_loaded[[m]],
    start_ym = val_start,
    end_ym = val_end,
    data_start_ym = "1980-01"
  )
  
  down_full <- subset_year_month_by_index(
    r = down_loaded[[m]],
    start_ym = val_start,
    end_ym = val_end,
    data_start_ym = "1980-01"
  )
  
  matched <- match_by_yearmonth(
    obs_full = obs_full,
    raw_full = raw_full,
    down_full = down_full
  )
  
  aligned <- align_three(
    obs = matched$obs,
    raw = matched$raw,
    down = matched$down
  )
  
  raw_list[[m]]  <- aligned$raw
  down_list[[m]] <- aligned$down
  obs_aligned <- aligned$obs
}

# ============================================================
# 4. SAVE VALUES
# ============================================================

vals <- make_value_tables(
  obs_r = obs_aligned,
  raw_list = raw_list,
  down_list = down_list
)

write.csv(
  vals$seasonal,
  file.path(value_root, paste0(var, "_seasonal_values.csv")),
  row.names = FALSE
)

write.csv(
  vals$annual,
  file.path(value_root, paste0(var, "_annual_values.csv")),
  row.names = FALSE
)

# ============================================================
# 5. MAKE PDF-STYLE PLOTS
# ============================================================

make_pdf(
  obs_r = obs_aligned,
  raw_list = raw_list,
  down_list = down_list,
  out_file = file.path(plot_root, paste0(var, "_seasonal_pdf_bhutan.png")),
  time_scale = "seasonal"
)

make_pdf(
  obs_r = obs_aligned,
  raw_list = raw_list,
  down_list = down_list,
  out_file = file.path(plot_root, paste0(var, "_annual_pdf_bhutan.png")),
  time_scale = "annual"
)

tmpFiles(remove = TRUE)

message("Done for variable: ", var)
message("Plots saved to: ", plot_root)
message("CSV values saved to: ", value_root)

