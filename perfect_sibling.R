# ============================================================
# Perfect sibling validation for modified delta downscaling
#
# This script evaluates the robustness of the delta-change method by:
#   1. Treating one GCM as the reference "sibling"
#   2. Treating another GCM as the donor model
#   3. Computing donor future/historical monthly anomalies
#   4. Applying those anomalies to the sibling historical climatology
#   5. Comparing raw donor future and delta-corrected future against
#      the sibling future climatology
#   6. Saving seasonal RMSE/correlation results and diagnostic plots
#
# Notes:
#   - This is a perfect sibling validation script.
#   - It uses GCM-to-GCM comparison, not observations.
#   - Precipitation uses multiplicative delta correction.
#   - Temperature uses additive delta correction.
#
# To run:
#   Sys.setenv(CLIM_ROOT = "path/to/project")
#   source("perfect_sibling.R")
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

terraOptions(progress = 1, memfrac = 0.6)

# ============================================================
# 1. SETTINGS
# ============================================================

root <- Sys.getenv("CLIM_ROOT")

if (root == "") {
  root <- "~/Climate_downscaling/Downscaling"
}

var <- "pr"   # options: "pr", "tasmax", "tasmin"

ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")
# ssps <- c("ssp126")   # use this line for testing only

gcm_root <- file.path(root, "data/gcm_monthly")

# Historical and future periods
hist_start <- as.Date("1986-01-01")
hist_end   <- as.Date("2014-12-01")

fut_start  <- as.Date("2041-01-01")
fut_end    <- as.Date("2060-12-01")

period_label <- paste0(
  format(hist_start, "%Y"), "-", format(hist_end, "%Y"),
  "_to_",
  format(fut_start, "%Y"), "-", format(fut_end, "%Y")
)

out_root <- file.path(
  root,
  "analysis/historical/perfect_sibling_check",
  period_label
)

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. FUNCTIONS
# ============================================================

assign_time <- function(r, start_date) {
  time(r) <- seq(as.Date(start_date), by = "month", length.out = nlyr(r))
  r
}

subset_dates <- function(r, start_date, end_date) {
  tt <- as.Date(time(r))
  
  if (is.null(tt)) {
    stop("Raster has no time metadata.")
  }
  
  idx <- which(tt >= start_date & tt <= end_date)
  
  if (length(idx) == 0) {
    stop("No layers found for requested date range.")
  }
  
  r[[idx]]
}

find_files <- function(model, type, ssp = NULL) {
  
  ff <- list.files(
    gcm_root,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (type == "historical") {
    
    ff <- ff[
      grepl(model, ff, fixed = TRUE) &
        grepl("historical", ff, ignore.case = TRUE) &
        grepl(var, ff, ignore.case = TRUE)
    ]
    
  } else {
    
    ff <- ff[
      grepl(model, ff, fixed = TRUE) &
        grepl(ssp, ff, ignore.case = TRUE) &
        grepl(var, ff, ignore.case = TRUE)
    ]
  }
  
  if (length(ff) == 0) return(NA_character_)
  ff[1]
}

monthly_clim <- function(r) {
  
  tt <- as.Date(time(r))
  months <- as.integer(format(tt, "%m"))
  
  out <- lapply(1:12, function(m) {
    mean(r[[which(months == m)]], na.rm = TRUE)
  })
  
  out <- rast(out)
  names(out) <- month.abb
  
  out
}

seasonal <- function(r, var_name) {
  
  if (nlyr(r) != 12) {
    stop("Seasonal function expects a 12-layer monthly climatology.")
  }
  
  if (var_name == "pr") {
    djf <- r[[12]] + r[[1]] + r[[2]]
    mam <- r[[3]]  + r[[4]] + r[[5]]
    jja <- r[[6]]  + r[[7]] + r[[8]]
    son <- r[[9]]  + r[[10]] + r[[11]]
  } else {
    djf <- (r[[12]] + r[[1]] + r[[2]]) / 3
    mam <- (r[[3]]  + r[[4]] + r[[5]]) / 3
    jja <- (r[[6]]  + r[[7]] + r[[8]]) / 3
    son <- (r[[9]]  + r[[10]] + r[[11]]) / 3
  }
  
  s <- c(djf, mam, jja, son)
  names(s) <- c("DJF", "MAM", "JJA", "SON")
  s
}

metrics <- function(obs, sim) {
  
  o <- values(obs)
  s <- values(sim)
  
  ok <- is.finite(o) & is.finite(s)
  
  o <- o[ok]
  s <- s[ok]
  
  if (length(o) == 0) {
    return(data.frame(
      rmse = NA_real_,
      mae  = NA_real_,
      bias = NA_real_,
      cor  = NA_real_
    ))
  }
  
  data.frame(
    rmse = sqrt(mean((s - o)^2)),
    mae  = mean(abs(s - o)),
    bias = mean(s - o),
    cor  = suppressWarnings(cor(o, s))
  )
}

get_models <- function() {
  
  files <- list.files(
    gcm_root,
    pattern = "\\.nc$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  files <- files[
    grepl("historical", files, ignore.case = TRUE) &
      grepl(var, files, ignore.case = TRUE)
  ]
  
  if (length(files) == 0) {
    stop("No historical files found for variable: ", var)
  }
  
  models <- unique(
    sub(".*_([A-Za-z0-9.-]+)_historical.*", "\\1", basename(files))
  )
  
  sort(models)
}

# ============================================================
# 3. GET MODELS
# ============================================================

models <- get_models()

# For testing only:
# models <- c("ACCESS-CM2", "BCC-CSM2-MR")

cat("\nModels found:\n")
print(models)

# ============================================================
# 4. PERFECT SIBLING LOOP
# ============================================================

all_results <- list()
fail_log <- list()

for (ssp in ssps) {
  
  cat("\n==============================\n")
  cat("SCENARIO:", ssp, "\n")
  cat("==============================\n")
  
  for (sib in models) {
    
    cat("\nSibling:", sib, "\n")
    
    sib_hist_file <- find_files(sib, "historical")
    sib_fut_file  <- find_files(sib, "future", ssp)
    
    if (is.na(sib_hist_file) || is.na(sib_fut_file)) {
      fail_log[[length(fail_log) + 1]] <- data.frame(
        ssp = ssp,
        sibling = sib,
        donor = NA_character_,
        reason = "Sibling historical or future file missing"
      )
      next
    }
    
    tryCatch({
      
      sib_hist <- rast(sib_hist_file)
      sib_fut  <- rast(sib_fut_file)
      
      sib_hist <- assign_time(sib_hist, "1980-01-01")
      sib_fut  <- assign_time(sib_fut,  "2015-01-01")
      
      sib_hist <- subset_dates(sib_hist, hist_start, hist_end)
      sib_fut  <- subset_dates(sib_fut,  fut_start, fut_end)
      
      sib_hist_clim <- monthly_clim(sib_hist)
      sib_fut_clim  <- monthly_clim(sib_fut)
      
      donors <- setdiff(models, sib)
      
      for (don in donors) {
        
        cat("  Donor:", don, "\n")
        
        d_hist_file <- find_files(don, "historical")
        d_fut_file  <- find_files(don, "future", ssp)
        
        if (is.na(d_hist_file) || is.na(d_fut_file)) {
          fail_log[[length(fail_log) + 1]] <- data.frame(
            ssp = ssp,
            sibling = sib,
            donor = don,
            reason = "Donor historical or future file missing"
          )
          next
        }
        
        tryCatch({
          
          d_hist <- rast(d_hist_file)
          d_fut  <- rast(d_fut_file)
          
          d_hist <- assign_time(d_hist, "1980-01-01")
          d_fut  <- assign_time(d_fut,  "2015-01-01")
          
          d_hist <- subset_dates(d_hist, hist_start, hist_end)
          d_fut  <- subset_dates(d_fut,  fut_start, fut_end)
          
          # Align donor grids to sibling grids
          d_hist <- resample(d_hist, sib_hist, method = "bilinear")
          d_fut  <- resample(d_fut,  sib_fut,  method = "bilinear")
          
          d_hist_clim <- monthly_clim(d_hist)
          d_fut_clim  <- monthly_clim(d_fut)
          
          # Delta correction
          if (var == "pr") {
            
            delta <- d_fut_clim / d_hist_clim
            delta[!is.finite(delta)] <- 1
            
            dc <- sib_hist_clim * delta
            
          } else {
            
            delta <- d_fut_clim - d_hist_clim
            dc <- sib_hist_clim + delta
          }
          
          ref_fut <- seasonal(sib_fut_clim, var)
          raw_fut <- seasonal(d_fut_clim, var)
          dc_fut  <- seasonal(dc, var)
          
          for (s in names(ref_fut)) {
            
            m_raw <- metrics(ref_fut[[s]], raw_fut[[s]])
            m_dc  <- metrics(ref_fut[[s]], dc_fut[[s]])
            
            all_results[[length(all_results) + 1]] <- data.frame(
              variable = var,
              ssp = ssp,
              sibling = sib,
              donor = don,
              season = s,
              raw_rmse = m_raw$rmse,
              dc_rmse  = m_dc$rmse,
              raw_mae  = m_raw$mae,
              dc_mae   = m_dc$mae,
              raw_bias = m_raw$bias,
              dc_bias  = m_dc$bias,
              raw_cor  = m_raw$cor,
              dc_cor   = m_dc$cor
            )
          }
          
        }, error = function(e) {
          
          cat("    ERROR donor:", e$message, "\n")
          
          fail_log[[length(fail_log) + 1]] <<- data.frame(
            ssp = ssp,
            sibling = sib,
            donor = don,
            reason = e$message
          )
        })
      }
      
    }, error = function(e) {
      
      cat("  ERROR sibling:", e$message, "\n")
      
      fail_log[[length(fail_log) + 1]] <<- data.frame(
        ssp = ssp,
        sibling = sib,
        donor = NA_character_,
        reason = e$message
      )
    })
  }
}

if (length(all_results) == 0) {
  stop("No perfect sibling results generated. Check file paths and model names.")
}

results <- bind_rows(all_results)

results_csv <- file.path(
  out_root,
  paste0(var, "_perfect_sibling_results_", period_label, ".csv")
)

write.csv(results, results_csv, row.names = FALSE)

# ============================================================
# 5. PLOT RMSE VS CORRELATION BY SSP
# ============================================================

plot_df <- results |>
  pivot_longer(
    cols = c(raw_rmse, dc_rmse, raw_cor, dc_cor),
    names_to = c("dataset", ".value"),
    names_pattern = "(raw|dc)_(rmse|cor)"
  ) |>
  mutate(
    dataset = recode(
      dataset,
      raw = "Raw future",
      dc  = "DC future"
    ),
    dataset = factor(dataset, levels = c("Raw future", "DC future")),
    ssp_label = recode(
      ssp,
      "ssp126" = "SSP1-2.6",
      "ssp245" = "SSP2-4.5",
      "ssp370" = "SSP3-7.0",
      "ssp585" = "SSP5-8.5"
    ),
    ssp_label = factor(
      ssp_label,
      levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    )
  )

p <- ggplot(plot_df, aes(x = rmse, y = cor, colour = dataset)) +
  geom_point(alpha = 0.35, size = 1.4) +
  facet_wrap(~ssp_label, ncol = 2) +
  scale_colour_manual(
    values = c(
      "Raw future" = "#D55E00",
      "DC future"  = "#0072B2"
    )
  ) +
  theme_bw(base_size = 14, base_family = "serif") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, family = "serif"),
    axis.text = element_text(size = 15, family = "serif"),
    axis.title = element_text(size = 16, family = "serif")
  ) +
  labs(
    x = if (var == "pr") {
      expression("RMSE (mm season"^{-1} * ")")
    } else {
      expression("RMSE (" * degree*C * ")")
    },
    y = "Pearson correlation"
  )

plot_png <- file.path(
  out_root,
  paste0(var, "_perfect_sibling_rmse_correlation_by_ssp_", period_label, ".png")
)

ggsave(
  plot_png,
  p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600
)

# ============================================================
# 6. SAVE FAIL LOG
# ============================================================

if (length(fail_log) > 0) {
  
  fail_tbl <- bind_rows(fail_log)
  
  write.csv(
    fail_tbl,
    file.path(out_root, paste0(var, "_perfect_sibling_fail_log_", period_label, ".csv")),
    row.names = FALSE
  )
}

cat("\n====================================================\n")
cat("DONE. Perfect sibling validation completed.\n")
cat("Results saved to:", results_csv, "\n")
cat("Plot saved to:", plot_png, "\n")
cat("Output folder:", out_root, "\n")
cat("====================================================\n")
