generate_weather <- function(precip_xts, temp_xts, nyears = 100,
                             prcp_range = c(0, Inf),
                             temp_range = c(-60, 60),
                             use_bootstrap = TRUE, # bootstrap vs param
                             roll_window = 31, # for param smoothing
                             prior_counts = 1, # pseudo-count for P01/P11
                             enforce_truncate = TRUE, # cap param draws at 99th pct
                             ar_phi = NULL, # if not NULL, apply AR1 on log-int
                             bootstrap_params = NULL)
{
  # Load required packages
  load_or_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' not found. Installing...", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Required packages
  pkgs <- c("xts", "dplyr", "lubridate")
  invisible(lapply(pkgs, load_or_install))
  
  # Optional package
  if (!requireNamespace("moments", quietly = TRUE)) {
    warning("Package 'moments' not available, AR(1) estimation may be limited")
  }
  
  # Function to estimate AR(1) parameter for log-intensities
  estimate_ar_phi <- function(precip_data) {
    wet_intensities <- precip_data[precip_data > 0]
    if (length(wet_intensities) < 10) return(0)
    
    log_intensities <- log(wet_intensities + 1e-5)
    
    if (length(log_intensities) > 2) {
      x1 <- log_intensities[-length(log_intensities)]
      x2 <- log_intensities[-1]
      phi <- cor(x1, x2, use = "complete.obs")
      phi <- max(-0.9, min(0.9, phi))
      return(phi)
    }
    return(0)
  }
  
  # Simplified bootstrap sampling - uniform method only
  sample_bootstrap <- function(pool, seasonal_factor = 1.0) {
    if (length(pool) == 0) return(0)
    if (length(pool) == 1) return(pool[1] * seasonal_factor)
    
    return(sample(pool, 1) * seasonal_factor)
  }
  
  stopifnot(inherits(index(precip_xts), "Date"),
            inherits(index(temp_xts), "Date"))
  
  # 1. Build historical df with day-of-year
  df <- merge(precip = precip_xts, temp = temp_xts)
  colnames(df) <- c("PRCP", "TEMP")
  df$month <- month(index(df))
  df$doy <- yday(index(df))
  
  # 2. Estimate AR(1) parameter if not provided
  if (is.null(ar_phi)) {
    ar_phi <- estimate_ar_phi(as.numeric(df$PRCP))
    message(sprintf("Estimated AR(1) parameter: %.3f", ar_phi))
  }
  
  # 3. Create enhanced seasonal intensity scaling
  # Convert xts to data frame for dplyr operations
  df_regular <- data.frame(
    date = index(df),
    PRCP = as.numeric(df$PRCP),
    TEMP = as.numeric(df$TEMP),
    month = as.numeric(df$month),
    doy = as.numeric(df$doy)
  )
  
  # Calculate daily mean intensities from observed data
  daily_intensity_stats <- df_regular %>%
    dplyr::filter(PRCP > 0) %>%
    dplyr::group_by(doy) %>%
    dplyr::summarise(
      daily_mean_intensity = mean(PRCP, na.rm = TRUE),
      daily_median_intensity = median(PRCP, na.rm = TRUE),
      daily_90p_intensity = quantile(PRCP, 0.9, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Smooth the daily statistics to avoid noise - use wider span for more stability
  smooth_daily_stats <- daily_intensity_stats %>%
    mutate(
      smooth_mean = fitted(loess(daily_mean_intensity ~ doy, span = 0.15)),
      smooth_median = fitted(loess(daily_median_intensity ~ doy, span = 0.15)),
      smooth_90p = fitted(loess(daily_90p_intensity ~ doy, span = 0.15))
    )
  
  # Handle missing DOYs (leap year days) by interpolation
  all_doys <- data.frame(doy = 1:366)
  smooth_daily_stats <- all_doys %>%
    left_join(smooth_daily_stats, by = "doy") %>%
    mutate(
      smooth_mean = ifelse(is.na(smooth_mean), 
                           approx(doy[!is.na(smooth_mean)], smooth_mean[!is.na(smooth_mean)], doy)$y,
                           smooth_mean),
      smooth_median = ifelse(is.na(smooth_median), 
                             approx(doy[!is.na(smooth_median)], smooth_median[!is.na(smooth_median)], doy)$y,
                             smooth_median),
      smooth_90p = ifelse(is.na(smooth_90p), 
                          approx(doy[!is.na(smooth_90p)], smooth_90p[!is.na(smooth_90p)], doy)$y,
                          smooth_90p)
    )
  
  # Calculate overall means for scaling
  overall_mean <- mean(df_regular$PRCP[df_regular$PRCP > 0], na.rm = TRUE)
  
  # 4. Precompute empirical wet-day intensities per month
  empirical_wet <- lapply(1:12, function(m) {
    month_subset <- df_regular[df_regular$month == m & df_regular$PRCP > 0, ]
    month_subset$PRCP
  })
  names(empirical_wet) <- month.abb
  
  # 5. Compute Markov transition probs
  prcp_params <- lapply(1:12, function(m) {
    sub <- df_regular[df_regular$month == m, ]
    wet <- as.integer(sub$PRCP > 0)
    if (length(wet) < 2) {
      return(list(P01 = 0.2, P11 = 0.8))
    }
    w0 <- wet[-length(wet)]; w1 <- wet[-1]
    n01 <- sum(w0==0 & w1==1); n0 <- sum(w0==0)
    n11 <- sum(w0==1 & w1==1); n1 <- sum(w0==1)
    P01 <- (n01 + prior_counts) / (n0 + 2*prior_counts)
    P11 <- (n11 + prior_counts) / (n1 + 2*prior_counts)
    list(P01 = P01, P11 = P11)
  })
  names(prcp_params) <- month.abb
  
  # 6. Calculate monthly log-intensity parameters
  monthly_log_params <- lapply(1:12, function(m) {
    month_data <- df_regular[df_regular$month == m, ]
    wet_log_intens <- log(month_data$PRCP[month_data$PRCP > 0] + 1e-5)
    if (length(wet_log_intens) > 1) {
      list(mu = mean(wet_log_intens, na.rm = TRUE),
           sd = sd(wet_log_intens, na.rm = TRUE))
    } else {
      list(mu = 0, sd = 1)
    }
  })
  names(monthly_log_params) <- month.abb
  
  # 7. Create synthetic date index
  start_date <- as.Date("2000-01-01")
  dates <- seq(start_date, by = "day", length.out = nyears*365 + 
                 sum(leap_year(year(start_date) + 0:(nyears-1))))
  n <- length(dates)
  mon <- month(dates)
  doy <- yday(dates)
  
  # 8. Simulate wet/dry sequence
  prcp <- numeric(n); wet <- integer(n)
  p_obs <- mean(df_regular$PRCP[df_regular$month == mon[1]] > 0, na.rm=TRUE)
  wet[1] <- rbinom(1, 1, p_obs)
  
  for (i in 2:n) {
    p_info <- prcp_params[[month.abb[mon[i]]]]
    p <- if (wet[i-1]==1) p_info$P11 else p_info$P01
    wet[i] <- rbinom(1, 1, p)
  }
  
  # 9. Simulate intensities with seasonal scaling
  for (i in which(wet==1)) {
    m <- month.abb[mon[i]]
    current_doy <- doy[i]
    
    if (use_bootstrap) {
      pool <- empirical_wet[[m]]
      
      # Get seasonal scaling factor with better handling
      current_doy_adj <- min(current_doy, 366)  # Handle leap years
      seasonal_stats <- smooth_daily_stats[smooth_daily_stats$doy == current_doy_adj, ]
      
      if (nrow(seasonal_stats) == 0 || is.na(seasonal_stats$smooth_mean)) {
        seasonal_factor <- 1.0  # Fallback
      } else {
        seasonal_factor <- seasonal_stats$smooth_mean / overall_mean
        # More aggressive scaling but with safety bounds
        seasonal_factor <- max(0.3, min(3.0, seasonal_factor))
      }
      
      # Uniform bootstrap sampling only
      prcp[i] <- sample_bootstrap(pool, seasonal_factor = seasonal_factor)
      
    } else {
      # Parametric approach with seasonal adjustment
      params <- monthly_log_params[[m]]
      
      # Apply seasonal scaling to the mean with better handling
      current_doy_adj <- min(current_doy, 366)
      seasonal_stats <- smooth_daily_stats[smooth_daily_stats$doy == current_doy_adj, ]
      
      if (nrow(seasonal_stats) == 0 || is.na(seasonal_stats$smooth_mean)) {
        seasonal_factor <- 1.0
      } else {
        seasonal_factor <- seasonal_stats$smooth_mean / overall_mean
        seasonal_factor <- max(0.3, min(3.0, seasonal_factor))
      }
      
      adjusted_mu <- params$mu + log(seasonal_factor)
      x <- rnorm(1, mean = adjusted_mu, sd = params$sd)
      
      if (enforce_truncate) {
        cap <- quantile(empirical_wet[[m]], 0.99, na.rm=TRUE)
        x <- min(x, log(cap + 1e-5))
      }
      prcp[i] <- exp(x)
    }
  }
  
  # 10. Apply AR(1) on log-intensity if ar_phi > 0
  if (!is.null(ar_phi) && abs(ar_phi) > 0.01) {
    logp <- log(prcp + 1e-5)
    first_wet <- which(wet==1)[1]
    
    if (!is.na(first_wet)) {
      logp_mean <- mean(logp[wet==1], na.rm = TRUE)
      logp_sd <- sd(logp[wet==1], na.rm = TRUE)
      
      for (i in 2:n) {
        if (wet[i]==1) {
          if (wet[i-1]==1) {
            innovation_sd <- sqrt(1 - ar_phi^2) * logp_sd
            logp[i] <- ar_phi * logp[i-1] + (1 - ar_phi) * logp_mean + 
              rnorm(1, 0, innovation_sd)
          }
          prcp[i] <- exp(logp[i])
        }
      }
    }
  }
  
  prcp <- pmax(pmin(prcp, prcp_range[2]), prcp_range[1])
  
  # 11. Simulate temperature (unchanged)
  temp <- numeric(n)
  temp_params <- t(sapply(1:12, function(m) {
    temps <- df_regular$TEMP[df_regular$month==m]
    c(mu = mean(temps, na.rm=TRUE), sd = sd(temps, na.rm=TRUE))
  }))
  rownames(temp_params) <- 1:12
  
  # Simple monthly temperature simulation (can be enhanced similarly)
  for (i in seq_len(n)) {
    m <- mon[i]
    temp[i] <- rnorm(1, mean = temp_params[m, "mu"], sd = temp_params[m, "sd"])
  }
  temp <- pmax(pmin(temp, temp_range[2]), temp_range[1])
  
  # 12. Return xts
  out <- xts::xts(cbind(PRECIP = prcp, TEMP = temp), order.by = dates)
  return(out)
}
