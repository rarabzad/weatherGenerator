generate_weather <- function(precip_xts, temp_xts, nyears = 100,
                                 prcp_range = c(0, Inf),
                                 temp_range = c(-60, 60),
                                 use_bootstrap = TRUE,      # bootstrap vs param
                                 roll_window = 31,          # for param smoothing
                                 prior_counts = 1,          # pseudo-count for P01/P11
                                 enforce_truncate = TRUE,   # cap param draws at 99th pct
                                 ar_phi = NULL              # if not NULL, apply AR1 on log-int
)
{
  stopifnot(inherits(index(precip_xts), "Date"),
            inherits(index(temp_xts),   "Date"))
  
  # 1. Build historical df
  df <- merge(precip = precip_xts, temp = temp_xts)
  colnames(df) <- c("PRCP", "TEMP")
  df$month <- month(index(df))
  
  # 2. Precompute empirical wet-day intensities per month
  empirical_wet <- lapply(1:12, function(m) {
    as.numeric(df$PRCP[df$month == m & df$PRCP > 0])
  })
  names(empirical_wet) <- month.abb
  
  # 3. Compute Markov transition probs with add-1 smoothing
  prcp_params <- lapply(1:12, function(m) {
    sub <- df[df$month == m, ]
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
  
  # 4. (Optional) Parametric log-mu/log-sd with rolling window
  log_intens <- log(df$PRCP + 1e-5)
  mu_roll <- rollapply(log_intens, width = roll_window,
                       FUN = function(x) mean(x[x>0], na.rm=TRUE),
                       fill = NA, align = "center")
  sd_roll <- rollapply(log_intens, width = roll_window,
                       FUN = function(x) sd(x[x>0], na.rm=TRUE),
                       fill = NA, align = "center")
  
  # 5. Create synthetic date index
  start_date <- as.Date("2000-01-01")
  dates <- seq(start_date, by = "day", length.out = nyears*365 + 
                 sum(leap_year(start_date + 0:(nyears*365-1))))
  n <- length(dates)
  mon <- month(dates)
  
  # 6. Simulate wet/dry sequence
  prcp <- numeric(n); wet <- integer(n)
  first_month <- as.character(month.abb[mon[1]])
  p_obs <- mean(df$PRCP[df$month == mon[1]] > 0, na.rm=TRUE)
  wet[1] <- rbinom(1,1,p_obs)
  
  for (i in 2:n) {
    p_info <- prcp_params[[ month.abb[mon[i]] ]]
    p <- if (wet[i-1]==1) p_info$P11 else p_info$P01
    wet[i] <- rbinom(1,1,p)
  }
  
  # 7. Simulate intensities
  for (i in which(wet==1)) {
    m <- month.abb[mon[i]]
    
    if (use_bootstrap) {
      # Empirical bootstrap
      pool <- empirical_wet[[m]]
      prcp[i] <- if (length(pool)>0) sample(pool,1) else 0
      
    } else {
      # Parametric log-normal with rolling mu/sd
      mu <- mu_roll[i]; sdv <- sd_roll[i]
      x <- rnorm(1, mean = mu, sd = sdv)
      if (enforce_truncate) {
        cap <- quantile(empirical_wet[[m]], 0.99, na.rm=TRUE)
        x <- min(x, log(cap + 1e-5))
      }
      prcp[i] <- exp(x)
    }
  }
  
  # 8. (Optional) AR(1) on log-intensity
  if (!is.null(ar_phi)) {
    logp <- log(prcp + 1e-5)
    for (i in 2:n) {
      if (wet[i]==1 && wet[i-1]==1) {
        logp[i] <- ar_phi * logp[i-1] + sqrt(1-ar_phi^2)*rnorm(1,0,sd(logp, na.rm=TRUE))
        prcp[i] <- exp(logp[i])
      }
    }
  }
  
  prcp <- pmax(pmin(prcp, prcp_range[2]), prcp_range[1])
  
  # 9. Simulate temperature (simple monthly normal)
  temp <- numeric(n)
  temp_params <- t(sapply(1:12, function(m) {
    temps <- df$TEMP[df$month==m]
    c(mu = mean(temps, na.rm=TRUE), sd = sd(temps, na.rm=TRUE))
  }))
  rownames(temp_params) <- month.abb
  
  for (i in seq_len(n)) {
    m <- month.abb[mon[i]]
    temp[i] <- rnorm(1, 
                     mean = temp_params[m,"mu"], 
                     sd   = temp_params[m,"sd"])
  }
  temp <- pmax(pmin(temp, temp_range[2]), temp_range[1])
  
  # 10. Return xts
  out <- xts::xts(cbind(PRECIP = prcp, TEMP = temp), order.by = dates)
  return(out)
}
