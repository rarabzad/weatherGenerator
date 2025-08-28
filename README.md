
# Weather Generator

This repository provides an R function to generate synthetic daily weather time series (precipitation and temperature) based on historical records, along with an example workflow that:

1. Loads the generator function and sample forcing data  
2. Converts data into `xts` objects and runs the simulation  
3. Computes daily & monthly climatologies and heatwave statistics  
4. Produces comparison plots and summary tables

---

## 1. Installation

You can install directly from GitHub using **devtools**, or simply source the function file:

```r
# Install from GitHub via devtools
source("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/generate_weather.R")
````

### Package dependencies

The function will attempt to **load packages silently**, and if any are missing, it will automatically install them from CRAN:

* Required: `xts`, `dplyr`, `lubridate`
* Optional (recommended): `moments` (used for AR(1) estimation, but function still works without it)

You can also install them manually:

```r
install.packages(c("xts", "dplyr", "lubridate", "moments"))
```

---

## 2. Function Documentation

**Generate Synthetic Daily Weather**

**Description**
Simulates a multi-year daily time series of precipitation and temperature by fitting a monthly two-state (wet/dry) Markov chain to historical data, sampling wet-day precipitation either by bootstrap or a log-normal parametric model (with optional AR(1) on log-intensity), and drawing daily temperature from a monthly normal distribution.

**Arguments**

* `precip_xts`: An `xts` series of daily precipitation (mm), indexed by Date.
* `temp_xts`: An `xts` series of daily temperature (°C), indexed by Date.
* `nyears` (default = 100): Number of years to simulate.
* `prcp_range` (default = c(0, Inf)): Clamp simulated precipitation.
* `temp_range` (default = c(-60, 60)): Clamp simulated temperature.
* `use_bootstrap` (default = TRUE): If TRUE, sample wet-day intensities from historical pool; else use parametric log-normal with seasonal adjustment.
* `roll_window` (default = 31): Window size for parametric smoothing.
* `prior_counts` (default = 1): Additive smoothing when estimating transition probabilities P01 and P11.
* `enforce_truncate` (default = TRUE): If TRUE, cap parametric draws at 99th percentile of historical intensities.
* `ar_phi` (default = NULL): If not supplied, estimated automatically. Applies AR(1) persistence to log-intensity of consecutive wet days.
* `bootstrap_params`: Reserved (not yet implemented).

**Mechanics**

1. **Package loading**: Ensures required libraries are loaded or installed.
2. **Data preparation**: Merge precipitation & temperature, extract month/day-of-year tags.
3. **AR(1) estimation**: If `ar_phi = NULL`, estimate wet-day log-intensity autocorrelation.
4. **Seasonal scaling**: Compute daily smoothed mean/median/90p precipitation intensities by day-of-year, fill missing days by interpolation.
5. **Empirical wet-day pool**: Store monthly wet-day intensities for bootstrap sampling.
6. **Markov chain**: Estimate monthly transition probabilities (P01, P11) with pseudo-counts.
7. **Log-normal parameters**: Fit monthly log-intensity mean/SD.
8. **Synthetic timeline**: Generate daily date sequence over `nyears`, including leap years.
9. **Wet/dry simulation**: Use monthly Markov chain to generate binary wet sequence.
10. **Precipitation simulation**:

    * **Bootstrap**: Randomly sample from empirical wet pool, scaled by daily seasonal factor.
    * **Parametric**: Draw from Normal(μ,σ) in log-space, adjust by seasonal factor, truncate if enabled.
11. **AR(1) persistence**: Optionally apply AR(1) to log-intensities across wet days.
12. **Temperature simulation**: Draw daily values from monthly Normal(μ,σ).
13. **Clamp values**: Enforce user-defined ranges for precipitation and temperature.

**Return Value**
An `xts` object with columns:

* `PRECIP` — simulated daily precipitation (mm)
* `TEMP` — simulated daily average temperature (°C)

---

## 3. Example Workflow

## 3. Example Workflow

Below is a step‑by‑step walkthrough of how to run the generator, compute summary statistics, and produce comparison plots.

1. **Load data & function**
   First, source the generator and load your historical forcing data into `xts` time series.

   ```r
   library(xts)
   library(zoo)
   library(lubridate)
   library(dplyr)
   library(ggplot2)
   library(knitr)
   
   source("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/generate_weather.R")

   data <- read.csv("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/ForcingFunctions.csv")

   precip_xts <- xts(data$precipitation..mm., order.by = as.Date(data$date))
   temp_xts   <- xts(data$temp..C.,           order.by = as.Date(data$date))
   colnames(precip_xts) <- "PRECIP"
   colnames(temp_xts)   <- "TEMP"
   ```

2. **Generate synthetic series**
   Call `generate_weather()` with default settings (100 years, bootstrap).

   ```r
   synthetic_xts <- generate_weather(precip_xts, temp_xts)
   ```

3. **Convert to data frames & add Year/DOY**
   Prepare observed and synthetic data frames for analysis.

   ```r
   obs_df <- data.frame(
     date   = index(precip_xts),
     PRECIP = coredata(precip_xts),
     TEMP   = coredata(temp_xts)
   ) %>%
     mutate(year = year(date), doy = yday(date))

   syn_df <- data.frame(
     date  = index(synthetic_xts),
     PRECIP = coredata(synthetic_xts[, "PRECIP"]),
     TEMP   = coredata(synthetic_xts[, "TEMP"])
   ) %>%
     mutate(year = year(date), doy = yday(date))
   ```

4. **Define summary routines**
   Functions for daily climatology, monthly stats, and heatwave counts.

   ```r
   daily_means <- function(df) {
     df %>%
       group_by(doy) %>%
       summarise(
         mean_temp   = mean(TEMP,   na.rm = TRUE),
         mean_precip = mean(PRECIP, na.rm = TRUE)
       )
   }

   monthly_stats <- function(df) {
     df %>%
       mutate(month = month(date, label = TRUE)) %>%
       group_by(month) %>%
       summarise(
         mean_temp   = mean(TEMP,   na.rm = TRUE),
         sd_temp     = sd(TEMP,     na.rm = TRUE),
         mean_precip = mean(PRECIP, na.rm = TRUE),
         sd_precip   = sd(PRECIP,   na.rm = TRUE)
       )
   }

   count_heatwaves <- function(df, threshold = 30, duration = 3) {
     df %>%
       mutate(hot = TEMP > threshold) %>%
       group_by(year) %>%
       summarise(
         heatwaves = sum(rle(hot)$lengths[rle(hot)$values] >= duration)
       )
   }
   ```

5. **Compute summaries**

   ```r
   obs_daily   <- daily_means(obs_df)
   syn_daily   <- daily_means(syn_df)
   obs_monthly <- monthly_stats(obs_df)
   syn_monthly <- monthly_stats(syn_df)
   obs_hw      <- count_heatwaves(obs_df, threshold = 19, duration = 3)
   syn_hw      <- count_heatwaves(syn_df, threshold = 19, duration = 3)
   ```

## 4. Visualization

1. **Daily climatology** (mean temp & precip by day of year)
2. **Monthly statistics** (mean and SD of temp and precip)
3. **Heatwave frequency** (boxplot of yearly counts)

---

### 📈 1. Daily Climatology: Mean Temp & Precip by DOY

```r
library(ggplot2)
library(dplyr)

# Merge for plotting
daily_combined <- bind_rows(
  obs_daily %>% mutate(Source = "Observed"),
  syn_daily %>% mutate(Source = "Synthetic")
)

# Temperature
ggplot(daily_combined, aes(x = doy, y = mean_temp, color = Source)) +
  geom_line() +
  labs(title = "Daily Mean Temperature", x = "Day of Year", y = "Temperature (°C)") +
  theme_minimal()

# Precipitation
ggplot(daily_combined, aes(x = doy, y = mean_precip, color = Source)) +
  geom_line() +
  labs(title = "Daily Mean Precipitation", x = "Day of Year", y = "Precipitation (mm)") +
  theme_minimal()
```

---

### 📊 2. Monthly Mean & SD of Temp and Precip

```r
# Merge monthly stats
monthly_combined <- bind_rows(
  obs_monthly %>% mutate(Source = "Observed"),
  syn_monthly %>% mutate(Source = "Synthetic")
)

# Mean temperature
ggplot(monthly_combined, aes(x = month, y = mean_temp, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Mean Temperature", x = "Month", y = "°C") +
  theme_minimal()

# Mean precipitation
ggplot(monthly_combined, aes(x = month, y = mean_precip, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Mean Precipitation", x = "Month", y = "mm") +
  theme_minimal()

# SD temperature
ggplot(monthly_combined, aes(x = month, y = sd_temp, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Temperature SD", x = "Month", y = "°C") +
  theme_minimal()

# SD precipitation
ggplot(monthly_combined, aes(x = month, y = sd_precip, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Precipitation SD", x = "Month", y = "mm") +
  theme_minimal()
```

---

### 🔥 3. Heatwave Frequency (Boxplot of Annual Counts)

```r
hw_combined <- bind_rows(
  obs_hw %>% mutate(Source = "Observed"),
  syn_hw %>% mutate(Source = "Synthetic")
)

ggplot(hw_combined, aes(x = Source, y = heatwaves, fill = Source)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Annual Heatwave Counts (≥19°C, ≥3 Days)", y = "Count", x = "") +
  theme_minimal()
```

---

### 🌡️ 4. Temperature Density Plot

```r
# Combine observed and synthetic temperatures
temp_density_df <- bind_rows(
  obs_df %>% select(TEMP) %>% mutate(Source = "Observed"),
  syn_df %>% select(TEMP) %>% mutate(Source = "Synthetic")
)

ggplot(temp_density_df, aes(x = TEMP, fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Temperature Density", x = "Temperature (°C)", y = "Density") +
  theme_minimal()
```

---

### 🌧️ 5. Precipitation Density Plot

#### A. **Raw Precipitation**

```r
precip_density_df <- bind_rows(
  obs_df %>% select(PRECIP) %>% mutate(Source = "Observed"),
  syn_df %>% select(PRECIP) %>% mutate(Source = "Synthetic")
)

ggplot(precip_density_df, aes(x = PRECIP, fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Precipitation Density", x = "Precipitation (mm)", y = "Density") +
  xlim(0, quantile(precip_density_df$PRECIP, 0.99, na.rm = TRUE)) +
  theme_minimal()
```

#### B. **Log-Transformed Precipitation**

```r
ggplot(precip_density_df, aes(x = log(PRECIP + 1e-5), fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Log-Transformed Precipitation Density", x = "log(Precipitation + 1e-5)", y = "Density") +
  theme_minimal()
```

---

### License

MIT © Rezgar Arabzadeh





# Weather Generator

This repository provides an R function to generate synthetic daily weather time series (precipitation and temperature) based on historical records, along with an example workflow that:

1. Loads the generator function and sample forcing data  
2. Converts data into `xts` objects and runs the simulation  
3. Computes daily & monthly climatologies and heatwave statistics  
4. Produces comparison plots and summary tables

---

## 1. Installation

You can install directly from GitHub using **devtools**, or simply source the function file:

```r
# Install from GitHub via devtools
source("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/generate_weather.R")
````

Required packages:

```r
install.packages(c("xts", "zoo", "lubridate", "dplyr", "ggplot2", "knitr"))
```

---

## 2. Function Documentation

**Generate Synthetic Daily Weather**

**Description**
Simulate a multi‑year daily time series of precipitation and temperature by fitting a monthly two‑state (wet/dry) Markov chain to historical data, sampling wet‑day precipitation either by bootstrap or a log‑normal parametric model (with optional AR(1) on log‑intensity), and drawing daily temperature from a monthly normal distribution.

**Arguments**

* `precip_xts`: An `xts` series of daily precipitation (mm), indexed by Date.
* `temp_xts`: An `xts` series of daily temperature (°C), indexed by Date.
* `nyears` (default = 100): Number of years to simulate.
* `prcp_range` (default = c(0, Inf)): Numeric vector of length 2; clamps simulated precipitation to \[min, max].
* `temp_range` (default = c(-60, 60)): Numeric vector of length 2; clamps simulated temperature to \[min, max].
* `use_bootstrap` (default = TRUE): If TRUE, sample precipitation intensities by drawing randomly from historical wet‑day values; if FALSE, draw from a log‑normal distribution fitted with a rolling window.
* `roll_window` (default = 31): Window size (days) for computing rolling‑window mean and SD of log‑intensity when `use_bootstrap = FALSE`.
* `prior_counts` (default = 1): Pseudo‑count (additive smoothing) for estimating monthly Markov probabilities P01 and P11.
* `enforce_truncate` (default = TRUE): If TRUE, cap log‑normal draws at the 99th percentile of historical wet‑day log‑intensities.
* `ar_phi` (default = NULL): Optional numeric; if provided, applies AR(1) on the log(precip+ε) series across consecutive wet days with coefficient φ.

**Mechanics**

1. **Data merge & month tagging**: Merge precipitation and temperature `xts` objects and extract each day’s month.
2. **Empirical wet‑day pool**: For each calendar month, collect all historical positive precipitation values for bootstrap sampling.
3. **Markov chain parameters**: Compute monthly probabilities P01 (dry→wet) and P11 (wet→wet) with add‑1 smoothing:

   $$
     P01 = \frac{n_{0\to1} + \alpha}{n_{0} + 2\,\alpha},\quad
     P11 = \frac{n_{1\to1} + \alpha}{n_{1} + 2\,\alpha}.
   $$
4. **Rolling log‑normal fit** (when `use_bootstrap = FALSE`): Calculate centered rolling‑window mean and SD of log(precip + 1e–5), then sample new values from Normal(μ, σ) and back‑transform.
5. **Synthetic timeline**: Build a sequence of Dates spanning `nyears` years (including leap days) starting at “2000‑01‑01”.
6. **Simulate wet/dry**: Initialize the first day’s wet probability from historical frequency, then use the monthly Markov chain to generate a Bernoulli sequence.
7. **Sample intensities**:
   * **Bootstrap**: Draw a random historic wet‑day precipitation for each wet day.
   * **Parametric**: Sample log‑intensity, optionally truncate at the 99th percentile, then exponentiate.
8. **AR(1) on log‑intensity**: If `ar_phi` is set, apply log pᵢ = ϕ · log pᵢ₋₁ + √(1 - ϕ²) · εᵢ across consecutive wet days.
9. **Clamp values**: Ensure simulated `PRECIP` and `TEMP` lie within user‑specified ranges.
11. **Temperature draw**: For each day, sample from Normal(μ\_m, σ\_m) estimated from historical temperatures in month *m*.

**Return Value**
An `xts` object with columns:

* `PRECIP` (simulated daily precipitation in mm)
* `TEMP_DAILY_AVE` (simulated daily average temperature in °C)

---

## 3. Example Workflow

Below is a step‑by‑step walkthrough of how to run the generator, compute summary statistics, and produce comparison plots.

1. **Load data & function**
   First, source the generator and load your historical forcing data into `xts` time series.

   ```r
   library(xts)
   library(zoo)
   library(lubridate)
   library(dplyr)
   library(ggplot2)
   library(knitr)
   
   source("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/generate_weather.R")

   data <- read.csv("https://github.com/rarabzad/weatherGenerator/raw/refs/heads/main/ForcingFunctions.csv")

   precip_xts <- xts(data$precipitation..mm., order.by = as.Date(data$date))
   temp_xts   <- xts(data$temp..C.,           order.by = as.Date(data$date))
   colnames(precip_xts) <- "PRECIP"
   colnames(temp_xts)   <- "TEMP"
   ```

2. **Generate synthetic series**
   Call `generate_weather()` with default settings (100 years, bootstrap).

   ```r
   synthetic_xts <- generate_weather(precip_xts, temp_xts)
   ```

3. **Convert to data frames & add Year/DOY**
   Prepare observed and synthetic data frames for analysis.

   ```r
   obs_df <- data.frame(
     date   = index(precip_xts),
     PRECIP = coredata(precip_xts),
     TEMP   = coredata(temp_xts)
   ) %>%
     mutate(year = year(date), doy = yday(date))

   syn_df <- data.frame(
     date  = index(synthetic_xts),
     PRECIP = coredata(synthetic_xts[, "PRECIP"]),
     TEMP   = coredata(synthetic_xts[, "TEMP"])
   ) %>%
     mutate(year = year(date), doy = yday(date))
   ```

4. **Define summary routines**
   Functions for daily climatology, monthly stats, and heatwave counts.

   ```r
   daily_means <- function(df) {
     df %>%
       group_by(doy) %>%
       summarise(
         mean_temp   = mean(TEMP,   na.rm = TRUE),
         mean_precip = mean(PRECIP, na.rm = TRUE)
       )
   }

   monthly_stats <- function(df) {
     df %>%
       mutate(month = month(date, label = TRUE)) %>%
       group_by(month) %>%
       summarise(
         mean_temp   = mean(TEMP,   na.rm = TRUE),
         sd_temp     = sd(TEMP,     na.rm = TRUE),
         mean_precip = mean(PRECIP, na.rm = TRUE),
         sd_precip   = sd(PRECIP,   na.rm = TRUE)
       )
   }

   count_heatwaves <- function(df, threshold = 30, duration = 3) {
     df %>%
       mutate(hot = TEMP > threshold) %>%
       group_by(year) %>%
       summarise(
         heatwaves = sum(rle(hot)$lengths[rle(hot)$values] >= duration)
       )
   }
   ```

5. **Compute summaries**

   ```r
   obs_daily   <- daily_means(obs_df)
   syn_daily   <- daily_means(syn_df)
   obs_monthly <- monthly_stats(obs_df)
   syn_monthly <- monthly_stats(syn_df)
   obs_hw      <- count_heatwaves(obs_df, threshold = 19, duration = 3)
   syn_hw      <- count_heatwaves(syn_df, threshold = 19, duration = 3)
   ```




## 4. Visualization

1. **Daily climatology** (mean temp & precip by day of year)
2. **Monthly statistics** (mean and SD of temp and precip)
3. **Heatwave frequency** (boxplot of yearly counts)

---

### 📈 1. Daily Climatology: Mean Temp & Precip by DOY

```r
library(ggplot2)
library(dplyr)

# Merge for plotting
daily_combined <- bind_rows(
  obs_daily %>% mutate(Source = "Observed"),
  syn_daily %>% mutate(Source = "Synthetic")
)

# Temperature
ggplot(daily_combined, aes(x = doy, y = mean_temp, color = Source)) +
  geom_line() +
  labs(title = "Daily Mean Temperature", x = "Day of Year", y = "Temperature (°C)") +
  theme_minimal()

# Precipitation
ggplot(daily_combined, aes(x = doy, y = mean_precip, color = Source)) +
  geom_line() +
  labs(title = "Daily Mean Precipitation", x = "Day of Year", y = "Precipitation (mm)") +
  theme_minimal()
```

---

### 📊 2. Monthly Mean & SD of Temp and Precip

```r
# Merge monthly stats
monthly_combined <- bind_rows(
  obs_monthly %>% mutate(Source = "Observed"),
  syn_monthly %>% mutate(Source = "Synthetic")
)

# Mean temperature
ggplot(monthly_combined, aes(x = month, y = mean_temp, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Mean Temperature", x = "Month", y = "°C") +
  theme_minimal()

# Mean precipitation
ggplot(monthly_combined, aes(x = month, y = mean_precip, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Mean Precipitation", x = "Month", y = "mm") +
  theme_minimal()

# SD temperature
ggplot(monthly_combined, aes(x = month, y = sd_temp, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Temperature SD", x = "Month", y = "°C") +
  theme_minimal()

# SD precipitation
ggplot(monthly_combined, aes(x = month, y = sd_precip, fill = Source)) +
  geom_col(position = "dodge") +
  labs(title = "Monthly Precipitation SD", x = "Month", y = "mm") +
  theme_minimal()
```

---

### 🔥 3. Heatwave Frequency (Boxplot of Annual Counts)

```r
hw_combined <- bind_rows(
  obs_hw %>% mutate(Source = "Observed"),
  syn_hw %>% mutate(Source = "Synthetic")
)

ggplot(hw_combined, aes(x = Source, y = heatwaves, fill = Source)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Annual Heatwave Counts (≥19°C, ≥3 Days)", y = "Count", x = "") +
  theme_minimal()
```








---

### 🌡️ 4. Temperature Density Plot

```r
# Combine observed and synthetic temperatures
temp_density_df <- bind_rows(
  obs_df %>% select(TEMP) %>% mutate(Source = "Observed"),
  syn_df %>% select(TEMP) %>% mutate(Source = "Synthetic")
)

ggplot(temp_density_df, aes(x = TEMP, fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Temperature Density", x = "Temperature (°C)", y = "Density") +
  theme_minimal()
```

---

### 🌧️ 5. Precipitation Density Plot

#### A. **Raw Precipitation**

```r
precip_density_df <- bind_rows(
  obs_df %>% select(PRECIP) %>% mutate(Source = "Observed"),
  syn_df %>% select(PRECIP) %>% mutate(Source = "Synthetic")
)

ggplot(precip_density_df, aes(x = PRECIP, fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Precipitation Density", x = "Precipitation (mm)", y = "Density") +
  xlim(0, quantile(precip_density_df$PRECIP, 0.99, na.rm = TRUE)) +
  theme_minimal()
```

#### B. **Log-Transformed Precipitation**

```r
ggplot(precip_density_df, aes(x = log(PRECIP + 1e-5), fill = Source)) +
  geom_density(alpha = 0.5) +
  labs(title = "Log-Transformed Precipitation Density", x = "log(Precipitation + 1e-5)", y = "Density") +
  theme_minimal()
```

---

---

## 5. Interpretation

Use these comparisons to evaluate how closely the synthetic series matches observed climatology, variability, and extreme‐event statistics. Adjust function arguments (e.g. `use_bootstrap = FALSE`, `ar_phi = 0.8`, `nyears = 50`) to test sensitivity.

---

### License

MIT © Rezgar Arabzadeh
