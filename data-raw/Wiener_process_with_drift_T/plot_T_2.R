# $Id: $
# Parameter settings
xd <- log(47 / 10)     # Log-distance from latest size to extinction threshold
mu <- 0.0023556        # Drift (mean growth rate)
sigma <- sqrt(0.01087) # Environmental variability standard deviation
time_limit <- 600      # Observation period

# Theoretical extinction time density function g(t | xd)
g_t <- function(t) {
  xd / sqrt(2 * pi * sigma^2 * t^3) *
    exp(- (xd + mu * t)^2 / (2 * sigma^2 * t))
}

ext_times <- read.csv("extinction_times.csv")

# Create histogram (not plotted, density scale)
histdata <- hist(ext_times$extinction_time, breaks = 50, plot = FALSE)

# Calculate total area of histogram (sum of density × bin width)
hist_area <- sum(histdata$density * diff(histdata$breaks))

# Calculate integral of theoretical density from 0 to time_limit (partial area)
p_ext_time_limit <- integrate(g_t, lower = 0, upper = time_limit)$value

# Scale factor: ratio of histogram area to theoretical density area
scale_factor <- hist_area / p_ext_time_limit

# Scaled theoretical density function
g_t_scaled <- function(t) {
  g_t(t) * scale_factor
}

# Plot histogram (density)
hist(ext_times$extinction_time, breaks = 50, col = "lightblue",
     main = "Histogram of Extinction Times",
     xlab = "Extinction time (years)",
     freq = FALSE)

# Overlay scaled theoretical density curve
curve(g_t_scaled, from = 0.01, to = time_limit,
      col = "red", lwd = 2, add = TRUE)
