ext_times <- read.csv("extinction_times.csv")
hist(ext_times$extinction_time, breaks = 50, col = "lightblue",
     main = "Histogram of Extinction Times",
     xlab = "Extinction time (years)")
