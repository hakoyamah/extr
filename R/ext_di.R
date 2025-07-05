# $Id: $
# Copyright (c) 2013-2025 Hiroshi Hakoyama
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#' @title Extinction Risk Estimation for a Density-Independent Population Model
#'
#' @description Estimates demographic parameters and extinction probabilities
#' over a specified time period based on a density-independent population model
#' formulated as the Wiener process with drift to approximate population
#' dynamics under environmental fluctuations. This function uses a time series
#' of population sizes to estimate extinction risk under this model and applies
#' the \eqn{w}-\eqn{z} method (Hakoyama, 2025), which provides more accurate
#' confidence intervals than conventional approaches such as the delta method.
#'
#' @param dat a data frame with two columns representing time and population
#' size, respectively. The column names are not restricted, but the columns must
#' be ordered by time first, followed by population size. Time intervals do not
#' need to be equally spaced.
#' @param ne numeric: Extinction threshold, \eqn{n_e \geq 1}. The default is 1
#' individual.
#' @param th numeric: Time horizon, \eqn{t_h}. The risk of extinction \eqn{G}
#' is defined as the probability that a population's size will drop below a
#' threshold \eqn{n_e} within a given time period \eqn{t_h}. The default is 100.
#' @param alpha numeric: Significance level, i.e., \eqn{\alpha} =
#' 1 - (confidence level). The default is 0.05.
#' @param unit character: The unit of time used in the dataset (e.g., "years",
#' "days", "generations"). The default is "years". This determines the time unit
#' for the growth rate \eqn{\mu}, environmental variance \eqn{\sigma^2}, and
#' extinction probability \eqn{G}.
#' @param qq_plot logical: If \code{TRUE}, shows a QQ-plot for
#' \eqn{\left(\log(n_{t_i} / n_{t_{i-1}}) - \hat{\mu} \tau_i\right) /
#' \left(\hat{\sigma} \sqrt{\tau_i}\right) \sim N(0, 1)}, where \eqn{\hat{\mu}}
#' and \eqn{\hat{\sigma}} denote the estimated parameters. Here,
#' \eqn{(n_{t_0}, n_{t_1}, \dots, n_{t_q})} represents the time series of
#' population sizes taken from the second column of \code{dat}. The default is
#' \code{FALSE}.
#' @param formatted logical: If \code{TRUE}, returns a user-friendly, formatted
#' output. Otherwise, returns a raw list. The default is \code{TRUE}.
#'
#' @details The \code{ext_di()} function estimates demographic parameters and
#' extinction probabilities over a specified time period. It is based on the
#' Wiener drift process model, which is one of the simplest frameworks for
#' assessing extinction risk (Lande & Orzack, 1988). Population dynamics are
#' modeled by the stochastic differential equation:
#' \deqn{
#' dX = \mu \, dt + \sigma \, dW,
#' }
#' where \eqn{X(t) = \log{N(t)}} is the logarithm of population size at time
#' \eqn{t}, \eqn{\mu} is the growth rate, \eqn{\sigma^2} is the environmental
#' variance, and \eqn{W(t)} is the Wiener process. The diffusion term on the
#' right-hand side is interpreted in the Ito sense.
#'
#' The function implements a novel method (the \eqn{w}–\eqn{z} method, Hakoyama,
#' 2025) that provides more accurate confidence intervals for extinction risk
#' than conventional methods such as the delta method.
#'
#' The estimation process involves three main steps:
#' 1. Estimating demographic parameters such as population growth rate and
#'    environmental variance (Dennis et al., 1991).
#' 2. Estimating extinction probability using analytical formulas based on the
#'    Wiener drift process developed by Lande and Orzack (1988), who modeled
#'    extinction risk under environmental fluctuations.
#' 3. Calculating confidence intervals for these estimates using a novel method
#'    (Hakoyama, 2025).
#'
#' For step 3, the function implements a novel method (the \eqn{w}–\eqn{z}
#' method, Hakoyama, 2025) that introduces two transformed parameters,
#' \eqn{w} and \eqn{z}, as functions of \eqn{\mu}, \eqn{\sigma}, the time
#' horizon \eqn{t_h}, and the log-distance to the extinction threshold
#' \eqn{x_d}. These transformations simplify the extinction probability
#' formula \eqn{G(w, z)} and yield confidence intervals with improved coverage
#' properties compared to existing methods.
#'
#' Key features:
#' - Accepts a time series of population sizes
#'   (\eqn{n_{t_0}, n_{t_1}, \dots, n_{t_q}}).
#' - Allows for irregular measurement intervals
#'   (\eqn{\tau_i = t_i - t_{i-1}}), permitting missing values in the time
#'   series (Dennis et al., 1991).
#' - Returns maximum likelihood estimates (MLEs) and confidence intervals for
#'   the growth rate, environmental variance, and extinction probability using
#'   the \eqn{w}–\eqn{z} method (Hakoyama, 2025).
#' - Can generate QQ-plots to assess model assumptions.
#'
#' @return If \code{formatted = TRUE}, returns the following components in a
#' formatted structure. If \code{formatted = FALSE}, returns a list:
#'
#' Estimates:
#'   - \code{Growth.rate}: MLE of growth rate, \eqn{\hat{\mu}}.
#'   - \code{Variance}: MLE of environmental variance, \eqn{\hat{\sigma}^2}.
#'   - \code{Unbiased.variance}: Unbiased estimate of environmental variance.
#'   - \code{AIC}: Akaike Information Criterion for the distribution of
#'     population size \eqn{N}.
#'   - \code{Extinction.probability}: MLE of the extinction probability within
#'     the time horizon \eqn{t_h}, denoted \eqn{G}.
#'   - \code{lower_cl_mu}: Lower \eqn{(1 - \alpha)\%} confidence limit of
#'     \eqn{\mu}.
#'   - \code{upper_cl_mu}: Upper \eqn{(1 - \alpha)\%} confidence limit of
#'     \eqn{\mu}.
#'   - \code{lower_cl_s}: Lower \eqn{(1 - \alpha)\%} confidence limit of
#'     \eqn{\sigma^2}.
#'   - \code{upper_cl_s}: Upper \eqn{(1 - \alpha)\%} confidence limit of
#'     \eqn{\sigma^2}.
#'   - \code{lower_cl_p}: Lower \eqn{(1 - \alpha)\%} confidence limit of
#'     extinction probability \eqn{G}.
#'   - \code{upper_cl_p}: Upper \eqn{(1 - \alpha)\%} confidence limit of
#'     extinction probability \eqn{G}.
#'
#' Data summary:
#'   - \code{nq}: Latest observed population size, \eqn{n_q}.
#'   - \code{xd}: Log-distance from latest size to extinction threshold,
#'     \eqn{x_d = \log(n_q / n_e)}.
#'   - \code{sample.size}: Sample size (number of observations), \eqn{q + 1}.
#'
#' Input parameters:
#'   - \code{unit}, \code{ne}, \code{th}, and \code{alpha}.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @references
#'
#' Lande, R. and Orzack, S.H. (1988) Extinction dynamics of age-structured
#' populations in a fluctuating environment. *Proceedings of the National
#' Academy of Sciences*, 85(19), 7418–7421.
#'
#' Dennis, B., Munholland, P.L., and Scott, J.M. (1991) Estimation of growth
#' and extinction parameters for endangered species. *Ecological Monographs*,
#' 61, 115–143.
#'
#' Hakoyama, H. (2025) Confidence interval for extinction risk: Revisiting the
#' misconception that the interval is too imprecise to be useful. Preprint,
#' available at \url{https://###/###.pdf} (accessed 2025-07-##).
#'
#' @examples
#' # Yellowstone grizzly bears (from Dennis et al., 1991)
#' dat <- data.frame(Time = c(1959, 1960, 1961, 1962, 1963, 1964, 1965,
#' 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976,
#' 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987),
#' Population = c(44, 47, 46, 44, 46, 45, 46, 40, 39, 39, 42, 44, 41, 40,
#' 33, 36, 34, 39, 35, 34, 38, 36, 37, 41, 39, 51, 47, 57, 47))
#'
#' # Probability of decline to 1 individual within 100 years
#' ext_di(dat, th = 100)
#'
#' # Probability of decline to 10 individuals within 100 years
#' ext_di(dat, th = 100, ne = 10)
#'
#' # With QQ-plot
#' ext_di(dat, th = 100, ne = 10, qq_plot = TRUE)
#'
#' @keywords models survival time-series methods
#'
#' @importFrom stats complete.cases qqnorm qqline qt qchisq
#'
#' @export
#'
#' @seealso [`extr::w_statistic`], [`extr::z_statistic`],
#' [`extr::extinction_probability_wz_di`],
#' [`extr::confidence_interval_wz_di`], [`extr::print.ext_di`]
#'
ext_di <- function(dat, ne = 1, th = 100, alpha = 0.05, unit = "years",
                   qq_plot = FALSE, formatted = TRUE) {
  if (ne < 1) stop("Extinction threshold 'ne' must be >= 1.")
  if (!is.data.frame(dat) || ncol(dat) != 2) {
    stop(paste0("Input 'dat' must be a data frame with two columns:",
                "time first, then population size."))
  }
  yr <- dat[, 1]
  n  <- dat[, 2]
  if (!is.numeric(yr) || !is.numeric(n)) {
    stop("Both columns of 'dat' must be numeric.")
  }
  if (any(yr < 0)) {
    stop("'time' column contains negative values, which is not allowed.")
  }
  if (any(n < 1)) {
    stop(paste0("Population size column contains values less than 1, ",
                "which is not allowed."))
  }
  if (length(yr) < 2) stop("At least two observations are required.")
  complete <- complete.cases(yr, n)
  yr <- yr[complete]
  n  <-  n[complete]
  ti <- yr - yr[1]
  tau <- diff(ti)
  delta_log_n <- diff(log(n))
  qq <- length(yr) - 1
  tq <- ti[length(ti)]
  nq <- n[length(n)]
  mu <- log(nq / n[1]) / tq
  s <- sum((delta_log_n - mu * tau)^2 / tau) / qq
  us <- qq * s / (qq - 1)
  xd <- log(nq / ne)
  lnl <- - sum(log(n[-1] * sqrt(2 * tau * pi))) - (qq / 2) * log(s) -
    sum((1 / tau) * (delta_log_n - mu * tau)^2) / (2 * s)
  aic <- - 2 * lnl + 2 * 2

  if (qq_plot) {
    z_vals <- (delta_log_n - mu * tau) / (sqrt(s) * sqrt(tau))
    qqnorm(z_vals)
    qqline(z_vals)
  }

  lower_cl_mu <- mu - qt(1 - alpha / 2, qq - 1) * sqrt(us / tq)
  upper_cl_mu <- mu + qt(1 - alpha / 2, qq - 1) * sqrt(us / tq)
  lower_cl_s <- qq * s / qchisq(1 - alpha / 2, qq - 1)
  upper_cl_s <- qq * s / qchisq(alpha / 2, qq - 1)

  ww <- w_statistic(mu, xd, s, th)
  zz <- z_statistic(mu, xd, s, th)
  pp <- extinction_probability_wz_di(ww, zz)

  cl_p <- confidence_interval_wz_di(mu, xd, s, th, tq, qq, alpha)

  lower_cl_p <- cl_p[[1]]
  upper_cl_p <- cl_p[[2]]

  results <- list(ne = ne, th = th,
                  alpha = alpha,
                  unit = unit,
                  aic = aic,
                  sample.size = qq + 1,
                  nq = nq,
                  xd = xd,
                  Growth.rate = mu,
                  lower_cl_mu = lower_cl_mu,
                  upper_cl_mu = upper_cl_mu,
                  Variance = s,
                  lower_cl_s = lower_cl_s,
                  upper_cl_s = upper_cl_s,
                  Unbiased.variance = us,
                  Extinction.probability = pp,
                  lower_cl_p = lower_cl_p,
                  upper_cl_p = upper_cl_p)
  if (formatted == TRUE) {
    class(results) <- "ext_di"
  }
  return(results)
}
