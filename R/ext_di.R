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
#' @title Extinction Risk Estimation for a Density-Independent Model
#'
#' @description
#' Estimates demographic parameters and extinction probability under a
#' density-independent (drifted Wiener) model. From a time series of
#' population sizes, it computes MLEs of growth rate and environmental
#' variance, then evaluates extinction risk over a horizon \eqn{t^{\ast}}.
#' Confidence intervals are constructed by the \eqn{w}-\eqn{z} method,
#' which achieve near-nominal coverage across the full parameter space.
#'
#' @param dat Data frame with two numeric columns: time (strictly increasing)
#'   and population size. Column names are not restricted.
#' @param ne Numeric. Extinction threshold \eqn{n_e \ge 1}. Default is 1.
#' @param th Numeric. Time horizon \eqn{t^{\ast} > 0}. Default is 100.
#' @param alpha Numeric. Significance level \eqn{\alpha \in (0,1)}.
#'   Default is 0.05.
#' @param unit Character. Unit of time (e.g., "years", "days", "generations").
#'   Default is "years".
#' @param qq_plot Logical. If \code{TRUE}, draws a QQ-plot of standardized
#'   increments to check model assumptions. Default is \code{FALSE}.
#' @param formatted Logical. If \code{TRUE}, returns an \code{"ext_di"} object;
#'   otherwise returns a raw list. Default is \code{TRUE}.
#' @param digits Integer. Preferred significant digits for printing. Affects
#'   display only. Default is \code{getOption("extr.digits", 5)}.
#'
#' @details
#' Population dynamics follow \deqn{dX=\mu\,dt+\sigma\,dW,} where
#' \eqn{X(t)=\log N(t)}, \eqn{\mu} is the growth rate, \eqn{\sigma^2} the
#' environmental variance, and \eqn{W} a Wiener process. Extinction risk is
#' \deqn{G=\Pr[T\le t^{\ast}\mid N(0)=n_0,n_e,\mu,\sigma],} the probability the
#' population falls below \eqn{n_e} within \eqn{t^{\ast}}. Irregular intervals
#' are allowed.
#'
#' The function:
#' \enumerate{
#'   \item estimates \eqn{\mu} and \eqn{\sigma^2} (Dennis et al., 1991),
#'   \item computes extinction probability \eqn{G(w,z)}
#'         (Lande and Orzack, 1988),
#'   \item constructs confidence intervals for \eqn{G} using the
#'         \eqn{w}-\eqn{z} method (Hakoyama, 2025).
#' }
#'
#' Numerical range. Probabilities are evaluated on \eqn{G}, \eqn{\log G}, and
#' \eqn{\log(1-G)} scales. The log-scale removes the \eqn{\approx
#' 4.94\times10^{-324}} lower bound of linear doubles and extends the safe
#' range down to \code{exp(-DBL_MAX)} (kept symbolically), avoiding
#' underflow/cancellation.
#'
#' @return A list (class \code{"ext_di"} if \code{formatted=TRUE}) with:
#' \itemize{
#'   \item \code{Growth.rate}, \code{Variance}, \code{Unbiased.variance};
#'   \item \code{AIC};
#'   \item \code{Extinction.probability} with confidence limits;
#'   \item data summary (\code{nq}, \code{xd}, \code{sample.size});
#'   \item input parameters (\code{unit}, \code{ne}, \code{th}, \code{alpha}).
#' }
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @references
#' Lande, R. and Orzack, S.H. (1988) Extinction dynamics of age-structured
#' populations in a fluctuating environment. *Proceedings of the National
#' Academy of Sciences*, 85(19), 7418–7421.
#'
#' Dennis, B., Munholland, P.L., and Scott, J.M. (1991) Estimation of growth
#' and extinction parameters for endangered species. *Ecological Monographs*,
#' 61, 115–143.
#'
#' Hakoyama, H. (2025) Confidence intervals for extinction risk: validating
#' population viability analysis with limited data.
#' Preprint, doi:10.48550/arXiv.2509.09965
#'
#' @examples
#' # Example from Dennis et al. (1991), Yellowstone grizzly bears
#' dat <- data.frame(Time = 1959:1987,
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
#' ext_di(dat, th = 100, qq_plot = TRUE)
#'
#' # Change digits
#' ext_di(dat, th = 100, ne = 10, digits = 9)
#'
#' @keywords models survival time-series methods
#'
#' @importFrom stats complete.cases qqnorm qqline qt qchisq
#'
#' @export
#'
#' @seealso [`extr::statistics_di`],
#' [`extr::extinction_probability_di`],
#' [`extr::confidence_interval_wz_di`], [`extr::print.ext_di`]
#'

ext_di <- function(dat, ne = 1, th = 100, alpha = 0.05, unit = "years",
                   qq_plot = FALSE, formatted = TRUE,
                   digits = getOption("extr.digits", 5L)) {
  if (ne < 1) stop("Extinction threshold 'ne' must be >= 1.")
  if (th <= 0) stop("Time horizon 'th' must be > 0.")
  if (!(is.numeric(alpha) && alpha > 0 && alpha < 1)) {
    stop("'alpha' must be in (0, 1).")
  }
  if (!is.data.frame(dat) || ncol(dat) != 2) {
    stop(paste0("Input 'dat' must be a data frame with two columns: ",
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
  if (length(yr) < 3) {
    stop("At least three observations are required to compute CIs (q+1 >= 3).")
  }

  complete <- complete.cases(yr, n)
  yr <- yr[complete]
  n  <-  n[complete]
  ti <- yr - yr[1]
  tau <- diff(ti)
  if (any(tau <= 0)) {
    stop("Time values must be strictly increasing (no ties or decreases).")
  }
  delta_log_n <- diff(log(n))
  qq <- length(yr) - 1
  tq <- ti[length(ti)]
  nq <- n[length(n)]
  mu <- log(nq / n[1]) / tq
  s <- sum((delta_log_n - mu * tau)^2 / tau) / qq
  if (s == 0) {
    stop(paste0(
      "Estimated variance s is zero; check data ",
      "(perfectly linear log-series)."
    ))
  }
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
  if (ww + zz <= 0) {
    stop("Invalid input: require n0 = nq > ne")
  }

  linear_g  <- ext_prob_di(ww, zz)
  log_g  <- log_ext_prob_di(ww, zz)
  log_q  <- log_ext_comp_di(ww, zz)

  ci_linear_g <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = ext_prob_di
  )
  ci_log_g <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = log_ext_prob_di
  )
  ci_log_q <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = log_ext_comp_di
  )

  repr_point <- repr_mode(linear_g)
  repr_lower <- repr_mode(ci_linear_g[[1]])
  repr_upper <- repr_mode(ci_linear_g[[2]])

  results <- list(ne = ne, th = th,
                  alpha = alpha,
                  unit = unit,
                  aic = aic,
                  tq = tq,
                  qq = qq,
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
                  linear_g = linear_g,
                  log_g = log_g,
                  log_q = log_q,
                  ci_linear_g = ci_linear_g,
                  ci_log_g = ci_log_g,
                  ci_log_q = ci_log_q,
                  repr_point = repr_point,
                  repr_lower = repr_lower,
                  repr_upper = repr_upper,
                  digits = as.integer(digits))
  if (formatted) {
    class(results) <- "ext_di"
  }
  return(results)
}
