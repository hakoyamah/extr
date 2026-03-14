# $Id: $
# Copyright (c) 2026 Hiroshi Hakoyama
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
#' @title OEAR Diffusion-Scale Estimator
#'
#' @description
#' Compute the OEAR (observation-error-and-autocovariance-robust)
#' diffusion-scale estimate \eqn{\tilde{\sigma}^2} by estimating the long-run
#' variance (LRV) of standardized growth increments using AR(1) pre-whitening
#' and a Bartlett (Newey-West) HAC estimator, followed by recoloring.
#'
#' @param mu Numeric scalar. Estimated growth rate \eqn{\hat{\mu}} used to
#' center increments.
#' @param delta_log_n Numeric vector of length \code{q}. Log-scale increments
#' \eqn{\Delta Y_i}.
#' @param tau Numeric vector of length \code{q}. Sampling intervals
#' \eqn{\tau_i > 0}.
#'
#' @details
#' The estimator targets the LRV of the standardized centered increment series
#' \deqn{U_i = (\Delta Y_i - \hat{\mu}\,\tau_i)/\sqrt{\tau_i},
#' \qquad i=1,\ldots,q,}
#' which is interpreted as an effective environmental variance per unit time.
#'
#' Implementation steps:
#' \enumerate{
#'   \item Center \eqn{u_i} by subtracting \eqn{\bar u} to obtain
#'     \eqn{\tilde u_i}.
#'   \item Estimate the AR(1) pre-whitening coefficient
#'     \eqn{\tilde\rho_{\mathrm{pw}}} by OLS in
#'     \eqn{\tilde u_i = \rho_{\mathrm{pw}}\,\tilde u_{i-1} + \varepsilon_i}
#'     and form pre-whitened residuals \eqn{\tilde\varepsilon_i}.
#'   \item Compute residual autocovariances and the Bartlett HAC LRV
#'     \eqn{\widetilde{\mathcal C}^{(\varepsilon)}_{\mathrm{NW}}(J)}.
#'   \item Choose the truncation lag \eqn{J} via the Andrews (1991) AR(1)
#'     plug-in rule specialized to the Bartlett window.
#'   \item Recolor to the original scale to obtain
#'     \deqn{\tilde{\sigma}^2
#'           = \widetilde{\mathcal C}_{\mathrm{NW}}(J)
#'           = \widetilde{\mathcal C}^{(\varepsilon)}_{\mathrm{NW}}(J)/
#'             (1-\tilde\rho_{\mathrm{pw}})^2.}
#' }
#'
#' This construction is designed to be robust to short-run autocovariance and,
#' under additive observation error on the log scale, to the cancellation
#' property of the LRV when increments are appropriately centered.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{sigma2_tilde}}{Numeric scalar. OEAR diffusion-scale estimate
#'     \eqn{\tilde{\sigma}^2}.}
#'   \item{\code{rho_tilde_pw}}{Numeric scalar. Estimated AR(1) pre-whitening
#'     coefficient \eqn{\tilde\rho_{\mathrm{pw}}}.}
#'   \item{\code{j}}{Integer. Selected Bartlett truncation lag $J$.}
#' }
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#'
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
#' heteroskedasticity and autocorrelation consistent covariance matrix.
#' \emph{Econometrica}, 55(3), 703--708.
#'
#' @keywords internal
#'
#' @name oear_sigma2_hac
NULL

oear_sigma2_hac <- function(mu, delta_log_n, tau) {
  qq <- length(tau)
  u <- (delta_log_n - mu * tau) / sqrt(tau)
  u_bar <- mean(u)
  u_tilde <- u - u_bar

  x <- u_tilde[2L:qq]
  y <- u_tilde[1L:(qq - 1L)]
  den <- sum(y * y)

  if (!is.finite(den) || den <= 0) {
    rho_tilde_pw <- 0
  } else {
    rho_tilde_pw <- sum(x * y) / den
    rho_tilde_pw <- max(min(rho_tilde_pw, 0.999), -0.999)
  }

  e_tilde <- numeric(qq)
  e_tilde[1L] <- u_tilde[1L]
  if (qq >= 2L) {
    e_tilde[2L:qq] <- u_tilde[2L:qq] - rho_tilde_pw * u_tilde[1L:(qq - 1L)]
  }

  j <- floor(((6 * rho_tilde_pw^2) / (1 - rho_tilde_pw^2)^2 * qq)^(1 / 3))
  j <- min(qq - 1L, max(0L, j))

  c_tilde <- function(lag) {
    if (lag == 0L) return(sum(e_tilde * e_tilde) / qq)
    idx <- (lag + 1L):qq
    sum(e_tilde[idx] * e_tilde[idx - lag]) / qq
  }

  c0 <- c_tilde(0L)
  if (j == 0L) {
    c_nw <- c0
  } else {
    bartlett_weighted_sum <- 0
    for (lag in 1L:j) {
      w_lag <- 1 - lag / (j + 1L)
      bartlett_weighted_sum <- bartlett_weighted_sum + w_lag * c_tilde(lag)
    }
    c_nw <- c0 + 2 * bartlett_weighted_sum
  }

  sigma2_tilde <- c_nw / (1 - rho_tilde_pw)^2

  list(
    sigma2_tilde = sigma2_tilde,
    rho_tilde_pw = rho_tilde_pw,
    j = j
  )
}
