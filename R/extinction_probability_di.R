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
#' @title Computes Extinction Probability for a Density-Independent Model
#'
#' @description
#' Compute \eqn{G(w,z)} and its complement \eqn{Q(w,z)=1-G(w,z)} for a
#' density-independent (drifted Wiener) model, on both linear and log scales.
#'
#' @details
#' For any \eqn{t>0} with \eqn{w+z>0},
#' \deqn{
#'   \Pr[T \leq t] = G(w,z)=\Phi(-w)+
#'   \exp\!\left(\tfrac{z^2-w^2}{2}\right)\Phi(-z), \qquad
#'   Q(w,z)=1-G(w,z).
#' }
#' Here \eqn{\Phi} and \eqn{\phi} denote the standard normal CDF and PDF.
#'
#' \strong{Stability strategy.}
#' (i) For large \eqn{z}, rewrite the product
#' \eqn{\exp((z^2-w^2)/2)\,\Phi(-z)} via the Mills ratio and replace it by an
#' 8-term asymptotic series when \eqn{z \ge 19}.
#' (ii) On the log scale, use log-sum-exp and a stable log-difference
#' (\code{log1mexp}) built from \code{log1p}/\code{expm1} to retain tail info.
#'
#' \strong{Domain.} Scalar inputs are assumed and require \eqn{w+z>0}.
#'
#' \strong{Functions.}
#' \itemize{
#'   \item \code{ext_prob_di(w,z)}: returns \eqn{G(w,z)} (linear scale).
#'   \item \code{log_ext_prob_di(w,z)}: returns \eqn{\log G(w,z)}.
#'   \item \code{log_ext_comp_di(w,z)}: returns \eqn{\log Q(w,z)}.
#'   \item \code{ext_prob_format_di(w,z,digits)}: formats a point estimate
#'     using \code{repr_mode()} and \code{format_by_mode()}.
#' }
#'
#' @param w Numeric; transformed parameter
#'   \eqn{w=(\mu t+x_d)/(\sigma\sqrt{t})}.
#' @param z Numeric; transformed parameter
#'   \eqn{z=(-\mu t+x_d)/(\sigma\sqrt{t})}.
#' @param digits Integer; significant digits for formatting (only for
#'   \code{ext_prob_format_di}).
#'
#' @return
#' For \code{ext_prob_di}: numeric scalar \eqn{G(w,z)}.\cr
#' For \code{log_ext_prob_di}: numeric scalar \eqn{\log G(w,z)}.\cr
#' For \code{log_ext_comp_di}: numeric scalar \eqn{\log Q(w,z)}.\cr
#' For \code{ext_prob_format_di}: character string formatted for display.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @importFrom stats pnorm dnorm
#'
#' @keywords internal
#'
#' @seealso \code{\link{statistics_di}}, \code{\link{repr_mode}},
#'   \code{\link{format_by_mode}}
#'
#' @name extinction_probability_di
NULL

# --- G(w, z) ---
#' @rdname extinction_probability_di
ext_prob_di <- function(w, z) {
  if (z < 19) {
    pnorm(-w) + exp(0.5 * (z + w) * (z - w)) * pnorm(-z)
  } else {
    pnorm(-w) + dnorm(w) * .mills_ratio_s7(z)
  }
}

# --- log G(w, z): log-sum-exp ---
#' @rdname extinction_probability_di
log_ext_prob_di <- function(w, z) {
  a <- pnorm(-w, log.p = TRUE)
  if (z < 19) {
    b <- 0.5 * (z + w) * (z - w) + pnorm(-z, log.p = TRUE)
  } else {
    b <- dnorm(w, log = TRUE) + log(.mills_ratio_s7(z))
  }
  m <- max(a, b)
  m + log1p(exp(-abs(a - b)))
}

# --- log Q(w, z): log-difference-exp ---
#' @rdname extinction_probability_di
log_ext_comp_di <- function(w, z) {
  c <- pnorm(w, log.p = TRUE)
  if (z < 19) {
    b <- 0.5 * (z + w) * (z - w) + pnorm(-z, log.p = TRUE)
  } else {
    b <- dnorm(w, log = TRUE) + log(.mills_ratio_s7(z))
  }
  c + log1p(-exp(b - c))
}

# --- 8-term Mills ratio series ---
.mills_ratio_s7 <- function(z) {
  1 / z - 1 / z^3 + 3 / z^5 - 15 / z^7 +
    105 / z^9  - 945 / z^11 + 10395 / z^13 - 135135 / z^15
}

# --- formatted G(w, z) ---
#' @rdname extinction_probability_di
ext_prob_format_di <- function(w, z, digits = 5L) {
  if (w + z <= 0) {
    stop("Invalid input: require w+z > 0")
  }
  linear_g <- ext_prob_di(w, z)
  log_g <- log_ext_prob_di(w, z)
  log_q <- log_ext_comp_di(w, z)
  mode <- repr_mode(linear_g)
  format_by_mode(
    mode = mode, kind = "point",
    linear_g = linear_g, log_g = log_g, log_q = log_q,
    ci_linear_g = c(NA_real_, NA_real_),
    ci_log_g = c(NA_real_, NA_real_),
    ci_log_q = c(NA_real_, NA_real_),
    digits = as.integer(digits)
  )
}
