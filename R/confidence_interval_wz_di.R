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
#' @title Computes Confidence Intervals of Extinction Risk for
#' a Density-Independent Population Model
#'
#' @description
#' Calculates the confidence interval for the extinction probability of a
#' density-independent population model using the \eqn{w}–\eqn{z} method
#' (Hakoyama, 2025).
#'
#' @param mu numeric: Estimated population growth rate, \eqn{\hat{\mu}}.
#' @param xd numeric: Distance to extinction threshold on a log scale,
#'   \eqn{x_d = \log(n_q / n_e)}.
#' @param s numeric: Estimated environmental variance, \eqn{\hat{\sigma}^2}.
#' @param th numeric: Time horizon for extinction probability evaluation,
#' \eqn{t_h}.
#' @param tq numeric: Total observation period (from first to last data point),
#' \eqn{t_q}.
#' @param qq integer: Number of time intervals (sample size minus 1), \eqn{q}.
#' @param alpha numeric: Significance level, i.e., \eqn{\alpha} =
#' 1 - (confidence level).
#'
#' @details
#' To improve extinction probability estimation, two transformed parameters,
#' \eqn{w} and \eqn{z}, are introduced as functions of the growth rate \eqn{\mu}
#' and environmental variance \eqn{\sigma^2}:
#' \deqn{
#' w = \frac{\mu t + x_d}{\sigma \sqrt{t}}, \quad
#' z = \frac{-\mu t + x_d}{\sigma \sqrt{t}},
#' }
#' where \eqn{t} is the time horizon and \eqn{x_d = \log(n_q / n_e)} is the
#' log-distance to the extinction threshold.
#'
#' The extinction probability within \eqn{t} is expressed as:
#' \deqn{
#' \Pr[T \leq t] = G(w, z) = \Phi(-w) +
#' \exp\left(\frac{z^2 - w^2}{2}\right) \Phi(-z),
#' }
#' where \eqn{\Phi(\cdot)} is the standard normal cumulative distribution
#' function.
#'
#' The maximum likelihood estimators \eqn{\widehat{w}} and \eqn{\widehat{z}}
#' follow non-central \eqn{t}-distributions with \eqn{q-1} degrees of freedom:
#' \deqn{
#' \widehat{w} \sqrt{\frac{q-1}{q}} \sqrt{\frac{t_q}{t}} \sim t(\delta_w, q-1),
#' \quad
#' \widehat{z} \sqrt{\frac{q-1}{q}} \sqrt{\frac{t_q}{t}} \sim t(\delta_z, q-1),
#' }
#' with drift parameters
#' \deqn{
#' \delta_w = w \sqrt{\frac{t_q}{t}}, \quad
#' \delta_z = z \sqrt{\frac{t_q}{t}}.
#' }
#'
#' Asymptotically exact \eqn{(1-\alpha)} confidence intervals for \eqn{w} and
#' \eqn{z} are obtained by numerically solving the non-central \eqn{t} quantile
#' equations corresponding to the observed statistics.
#'
#' The confidence interval for the extinction probability \eqn{G(w, z)} is then
#' approximated by combining these intervals as
#' \deqn{
#' \left( G(\overline{\Psi_w}, \underline{\Psi_z}), \quad
#' G(\underline{\Psi_w}, \overline{\Psi_z}) \right),
#' }
#' where \eqn{\overline{\Psi_w}}, \eqn{\underline{\Psi_w}},
#' \eqn{\overline{\Psi_z}}, and \eqn{\underline{\Psi_z}} are the upper and lower
#' confidence limits for \eqn{w} and \eqn{z}, respectively.
#'
#' This \emph{w-z method} yields exact confidence intervals when \eqn{z \gg 0}
#' and robust approximations when \eqn{z \ll 0}.
#'
#' @return A numeric vector of length two containing the lower and upper
#'   confidence limits of the extinction probability.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @importFrom stats pt uniroot
#'
#' @keywords internal
#'
confidence_interval_wz_di <- function(mu, xd, s, th, tq, qq, alpha) {
  den1 <- sqrt(s * th)
  w_est <- (mu * th + xd) / den1
  z_est <- (- mu * th + xd) / den1
  df <- qq - 1
  const1 <- sqrt((qq - 1) * tq / (qq * th))
  const2 <- sqrt(tq / th)
  find_cl <- function(tq, qq, th, est, a, width = 10) {
    t_obs <- const1 * est
    f <- function(x) pt(t_obs, df, const2 * x) - a
    d_est <- est / width + 1
    uniroot(f, c(- d_est + est, d_est + est), extendInt = "yes")$root
  }
  lower_cl_w <- find_cl(tq, qq, th, w_est, 1 - alpha / 2)
  upper_cl_w <- find_cl(tq, qq, th, w_est, alpha / 2)
  lower_cl_z <- find_cl(tq, qq, th, z_est, 1 - alpha / 2)
  upper_cl_z <- find_cl(tq, qq, th, z_est, alpha / 2)
  lower_cl_p <- extinction_probability_wz_di(upper_cl_w, lower_cl_z)
  upper_cl_p <- extinction_probability_wz_di(lower_cl_w, upper_cl_z)
  c(lower_cl_p, upper_cl_p)
}
