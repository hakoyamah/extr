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
#' @title Confidence Intervals for Extinction Probability (w–z, DI model)
#'
#' @description
#' Computes confidence intervals (CIs) for extinction probability in a
#' density-independent (drifted Wiener) model using the w–z method, and
#' provides a formatter for display-ready CI strings.
#'
#' @param mu Numeric. Estimated growth rate \eqn{\hat{\mu}}.
#' @param xd Numeric. Log-distance to threshold \eqn{x_d=\log(n_q/n_e)}.
#' @param s Numeric. Estimated environmental variance \eqn{\hat{\sigma}^2}.
#' @param th Numeric. Time horizon \eqn{t^{\ast}} for evaluation.
#' @param tq Numeric. Observation span \eqn{t_q} (first to last time).
#' @param qq Integer. Number of intervals \eqn{q} (sample size minus 1).
#' @param alpha Numeric. Significance level \eqn{\alpha}\,=\,1\,-\,CL.
#' @param prob_fun Function. One of \code{ext_prob_di},
#'   \code{log_ext_prob_di}, \code{log_ext_comp_di}.
#' @param digits Integer. Significant digits for formatting (used only by
#'   \code{ci_wz_format_di}).
#'
#' @details
#' The w–z method derives CIs by inverting noncentral-\eqn{t}
#' distributions for the transformed statistics \eqn{w} and \eqn{z},
#' then combining them to form bounds on \eqn{G(w,z)}.
#'
#' Exact confidence intervals for \eqn{w} and \eqn{z} are obtained by
#' numerically solving the noncentral-\eqn{t} quantile equations
#' corresponding to the observed statistics.
#'
#' The CI for \eqn{G(w,z)} is then approximated by
#' \deqn{
#'   \bigl( G(\overline{w},\,\underline{z}),\;
#'          G(\underline{w},\,\overline{z}) \bigr),
#' }
#' where \eqn{\overline{w}}, \eqn{\underline{w}},
#' \eqn{\overline{z}}, and \eqn{\underline{z}} are the upper and lower
#' confidence limits for \eqn{w} and \eqn{z}, respectively.
#' Across the full parameter space, this approach achieves near-nominal
#' coverage. When $z$ is large and positive, $G$ depends only on $w$, so
#' exact CIs are available.
#'
#' The argument \code{prob_fun} selects which probability is evaluated:
#' \itemize{
#'   \item \code{ext_prob_di}, \code{log_ext_prob_di}:
#'     CIs for \eqn{G(w,z)} and \eqn{\log G(w,z)};
#'     returned as (lower, upper).
#'   \item \code{log_ext_comp_di}:
#'     CIs for \eqn{\log Q(w,z)= \log (1-G(w,z))};
#'     returned as (lower, upper).
#' }
#'
#' @return
#' For \code{confidence_interval_wz_di}: numeric vector \code{c(lower, upper)}
#' on the chosen scale (natural log if \code{log_*}).\cr
#' For \code{ci_wz_format_di}: named character vector \code{c(lower, upper)}
#' with values preformatted for display.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @importFrom stats pt uniroot
#'
#' @keywords internal
#'
#' @name wz_ci_di
#' @aliases confidence_interval_wz_di ci_wz_format_di
NULL

#' @rdname wz_ci_di
confidence_interval_wz_di <- function(mu, xd, s, th, tq, qq, alpha,
                                      prob_fun = ext_prob_di) {
  den1 <- sqrt(s * th)
  w_est <- (mu * th + xd) / den1
  z_est <- (- mu * th + xd) / den1
  df <- qq - 1
  const1 <- sqrt((qq - 1) * tq / (qq * th))
  const2 <- sqrt(tq / th)
  find_cl <- function(est, a, lower_tail, width = 10) {
    t_obs <- const1 * est
    f <- function(x) stats::pt(t_obs, df, const2 * x, lower_tail) - a
    d_est <- est / width + 1
    stats::uniroot(f, c(- d_est + est, d_est + est), extendInt = "yes",
                   tol = 1e-14)$root
  }
  lower_cl_w <- find_cl(w_est, alpha / 2, lower_tail = FALSE)
  upper_cl_w <- find_cl(w_est, alpha / 2, lower_tail = TRUE)
  lower_cl_z <- find_cl(z_est, alpha / 2, lower_tail = FALSE)
  upper_cl_z <- find_cl(z_est, alpha / 2, lower_tail = TRUE)
  if (identical(prob_fun, log_ext_comp_di)) {
    lower_cl_p <- prob_fun(lower_cl_w, upper_cl_z)
    upper_cl_p <- prob_fun(upper_cl_w, lower_cl_z)
  } else {
    lower_cl_p <- prob_fun(upper_cl_w, lower_cl_z)
    upper_cl_p <- prob_fun(lower_cl_w, upper_cl_z)
  }
  c(lower_cl_p, upper_cl_p)
}

#' @rdname wz_ci_di
ci_wz_format_di <- function(mu, xd, s, th, tq, qq, alpha,
                            digits = 5L) {
  ci_linear_g <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = ext_prob_di
  )
  ci_log_g <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = log_ext_prob_di
  )
  ci_log_q <- confidence_interval_wz_di(
    mu, xd, s, th, tq, qq, alpha, prob_fun = log_ext_comp_di
  )

  repr_lower <- repr_mode(ci_linear_g[[1]])
  repr_upper <- repr_mode(ci_linear_g[[2]])

  digits <- as.integer(digits)

  lower_str <- format_by_mode(
    mode = repr_lower, kind = "lower",
    linear_g = NA_real_, log_g = NA_real_, log_q = NA_real_,
    ci_linear_g = ci_linear_g, ci_log_g = ci_log_g,
    ci_log_q = ci_log_q, digits = digits
  )
  upper_str <- format_by_mode(
    mode = repr_upper, kind = "upper",
    linear_g = NA_real_, log_g = NA_real_, log_q = NA_real_,
    ci_linear_g = ci_linear_g, ci_log_g = ci_log_g,
    ci_log_q = ci_log_q, digits = digits
  )

  c(lower = lower_str, upper = upper_str)
}
