# $Id: $
# Copyright (c) 2013-2026 Hiroshi Hakoyama
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
#' @title Print Extinction Risk Estimates
#'
#' @description
#' S3 method that prints a formatted summary for an \code{ext_di} object.
#' The output includes the plug-in extinction probability and its CI,
#' method-specific parameter estimates, a brief data summary,
#' and input settings.
#'
#' @param x An object of class \code{"ext_di"} as returned by
#'   \code{\link{ext_di}}.
#' @param digits Integer. Significant digits to display; defaults to
#'   \code{x$digits} when \code{NULL}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x} after printing the formatted results.
#'
#' @keywords internal
#'
#' @method print ext_di
#'
#' @export
#'

print.ext_di <- function(x, digits = NULL, ...) {
  digits <- as.integer(if (is.null(digits)) x$digits else digits)

  fmt_num <- function(val, k = digits) formatC(val, digits = k, format = "g")

  # ---------- strings ----------
  w <- w_statistic(x$Growth.rate, x$xd, x$Variance, x$th)
  z <- z_statistic(x$Growth.rate, x$xd, x$Variance, x$th)
  p_str <- ext_prob_format_di(w, z, digits = digits)

  # CI only (no point here)
  ci_fmt <- ci_wz_format_di(
    mu = x$Growth.rate, xd = x$xd, s = x$Variance,
    th = x$th, tq = x$tq, qq = x$qq, alpha = x$alpha,
    digits = digits
  )
  pr_ci <- paste0("(", ci_fmt["lower"], ", ", ci_fmt["upper"], ")")

  method <- x$method
  is_oear <- (method == "oear")

  # --- table ---
  if (is_oear) {
    est_table <- data.frame(
      Estimate = c(
        p_str,
        fmt_num(x$Growth.rate, digits),
        fmt_num(x$Variance, digits),
        fmt_num(x$method_diag$rho_tilde_pw, digits),
        fmt_num(x$method_diag$j, digits)
      ),
      CI = c(
        pr_ci,
        paste0("(", fmt_num(x$lower_cl_mu, digits), ", ",
               fmt_num(x$upper_cl_mu, digits), ")"),
        paste0("(", fmt_num(x$lower_cl_s, digits), ", ",
               fmt_num(x$upper_cl_s, digits), ")"),
        "-",
        "-"
      ),
      row.names = c(
        paste0("Probability of decline to ", x$ne,
               " within ", x$th, " ", x$unit, " (OEAR plug-in):"),
        "Growth rate (MLE):",
        "Process variance (OEAR):",
        "AR(1) pre-whitening rho:",
        "Bartlett lag truncation (j):"
      ),
      check.names = FALSE
    )
  } else {
    est_table <- data.frame(
      Estimate = c(
        p_str,
        fmt_num(x$Growth.rate, digits),
        fmt_num(x$Variance, digits),
        fmt_num(x$Unbiased.variance, digits),
        fmt_num(x$aic, digits)
      ),
      CI = c(
        pr_ci,
        paste0("(", fmt_num(x$lower_cl_mu, digits), ", ",
               fmt_num(x$upper_cl_mu, digits), ")"),
        paste0("(", fmt_num(x$lower_cl_s, digits), ", ",
               fmt_num(x$upper_cl_s, digits), ")"),
        "-",
        "-"
      ),
      row.names = c(
        paste0("Probability of decline to ", x$ne,
               " within ", x$th, " ", x$unit, " (MLE):"),
        "Growth rate (MLE):",
        "Environmental variance (MLE):",
        "Unbiased variance:",
        "AIC for the distribution of N:"
      ),
      check.names = FALSE
    )
  }

  data_table <- data.frame(
    Value = c(
      fmt_num(x$nq, digits),
      formatC(x$xd, digits = 4, format = "f"),
      fmt_num(x$sample.size, digits)
    ),
    row.names = c(
      "Current population size, nq:",
      "xd = ln(nq / ne):",
      "Sample size, q + 1:"
    ),
    check.names = FALSE
  )

  input_params <- data.frame(
    Parameter = c(
      x$unit,
      formatC(x$ne, digits = 0, format = "f"),
      formatC(x$th, digits = 1, format = "f"),
      fmt_num(x$alpha, digits)
    ),
    row.names = c(
      "Time unit:",
      "Extinction threshold of population size, ne:",
      paste0("Time window for extinction risk evaluation (",
             x$unit, "), th:"),
      "Significance level, alpha:"
    ),
    check.names = FALSE
  )

  cat("--- Estimates ---\n")
  print(est_table)
  cat("\n--- Data Summary ---\n")
  print(data_table)
  cat("\n--- Input Parameters ---\n")
  print(input_params)
  invisible(x)
}
