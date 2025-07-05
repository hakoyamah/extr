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
#' @title Print Extinction Risk Estimations
#'
#' @param x an object of class \code{"ext_di"} produced by
#' \code{\link{ext_di}}.
#' @param digits integer: Number of significant digits to display. Default is 5.
#'
#' @return Invisibly returns \code{obj} after printing formatted extinction risk
#'   estimation results, including estimates, confidence intervals, data
#'   summary, and input parameters.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @keywords internal
#'
#' @export
#'
print.ext_di <- function(x, digits = 5, ...) {
  est_table <- data.frame(
    Estimate = c(
      formatC(x$Extinction.probability, digits = digits, format = "g"),
      formatC(x$Growth.rate, digits = digits, format = "g"),
      formatC(x$Variance, digits = digits, format = "g"),
      formatC(x$Unbiased.variance, digits = digits, format = "g"),
      formatC(x$aic, digits = digits, format = "g")
    ),
    CI = c(
      paste0("(", formatC(x$lower_cl_p, digits = digits, format = "g"), ", ",
             formatC(x$upper_cl_p, digits = digits, format = "g"), ")"),
      paste0("(", formatC(x$lower_cl_mu, digits = digits, format = "g"), ", ",
             formatC(x$upper_cl_mu, digits = digits, format = "g"), ")"),
      paste0("(", formatC(x$lower_cl_s, digits = digits, format = "g"), ", ",
             formatC(x$upper_cl_s, digits = digits, format = "g"), ")"),
      "-",
      "-"
    )
  )
  rownames(est_table) <- c(
    paste0("Probability of decline to ", x$ne, " within ", x$th, " ",
           x$unit, " (MLE):"),
    "Growth rate (MLE):",
    "Variance (MLE):",
    "Unbiased variance:",
    "AIC for the distribution of N:"
  )

  data_table <- data.frame(
    Value = c(
      formatC(x$nq, digits = digits, format = "g"),
      formatC(x$xd, digits = 4, format = "f"),
      formatC(x$sample.size, digits = digits, format = "g")
    )
  )
  rownames(data_table) <- c(
    "Current population size, nq:",
    "xd = ln(nq / ne):",
    "Sample size, q + 1:"
  )

  input_params <- data.frame(
    Parameter = c(
      x$unit,
      formatC(x$ne, digits = 0, format = "f"),
      formatC(x$th, digits = 0, format = "f"),
      formatC(x$alpha, digits = digits, format = "g")
    )
  )
  rownames(input_params) <- c(
    "Time unit:",
    "Extinction threshold of population size, ne:",
    paste0("Time window for extinction risk evaluation (", x$unit, "), th:"),
    "Significance level, alpha:"
  )

  cat("--- Estimates ---\n")
  print(est_table)
  cat("\n--- Data Summary ---\n")
  print(data_table)
  cat("\n--- Input Parameters ---\n")
  print(input_params)
}
