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
#' @title Display Modes for Extinction Probabilities
#'
#' @description
#' A pair of internal functions for formatting extinction probabilities:
#' \code{repr_mode()} decides how to represent a probability \eqn{G} based on
#' its magnitude, and \code{format_by_mode()} formats numbers accordingly for
#' presentation. These functions affect display only and do not alter
#' underlying calculations.
#'
#' @details
#' Probabilities very close to 0 or 1 are displayed on a log scale to improve
#' readability. The thresholds are fixed at \code{1e-300} ("near zero") and
#' \code{1 - 1e-15} ("near one"). Non-finite inputs are reported as
#' \code{"undefined"}.
#'
#' @param g Numeric scalar, a probability in the unit interval. May be
#'   non-finite.
#'
#' @return
#' \code{repr_mode()} returns a character scalar, one of:
#' \itemize{
#'   \item \code{"linear_G"} – display \eqn{G} on the linear scale;
#' \item \code{"log_G"} – display \eqn{G} on log10 scale for tiny
#'   probabilities;
#' \item \code{"log_Q"} – display \eqn{Q=1-G} on log10 scale when
#'   \eqn{G \approx 1};
#'   \item \code{"undefined"} – input was not finite.
#' }
#'
#' \code{format_by_mode()} returns a character scalar formatted according to
#' the selected mode and the kind of value (point estimate, lower, or upper
#' CI bound).
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @keywords internal
#'
repr_mode <- function(g) {
  eps_low <- 1e-300
  eps_high <- 1 - 1e-15
  if (!is.finite(g)) return("undefined")
  if (g <= eps_low)  return("log_G")
  if (g >= eps_high) return("log_Q")
  "linear_G"
}

#' @title Format a Probability using a Display Mode
#'
#' @description
#' Format a probability (point estimate or CI bound) according to a chosen
#' representation mode. Allowed modes are \code{"linear_G"}, \code{"log_G"},
#' \code{"log_Q"}, and \code{"undefined"}. The function performs no inference;
#' it only formats the numeric values already computed (linear and log).
#'
#' Notes on arguments:
#' \itemize{
#'   \item \code{log_g}, \code{log_q}, \code{ci_log_g}, \code{ci_log_q} are
#'     natural logarithms (from \code{log()} in R).
#'   \item For \code{"log_Q"}: the upper bound for \eqn{G} uses
#'     \code{ci_log_q[1]} and the lower bound uses \code{ci_log_q[2]}.
#'     This ordering reflects that \eqn{Q = 1-G} decreases as \eqn{G} increases.
#' }
#'
#' @param mode Character scalar. One of \code{"linear_G"}, \code{"log_G"},
#'   \code{"log_Q"}, \code{"undefined"}.
#' @param kind Character scalar. One of \code{"point"}, \code{"lower"},
#'   \code{"upper"}.
#' @param linear_g Numeric scalar. Point estimate on the linear scale.
#' @param log_g Numeric scalar. \code{log(G)} (natural log). May be non-finite.
#' @param log_q Numeric scalar. \code{log(Q)} with \eqn{Q=1-G} (natural log).
#'   May be non-finite.
#' @param ci_linear_g Numeric length-2: \code{c(lower, upper)} for \eqn{G}.
#' @param ci_log_g Numeric length-2: \code{c(lower_logG, upper_logG)}.
#' @param ci_log_q Numeric length-2: \code{c(upper_logQ, lower_logQ)}.
#' @param digits Integer. Significant digits for numeric formatting.
#'
#' @return A character scalar formatted for display. Returns \code{"NA"} if the
#'   required value is non-finite or missing.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @keywords internal
#'
format_by_mode <- function(mode,
                           kind,
                           linear_g,
                           log_g,
                           log_q,
                           ci_linear_g,
                           ci_log_g,
                           ci_log_q,
                           digits) {
  fmt_num <- function(val, k = digits) formatC(val, digits = k, format = "g")
  to_log10 <- function(v) if (is.finite(v)) v / log(10) else NA_real_
  fmt_log_g <- function(log10_g) {
    if (is.finite(log10_g)) paste0("10^(", fmt_num(log10_g, digits), ")")
    else "NA"
  }
  fmt_log_q <- function(log10_q) {
    if (is.finite(log10_q)) paste0("1 - 10^(", fmt_num(log10_q, digits), ")")
    else "NA"
  }

  if (identical(mode, "linear_G")) {
    val <- switch(
      kind,
      point = linear_g,
      lower = ci_linear_g[[1]],
      upper = ci_linear_g[[2]]
    )
    return(if (is.finite(val)) fmt_num(val, digits) else "NA")
  }

  if (identical(mode, "log_G")) {
    logv <- switch(
      kind,
      point = log_g,
      lower = if (length(ci_log_g) >= 1) ci_log_g[[1]] else NA_real_,
      upper = if (length(ci_log_g) >= 2) ci_log_g[[2]] else NA_real_
    )
    return(fmt_log_g(to_log10(logv)))
  }

  if (identical(mode, "log_Q")) {
    # lower bound for G -> ci_log_q[2], upper bound for G -> ci_log_q[1]
    logv <- switch(
      kind,
      point = log_q,
      lower = if (length(ci_log_q) >= 2) ci_log_q[[2]] else NA_real_,
      upper = if (length(ci_log_q) >= 1) ci_log_q[[1]] else NA_real_
    )
    return(fmt_log_q(to_log10(logv)))
  }

  "NA"
}
