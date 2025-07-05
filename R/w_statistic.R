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
#' @title Computes the w Statistic
#'
#' @description
#' Computes the w statistic used for extinction probability estimation
#' in a density-independent population model.
#'
#' @param mu numeric: Estimated population growth rate, \eqn{\hat{\mu}}.
#' @param xd numeric: Distance to extinction threshold on a log scale,
#'   \eqn{x_d = \log(n_q / n_e)}.
#' @param s numeric: Estimated environmental variance, \eqn{\hat{\sigma}^2}.
#' @param th numeric: Time horizon for extinction probability evaluation,
#'   denoted \eqn{t_h}.
#'
#' @details
#' The \eqn{w} statistic is defined as
#' \deqn{
#' w = \frac{\mu t_h + x_d}{\sqrt{\sigma^2 t_h}},
#' }
#' where \eqn{\mu} is the population growth rate, \eqn{\sigma^2} is the
#' environmental variance, \eqn{t_h} is the time horizon, and \eqn{x_d} is
#' the log-distance to the extinction threshold.
#'
#' @return numeric: Value of the w statistic.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @keywords internal
#'
w_statistic <- function(mu, xd, s, th) (mu * th + xd) / sqrt(s * th)
