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
#' @title Computes the w and z Statistics
#'
#' @description
#' Estimators for the \eqn{w} and \eqn{z} statistics used in extinction
#' probability calculations under a density-independent population model.
#'
#' @param mu numeric: Estimated population growth rate, \eqn{\hat{\mu}}.
#' @param xd numeric: Distance to extinction threshold on a log scale,
#'   \eqn{x_d = \log(n_q / n_e)}.
#' @param s numeric: Estimated environmental variance, \eqn{\hat{\sigma}^2}.
#' @param th numeric: Time horizon for extinction probability evaluation,
#'   denoted \eqn{t^{\ast}}.
#'
#' @details
#' The statistics are defined as
#' \deqn{
#' \hat w = \frac{\hat \mu t^{\ast} + x_d}{\sqrt{\hat \sigma^2 t^{\ast}}},
#' \qquad
#' \hat z = \frac{- \hat \mu t^{\ast} + x_d}{\sqrt{\hat \sigma^2 t^{\ast}}}.
#' }
#'
#' @return numeric: Value of the statistic.
#'
#' @name statistics_di
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @keywords internal
#'
NULL

#' @rdname statistics_di
w_statistic <- function(mu, xd, s, th) {
  (mu * th + xd) / sqrt(s * th)
}

#' @rdname statistics_di
z_statistic <- function(mu, xd, s, th) {
  (-mu * th + xd) / sqrt(s * th)
}
