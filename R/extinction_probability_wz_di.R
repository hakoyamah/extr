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
#' Computes the extinction probability of a density-independent stochastic
#' population model based on the given w and z statistics.
#'
#' @param w numeric: w statistic for extinction probability calculation.
#' @param z numeric: z statistic for extinction probability calculation.
#'
#' @details
#' Let \eqn{T} be a random variable representing the time at which the
#' population becomes extinct. The extinction probability within a given
#' time horizon \eqn{t_h} is defined as \eqn{\Pr[T \leq t_h]}.
#'
#' Hakoyama (2025) expressed this extinction probability as a function of two
#' transformed statistics, \eqn{w} and \eqn{z}, yielding the following formula:
#' \deqn{
#' \Pr[T \leq t_h] = G(w, z) = \Phi(-w) + \exp{\left(\frac{z^2 - w^2}{2}\right)}
#' \Phi(-z),
#' }
#' where \eqn{\Phi(\cdot)} is the cumulative distribution function of the
#' standard normal distribution.
#'
#' For large values of \eqn{z} (typically \eqn{z \geq 35}), this exact formula
#' may suffer from numerical inaccuracies. Dennis et al. (1991) derived an
#' asymptotic approximation for extinction probabilities expressed directly in
#' terms of \eqn{\mu}, \eqn{\sigma^2}, \eqn{n_q}, \eqn{n_e}, and \eqn{t_h}.
#'
#' This package reformulates their approximation in terms of \eqn{w} and
#' \eqn{z}, resulting in the following asymptotic expansion:
#' \deqn{
#' G(w, z) \approx \Phi(-w) + \exp{\left(-\frac{w^2}{2}\right)}
#' \frac{\sqrt{2}}{2\sqrt{\pi}}
#' \left(\frac{1}{z} - \frac{1}{z^3} + \frac{3}{z^5} - \frac{15}{z^7}
#' + \frac{105}{z^9} - \frac{945}{z^{11}} + \frac{10395}{z^{13}}\right).
#' }
#'
#' This approximation improves numerical stability while maintaining high
#' accuracy for large \eqn{z}.
#'
#' @return numeric: Estimated extinction probability.
#'
#' @author Hiroshi Hakoyama, \email{hiroshi.hakoyama@gmail.com}
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
#'
#' @seealso \code{\link{w_statistic}}, \code{\link{z_statistic}}
#'
extinction_probability_wz_di <- function(w, z) {
  if (z < 35) {
    pnorm(-w) + exp((z^2 - w^2) / 2) * pnorm(-z)
  } else {
    pnorm(-w) + exp(- w^2 / 2) * (sqrt(2) / (2 * sqrt(pi))) *
      (1 / z - 1 / z^3 + 3 / z^5 - 15 / z^7
       + 105 / z^9 - 945 / z^11 + 10395 / z^13)
  }
}
