# $Id: $
# Copyright (c) 2025 Hiroshi Hakoyama
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

library("extr")

extr_test_1 <- function() {
  time <- 1:10
  population <- seq(100, 55, length.out = length(time))
  dat <- data.frame(time = time, population = population)
  res <- ext_di(dat, th = 100, ne = 10)
  stopifnot(inherits(res, "ext_di"))
}

extr_test_2 <- function(tol = 1e-8) {
  time <- 1:20
  population <- seq(100, 50, length.out = length(time))
  dat <- data.frame(time = time, population = population)
  res <- ext_di(dat, th = 50, ne = 10)
  stopifnot(is.numeric(res$Extinction.probability))
  stopifnot(res$Extinction.probability >= 0 && res$Extinction.probability <= 1)
  stopifnot(is.numeric(res$Growth.rate))
  stopifnot(is.numeric(res$Variance))
  stopifnot(is.numeric(res$Unbiased.variance))
  stopifnot(is.numeric(res$aic))
}

extr_test_3 <- function() {
  time <- 1:15
  population <- seq(120, 80, length.out = length(time))
  dat <- data.frame(time = time, population = population)
  res <- ext_di(dat, th = 30, ne = 10)
  out <- tryCatch(print(res, digits = 4), error = function(e) e)
  stopifnot(!inherits(out, "error"))
}

extr_test_4 <- function() {
  # Non-data.frame input
  stopifnot(inherits(tryCatch(ext_di("invalid"), error = function(e) e),
                     "error"))

  # Data frame with wrong columns
  dat_bad <- data.frame(x = 1:10)
  stopifnot(inherits(tryCatch(ext_di(dat_bad), error = function(e) e),
                     "error"))

  # Negative time values
  time <- -1:10
  population <- seq(100, 110, length.out = length(time))
  dat_neg <- data.frame(time = time, population = population)
  stopifnot(inherits(tryCatch(ext_di(dat_neg), error = function(e) e),
                     "error"))

  # ne or th invalid
  time <- 1:10
  population <- seq(100, 50, length.out = length(time))
  dat <- data.frame(time = time, population = population)
  stopifnot(inherits(tryCatch(ext_di(dat, ne = -1), error = function(e) e),
                     "error"))
  stopifnot(inherits(tryCatch(ext_di(dat, th = 0), error = function(e) e),
                     "error"))
}

# Ensure the function can handle uneven time intervals
extr_test_5 <- function() {
  # Create data with uneven time intervals
  time <- c(1, 2, 4, 7, 11, 16, 22, 29)
  population <- seq(100, 50, length.out = length(time))
  dat <- data.frame(time = time, population = population)

  # Run estimation and check no errors occur
  res <- tryCatch(ext_di(dat, th = 30, ne = 10), error = function(e) e)
  stopifnot(!inherits(res, "error"))

  # Check the result inherits the correct class
  stopifnot(inherits(res, "ext_di"))
}

# Verify that providing a non-"years" unit string works correctly
extr_test_6 <- function() {
  time <- 1:10
  population <- seq(100, 50, length.out = length(time))
  dat <- data.frame(time = time, population = population)

  # Use a custom unit string (e.g., "months")
  res <- tryCatch(ext_di(dat, th = 30, ne = 10, unit = "months"),
                  error = function(e) e)
  stopifnot(!inherits(res, "error"))
  stopifnot(inherits(res, "ext_di"))

  # Check that the output contains the custom unit string
  printed <- capture.output(print(res))
  found_unit <- any(grepl("months", printed))
  stopifnot(found_unit)
}


# Run all tests
extr_test_1()
extr_test_2()
extr_test_3()
extr_test_4()
extr_test_5()
extr_test_6()
