# $Id: $
# Copyright (c) 2025-2026 Hiroshi Hakoyama
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
repr_mode      <- extr:::repr_mode
format_by_mode <- extr:::format_by_mode

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
  stopifnot(is.numeric(res$linear_g))
  stopifnot(res$linear_g >= 0 && res$linear_g <= 1)
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

extr_test_7_oear <- function() {
  time <- 1:30
  population <- round(seq(120, 60, length.out = length(time)))
  dat <- data.frame(time = time, population = population)

  res <- tryCatch(
    ext_di(dat, th = 100, ne = 10, method = "oear"),
    error = function(e) e
  )
  stopifnot(!inherits(res, "error"))
  stopifnot(inherits(res, "ext_di"))
  stopifnot(identical(res$method, "oear"))

  stopifnot(
    is.numeric(res$Variance),
    is.finite(res$Variance),
    res$Variance >= 0
  )
  stopifnot(is.list(res$method_diag))
  stopifnot(
    is.numeric(res$method_diag$rho_tilde_pw),
    length(res$method_diag$rho_tilde_pw) == 1L
  )
  stopifnot(
    is.numeric(res$method_diag$j),
    length(res$method_diag$j) == 1L
  )
  stopifnot(
    res$method_diag$j >= 0,
    res$method_diag$j <= (res$qq - 1)
  )

  stopifnot(is.finite(res$Unbiased.variance), res$Unbiased.variance > 0)
  stopifnot(is.na(res$aic))
  stopifnot(is.finite(res$lower_cl_s), is.finite(res$upper_cl_s))
  stopifnot(res$lower_cl_s > 0, res$upper_cl_s > res$lower_cl_s)
}

extr_test_8_oear_nonpositive_s <- function() {
  dat <- data.frame(
    time = 1:6,
    population = 100 * (1.1)^(0:5)
  )
  err <- tryCatch(
    ext_di(dat, th = 100, ne = 1, method = "oear"),
    error = function(e) e
  )
  stopifnot(inherits(err, "error"))
  stopifnot(grepl("Non-positive or non-finite variance estimate",
                  conditionMessage(err), fixed = TRUE))
}

# --- runner with summary ---
all_extr_tests <- function() {
  passed <- TRUE
  tryCatch({
    extr_test_1()
    extr_test_2()
    extr_test_3()
    extr_test_4()
    extr_test_5()
    extr_test_6()
    extr_test_7_oear()
    extr_test_8_oear_nonpositive_s()
  }, error = function(e) {
    passed <<- FALSE
    message("Test failed: ", conditionMessage(e))
  })

  if (passed) {
    emojis <- c("\U1F389", "\U1F638", "\U1F3C5", "\U1F308")
    cat("All base-R package tests passed", sample(emojis, 1), "\n")
  } else {
    stop("Some package tests failed \u274C")
  }
  invisible(passed)
}

# run unconditionally
all_extr_tests()
