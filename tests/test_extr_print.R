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

library(extr)
repr_mode <- extr:::repr_mode
format_by_mode <- extr:::format_by_mode

message("=== extr print/format tests start ===")

# --- helper assertions ---
aeq <- function(x, y, tol = 1e-12, msg = "not equal") {
  if (is.na(x) && is.na(y)) return(invisible(TRUE))
  if (!isTRUE(all.equal(x, y, tolerance = tol))) {
    stop(sprintf("AEQ FAIL: %s\nx=%s\ny=%s", msg, deparse(x), deparse(y)))
  }
  invisible(TRUE)
}

atrue <- function(cond, msg = "assertion failed") {
  if (!isTRUE(cond)) stop(sprintf("ASSERT FAIL: %s", msg))
  invisible(TRUE)
}

acontains <- function(lines, pattern, msg = "pattern not found") {
  if (!any(grepl(pattern, lines))) {
    stop(sprintf("NOT FOUND: %s\npattern: %s", msg, pattern))
  }
  invisible(TRUE)
}

# --- 1) repr_mode ---
test_repr_mode <- function() {
  aeq(repr_mode(0.2), "linear_G", msg = "linear case")
  aeq(repr_mode(1e-310), "log_G", msg = "tiny probability")
  aeq(repr_mode(1 - 1e-20), "log_Q", msg = "near one")
  aeq(repr_mode(NaN), "undefined", msg = "NaN handling")
  message("ok: repr_mode")
}

# --- 2) format_by_mode ---
test_format_by_mode <- function() {
  out_p <- format_by_mode(
    mode = "linear_G", kind = "point",
    linear_g = 0.123456, log_g = NA_real_, log_q = NA_real_,
    ci_linear_g = c(0.01, 0.9), ci_log_g = c(NA_real_, NA_real_),
    ci_log_q = c(NA_real_, NA_real_), digits = 5
  )
  acontains(out_p, "^0\\.12346", "linear point")

  out_lo <- format_by_mode(
    "linear_G", "lower", 0.2, NA, NA, c(0.0034017, 0.99), c(NA, NA),
    c(NA, NA), digits = 6
  )
  acontains(out_lo, "^0\\.0034017(0)?$", "linear lower")

  out_lg <- format_by_mode(
    "log_G", "point", NA, log(1e-20), NA,
    c(NA, NA), c(log(1e-25), log(1e-15)), c(NA, NA), digits = 5
  )
  atrue(grepl("^10\\^\\([^()]+\\)$", out_lg), "log_G shape")
  exp_str <- sub("^10\\^\\(([^()]+)\\)$", "\\1", out_lg)
  exp_num <- suppressWarnings(as.numeric(exp_str))
  atrue(is.finite(exp_num), "log_G exponent numeric")
  aeq(exp_num, -20, tol = 1e-8, msg = "log_G exponent ~= -20")

  g_near1 <- 1 - 1e-9
  out_lq <- format_by_mode(
    "log_Q", "point", NA, NA, log(1 - g_near1),
    c(NA, NA), c(NA, NA),
    c(log(1 - (1 - 1e-10)), log(1 - (1 - 1e-8))),
    digits = 3
  )
  atrue(grepl("^1 - 10\\^\\(", out_lq), "log_Q point")

  out_nf <- format_by_mode(
    "log_G", "upper", NA, NA, NA,
    c(NA, NA), c(NA, NA), c(NA, NA), digits = 5
  )
  aeq(out_nf, "NA", "non-finite returns NA")

  message("ok: format_by_mode")
}

# --- 3) ext_di core ---
test_ext_di_core <- function() {
  set.seed(1)
  time <- 1:30
  pop <- round(exp(seq(log(200), log(80), length.out = length(time))))
  dat <- data.frame(time = time, population = pop)

  res <- ext_di(dat, th = 50, ne = 10, digits = 6)
  atrue(inherits(res, "ext_di"), "class ext_di")
  needed <- c(
    "linear_g", "log_g", "log_q", "ci_linear_g", "ci_log_g", "ci_log_q",
    "repr_point", "repr_lower", "repr_upper", "digits"
  )
  atrue(all(needed %in% names(res)), "fields present")

  atrue(is.numeric(res$Growth.rate), "mu numeric")
  atrue(is.numeric(res$Variance), "s numeric")
  atrue(is.numeric(res$Unbiased.variance), "us numeric")
  atrue(is.numeric(res$aic), "aic numeric")

  atrue(length(res$ci_linear_g) == 2, "ci_linear length 2")
  atrue(length(res$ci_log_g) == 2, "ci_log_g length 2")
  atrue(length(res$ci_log_q) == 2, "ci_log_q length 2")

  message("ok: ext_di core")
}

# --- 4) print.ext_di behavior ---
test_print_ext_di <- function() {
  time <- 1:20
  pop <- round(seq(100, 50, length.out = length(time)))
  dat <- data.frame(time = time, population = pop)

  res_small <- ext_di(dat, th = 4, ne = 1, digits = 5)
  lines_small <- capture.output(print(res_small))
  atrue(any(grepl("10\\^\\(", lines_small)), "log form appears")

  res_large <- ext_di(dat, th = 250, ne = 1, digits = 10)
  lines_large <- capture.output(print(res_large))
  atrue(any(grepl("0\\.99999|1 - 10\\^\\(", lines_large)),
        "high precision near 1")

  lines_7 <- capture.output(print(res_large, digits = 7))
  atrue(any(grepl("0\\.999999|1 - 10\\^\\(", lines_7)),
        "print(digits=7) respected")

  message("ok: print.ext_di")
}

# --- 5) error handling ---
test_errors <- function() {
  bad <- tryCatch(ext_di("bad"), error = identity)
  atrue(inherits(bad, "error"), "non-df error")

  bad2 <- tryCatch(
    ext_di(data.frame(a = 1:3, b = c(1, 0, 2))),
    error = identity
  )
  atrue(inherits(bad2, "error"), "n<1 error")

  dat <- data.frame(
    time = c(0, 1, 1, 2),
    population = c(10, 11, 12, 13)
  )
  bad3 <- tryCatch(ext_di(dat), error = identity)
  atrue(inherits(bad3, "error"), "non-increasing time")

  message("ok: error handling")
}

# --- runner with summary ---
run_all_extr_base_tests <- function() {
  passed <- TRUE
  tryCatch({
    test_repr_mode()
    test_format_by_mode()
    test_ext_di_core()
    test_print_ext_di()
    test_errors()
  }, error = function(e) {
    passed <<- FALSE
    message("Test failure: ", conditionMessage(e))
  })

  if (passed) {
    emojis <- c("\U1F389", "\U1F638", "\U1F3C5", "\U1F308")
    cat("All base-R numeric tests passed", sample(emojis, 1), "\n")
  } else {
    stop("Some numeric tests failed \u274C")
  }
  invisible(passed)
}

run_all_extr_base_tests()
