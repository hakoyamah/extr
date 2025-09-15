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
repr_mode      <- extr:::repr_mode
format_by_mode <- extr:::format_by_mode

run_numeric_tests <- function() {
  tol <- 1e-12
  passed <- TRUE
  check <- function(cond, msg) {
    if (!cond) {
      passed <<- FALSE
      message("\u274C Test failed: ", msg)
    }
  }

  # ---- G + Q = 1 (Q via log_ext_comp_di) ----
  pts <- list(
    c(w = 0.5, z = 1.2),
    c(w = -3.0, z = 6.0),
    c(w = 3.0, z = 6.0),
    c(w = 0.0, z = 25.0)
  )
  for (p in pts) {
    w <- p["w"]
    z <- p["z"]
    g_val <- extr:::ext_prob_di(w, z)
    log_qval <- extr:::log_ext_comp_di(w, z)
    q_val <- exp(log_qval)
    check(is.finite(g_val) && is.finite(q_val), "G or Q not finite")
    check(
      abs((g_val + q_val) - 1) < 1e-14,
      sprintf("G+Q != 1 at (w=%.2f, z=%.2f)", w, z)
    )
  }

  # ---- log(G)/log(Q) vs linear (skip extreme tails) ----
  for (p in pts) {
    w <- p["w"]
    z <- p["z"]
    g_val <- extr:::ext_prob_di(w, z)
    q_lin <- 1 - g_val
    if (g_val > 1e-300 && g_val < 1 - 1e-15) {
      log_g <- extr:::log_ext_prob_di(w, z)
      check(
        abs(log_g - log(g_val)) < tol,
        sprintf("log(G) mismatch at (%.2f, %.2f)", w, z)
      )
    }
    if (q_lin > 1e-300 && q_lin < 1 - 1e-15) {
      log_q <- extr:::log_ext_comp_di(w, z)
      check(
        abs(log_q - log(q_lin)) < tol,
        sprintf("log(Q) mismatch at (%.2f, %.2f)", w, z)
      )
    }
  }

  # ---- Mills switch continuity at z = 19 ----
  w <- 0.7
  z1 <- 19 - 1e-8
  z2 <- 19 + 1e-8
  check(
    abs(extr:::ext_prob_di(w, z1) - extr:::ext_prob_di(w, z2)) < 1e-10,
    "Continuity of G at z=19"
  )
  check(
    abs(
      extr:::log_ext_comp_di(w, z1) - extr:::log_ext_comp_di(w, z2)
    ) < 1e-10,
    "Continuity of log(Q) at z=19"
  )
  check(
    abs(
      extr:::log_ext_prob_di(w, z1) - extr:::log_ext_prob_di(w, z2)
    ) < 1e-10,
    "Continuity of log(G) at z=19"
  )

  if (passed) {
    emojis <- c("\U1F389", "\U1F638", "\U1F3C5", "\U1F308")
    cat("All base-R numeric tests passed", sample(emojis, 1), "\n")
  } else {
    stop("Some numeric tests failed \u274C")
  }
}

run_numeric_tests()
