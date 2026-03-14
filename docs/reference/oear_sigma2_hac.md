# OEAR Diffusion-Scale Estimator

Compute the OEAR (observation-error-and-autocovariance-robust)
diffusion-scale estimate \\\tilde{\sigma}^2\\ by estimating the long-run
variance (LRV) of standardized growth increments using AR(1)
pre-whitening and a Bartlett (Newey-West) HAC estimator, followed by
recoloring.

## Arguments

- mu:

  Numeric scalar. Estimated growth rate \\\hat{\mu}\\ used to center
  increments.

- delta_log_n:

  Numeric vector of length `q`. Log-scale increments \\\Delta Y_i\\.

- tau:

  Numeric vector of length `q`. Sampling intervals \\\tau_i \> 0\\.

## Value

A list with components:

- `sigma2_tilde`:

  Numeric scalar. OEAR diffusion-scale estimate \\\tilde{\sigma}^2\\.

- `rho_tilde_pw`:

  Numeric scalar. Estimated AR(1) pre-whitening coefficient
  \\\tilde\rho\_{\mathrm{pw}}\\.

- `j`:

  Integer. Selected Bartlett truncation lag \$J\$.

## Details

The estimator targets the LRV of the standardized centered increment
series \$\$U_i = (\Delta Y_i - \hat{\mu}\\\tau_i)/\sqrt{\tau_i}, \qquad
i=1,\ldots,q,\$\$ which is interpreted as an effective environmental
variance per unit time.

Implementation steps:

1.  Center \\u_i\\ by subtracting \\\bar u\\ to obtain \\\tilde u_i\\.

2.  Estimate the AR(1) pre-whitening coefficient
    \\\tilde\rho\_{\mathrm{pw}}\\ by OLS in \\\tilde u_i =
    \rho\_{\mathrm{pw}}\\\tilde u\_{i-1} + \varepsilon_i\\ and form
    pre-whitened residuals \\\tilde\varepsilon_i\\.

3.  Compute residual autocovariances and the Bartlett HAC LRV
    \\\widetilde{\mathcal C}^{(\varepsilon)}\_{\mathrm{NW}}(J)\\.

4.  Choose the truncation lag \\J\\ via the Andrews (1991) AR(1) plug-in
    rule specialized to the Bartlett window.

5.  Recolor to the original scale to obtain \$\$\tilde{\sigma}^2 =
    \widetilde{\mathcal C}\_{\mathrm{NW}}(J) = \widetilde{\mathcal
    C}^{(\varepsilon)}\_{\mathrm{NW}}(J)/
    (1-\tilde\rho\_{\mathrm{pw}})^2.\$\$

This construction is designed to be robust to short-run autocovariance
and, under additive observation error on the log scale, to the
cancellation property of the LRV when increments are appropriately
centered.

## References

Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
consistent covariance matrix estimation. *Econometrica*, 59(3), 817–858.

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
heteroskedasticity and autocorrelation consistent covariance matrix.
*Econometrica*, 55(3), 703–708.
