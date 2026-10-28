##' Integrate the leaf's thermal-damage states over a temperature series
##'
##' The leaf's lasting damage model, run on its own: three fractions of leaf
##' area, functional ⇌ impaired → damaged, plus an acclimation state, driven by a
##' leaf-temperature time series. This is the route for fitting the layer's
##' parameters to a leaf experiment without a whole-plant model; the comparison
##' with measurements (Fv/Fm, say) belongs to the caller.
##'
##' The dynamics, per day, with \eqn{T} the leaf temperature and \eqn{F = 1 - I - D}:
##'
##' \preformatted{
##'   dA/dt = acclim_alpha * softplus(T - acclim_T; acclim_s) - acclim_beta * A
##'   dI/dt = k_imp(T) * F - (k_rec(T) + k_mat) * I
##'   dD/dt = k_mat * I
##'   k_imp(T) = k_i * exp[(k_i_Ea/R)(1/T_ref - 1/T)] * K/(1 + K),
##'   K = exp[(impair_H/R)(1/Tmid - 1/T)],  Tmid shifted by the acclimation
##' }
##'
##' Every parameter is a field of `x` (`$k_i_`, `$impair_Tmid_`, `$acclim_beta_`,
##' ...); set them on the leaf before calling. Neither the gate
##' `use_thermal_damage_` nor the drivers are read. The rates themselves are
##' `x$damage_rates(leaf_temp, A, I, D)`.
##'
##' Temperature is linear between the given times. The integrator is fixed-step
##' RK4 with each interval cut so the fastest rate times the step stays below 0.1
##' and the temperature moves at most 0.25 K per step, so a coarse output grid does
##' not coarsen the solution.
##'
##' @param x a `Leaf`, from [leaf_model()]
##' @param times days, strictly increasing
##' @param leaf_temp leaf temperature at each time, deg C
##' @param A0,I0,D0 the acclimation state and the impaired and damaged fractions
##'   at `times[1]`
##' @return A data.frame with `time`, `leaf_temp`, `A`, `I`, `D` and the
##'   functional fraction `F = 1 - I - D`.
##' @examples
##' l <- leaf_model()
##' # a 6-hour 46 C heatwave, then two days at 25 C
##' t <- seq(0, 2.25, by = 1 / 24)
##' temp <- ifelse(t < 0.25, 46, 25)
##' d <- leaf_damage_trajectory(l, t, temp)
##' tail(d, 3)
##' @export
leaf_damage_trajectory <- function(x, times, leaf_temp, A0 = 0, I0 = 0, D0 = 0) {
  if (!inherits(x, "Leaf")) {
    stop("`x` must be a Leaf, from leaf_model()", call. = FALSE)
  }
  times <- as.numeric(times)
  leaf_temp <- as.numeric(leaf_temp)
  if (length(leaf_temp) == 1L) {
    leaf_temp <- rep(leaf_temp, length(times))
  }
  y <- x$damage_trajectory(times, leaf_temp, A0, I0, D0)
  n <- length(times)
  out <- data.frame(time = times, leaf_temp = leaf_temp,
                    A = y[seq_len(n)], I = y[n + seq_len(n)],
                    D = y[2L * n + seq_len(n)])
  out$F <- 1 - out$I - out$D
  out
}
