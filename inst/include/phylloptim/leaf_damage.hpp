// -*-c++-*-
#ifndef PHYLLOPTIM_LEAF_DAMAGE_HPP_
#define PHYLLOPTIM_LEAF_DAMAGE_HPP_

// The leaf's side of the ATLS thermal-damage layer, off unless
// use_thermal_damage_: its parameters, rates and trajectory, which read the
// dynamics in damage.hpp.
//
// Part of `Leaf`: include <phylloptim/leaf_model.hpp>, which includes every
// part in order.

#include <phylloptim/leaf.hpp>

namespace phylloptim {

// --- the ATLS thermal-damage layer ------------------------------------------


inline double Leaf::softplus_(double x, double s) {
  return damage::softplus(x, s);
}

// Inverts the peaked Arrhenius optimum Topt = H_d/(d_S - corr), corr =
// R*ln(Ea/(H_d - Ea)), for the d_S that puts it delta_topt higher. The optimum
// exists only for 0 < Ea < H_d, so outside it is refused, not a NaN capacity.
template <class T>
inline T Leaf::d_S_shifted(double Ea, double H_d, double d_S_base,
                           const T& delta_topt) const {
  if (!(Ea > 0.0 && H_d > Ea)) {
    util::stop("the thermal-damage layer moves the peaked Arrhenius optimum, "
               "which exists only for 0 < Ea < H_d; got Ea = " +
               util::format_double(Ea) + ", H_d = " + util::format_double(H_d));
  }
  const double corr = gas_constant * std::log(Ea / (H_d - Ea));
  const double Topt0_K = H_d / (d_S_base - corr);
  // Refused outside [250, 350] K rather than clamped there, as the impairment's
  // isokinetic temperature is (damage::k_imp).
  using odelia::util::to_passive;
  const T Topt_K = Topt0_K + delta_topt;
  if (to_passive(Topt_K) < 250.0 || to_passive(Topt_K) > 350.0) {
    util::stop("the thermal-damage layer moves the peaked Arrhenius optimum to " +
               util::format_double(to_passive(Topt_K)) + " K, outside [250, 350] K; "
               "topt_offset_ + dTopt_max_ * A/(A + K_A_) = " +
               util::format_double(to_passive(delta_topt)) + " K");
  }
  return H_d / Topt_K + corr;
}

// The layer's helpers on this leaf: its parameters, packed for the templates in
// `damage` (damage.hpp).
inline damage::dyn_pars<double> Leaf::damage_pars() const {
  damage::dyn_pars<double> p{};
  p[damage::dp_k_i] = k_i_;
  p[damage::dp_k_i_Ea] = k_i_Ea_;
  p[damage::dp_k_i_Tiso] = k_i_Tiso_;
  p[damage::dp_k_rec] = k_rec_;
  p[damage::dp_k_rec_Ea] = k_rec_Ea_;
  p[damage::dp_k_mat] = k_mat_;
  p[damage::dp_k_mat_Ea] = k_mat_Ea_;
  p[damage::dp_acclim_alpha] = acclim_alpha_;
  p[damage::dp_acclim_beta] = acclim_beta_;
  p[damage::dp_acclim_T] = acclim_T_;
  p[damage::dp_acclim_s] = acclim_s_;
  p[damage::dp_dTimpair_max] = dTimpair_max_;
  p[damage::dp_K_A] = K_A_;
  p[damage::dp_topt_offset] = topt_offset_;
  return p;
}

// A negative k_i_Ea_ is refused: impairment would slow as the leaf heats. Zero
// is allowed. A NaN passes, poisoning the outputs where a test can see it.
inline void Leaf::check_damage_pars() const {
  if (k_i_Ea_ < 0.0) {
    util::stop("the thermal-damage layer needs k_i_Ea_ >= 0, so that "
               "impairment does not slow with temperature; got k_i_Ea_ = " +
               util::format_double(k_i_Ea_));
  }
}

// Its own parameters only: reading jmax_H_d_ or jmax_d_S_ here would make
// recalibrating photosynthesis silently move the impairment rate.
inline double Leaf::k_imp_eff(double leaf_temp, double delta_impair) const {
  check_damage_pars();
  return damage::k_imp(leaf_temp, delta_impair, damage_pars());
}

inline double Leaf::delta_impair_at(double A) const {
  return damage::delta_impair(A, damage_pars());
}

inline double Leaf::k_rec_eff(double leaf_temp) const {
  return damage::k_rec_eff(leaf_temp, damage_pars());
}

inline std::vector<double> Leaf::damage_rates(double leaf_temp, double A,
                                              double I, double D) const {
  check_damage_pars();
  if (A <= -K_A_) {
    util::stop("damage_rates: the acclimation shift A/(A + K_A_) has its pole "
               "at A = -K_A_; got A = " + util::format_double(A) +
               ", K_A_ = " + util::format_double(K_A_));
  }
  const std::array<double, 3> r = damage::rates(leaf_temp, A, I, D, damage_pars());
  return {r[0], r[1], r[2]};
}

inline void Leaf::check_damage_initial_states(double A0, double I0,
                                              double D0) const {
  if (I0 < 0.0 || D0 < 0.0 || A0 <= -K_A_) {
    util::stop("damage_trajectory: the initial states need I0 >= 0, D0 >= 0 "
               "and A0 > -K_A_; got A0 = " + util::format_double(A0) +
               ", I0 = " + util::format_double(I0) + ", D0 = " +
               util::format_double(D0));
  }
}

inline void Leaf::check_damage_series(const std::vector<double>& times,
                                      const std::vector<double>& leaf_temp) const {
  const std::size_t n = times.size();
  if (n == 0 || leaf_temp.size() != n) {
    util::stop("damage_trajectory: times and leaf_temp must be non-empty and "
               "the same length; got " + std::to_string(n) + " and " +
               std::to_string(leaf_temp.size()));
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(times[i]) || !std::isfinite(leaf_temp[i]) ||
        (i > 0 && !(times[i] > times[i - 1]))) {
      util::stop("damage_trajectory: times must be finite and strictly "
                 "increasing, and leaf_temp finite; failed at index " +
                 std::to_string(i));
    }
  }
  check_damage_pars();
}

inline std::vector<double> Leaf::damage_trajectory(
    const std::vector<double>& times, const std::vector<double>& leaf_temp,
    double A0, double I0, double D0) const {
  check_damage_series(times, leaf_temp);
  check_damage_initial_states(A0, I0, D0);
  const std::size_t n = times.size();
  std::vector<double> out(3 * n);
  damage::trajectory<double>(times, leaf_temp, {A0, I0, D0}, damage_pars(),
                             [&](std::size_t i, const std::array<double, 3>& y) {
                               out[i] = y[0];
                               out[n + i] = y[1];
                               out[2 * n + i] = y[2];
                             });
  return out;
}

inline std::vector<std::string> Leaf::damage_par_names() const {
  return {"k_i_", "k_i_Ea_", "k_i_Tiso_", "k_rec_", "k_rec_Ea_", "k_mat_",
          "k_mat_Ea_", "acclim_alpha_", "acclim_beta_", "acclim_T_", "acclim_s_",
          "dTimpair_max_", "K_A_", "topt_offset_", "A0", "I0", "D0"};
}

// One forward-mode pass through the double trajectory's RK4, carrying every
// input's direction at once (a parameter or initial state), from which the
// requested columns are read: exact for the discrete map, at the double run's
// sub-step count.
inline std::vector<double> Leaf::damage_trajectory_gradient(
    const std::vector<double>& times, const std::vector<double>& leaf_temp,
    double A0, double I0, double D0, const std::vector<int>& which) const {
  check_damage_series(times, leaf_temp);
  check_damage_initial_states(A0, I0, D0);
  constexpr int n_inputs = damage::n_dyn_pars + 3;
  for (const int w : which) {
    if (w < 0 || w >= n_inputs) {
      util::stop("damage_trajectory_gradient: input index " + std::to_string(w) +
                 " is out of range; see damage_par_names()");
    }
  }
  const std::size_t n = times.size(), np = which.size();
  std::vector<double> out(3 * n * np);
  if (np == 0) {
    return out;
  }
  using T = xad::fwd<double, n_inputs>::active_type;
  const damage::dyn_pars<double> base = damage_pars();
  damage::dyn_pars<T> p;
  for (int j = 0; j < damage::n_dyn_pars; ++j) {
    p[std::size_t(j)] = base[std::size_t(j)];
    xad::derivative(p[std::size_t(j)])[std::size_t(j)] = 1.0;
  }
  std::array<T, 3> y0{T(A0), T(I0), T(D0)};
  for (int s = 0; s < 3; ++s) {
    xad::derivative(y0[std::size_t(s)])[std::size_t(damage::n_dyn_pars + s)] = 1.0;
  }
  damage::trajectory<T>(
      times, leaf_temp, y0, p, [&](std::size_t i, const std::array<T, 3>& y) {
        for (int st = 0; st < 3; ++st) {
          const auto& d = xad::derivative(y[std::size_t(st)]);
          for (std::size_t k = 0; k < np; ++k) {
            out[i + n * (std::size_t(st) + 3 * k)] = d[std::size_t(which[k])];
          }
        }
      });
  return out;
}

} // namespace phylloptim

#endif
