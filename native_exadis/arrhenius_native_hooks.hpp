#pragma once

// Native Arrhenius/TST hook interface scaffold for ExaDiS/ParaDiS integration.
//
// This header is deliberately self-contained.  It does not include ExaDiS
// internal headers yet because the Taylor_DDD repository is an adapter/testing
// repository rather than a vendored ExaDiS source tree.  The purpose is to fix
// the physics interface and audit outputs before these classes are moved into
// the native ExaDiS mobility, topology, cross-slip, and collision modules.
//
// Non-negotiable convention:
//   activated mechanisms are represented by hazards
//       R = eta0 exp[-G(tau_eff,T)/(kB T)]
//       P(dt) = 1 - exp[-R dt]
//   and the stress argument is the event-conjugate local stress.  If the event
//   is naturally force-work based, use
//       tau_eff = F_event * x_dagger / v_star
//   and do not use a diagnostic average stress as the kinetic input.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace arrhenius_exadis {

constexpr double kB_eV_per_K = 8.617333262145e-5;
constexpr double eV_J = 1.602176634e-19;

struct ExpFloorBarrierParams {
    double H_eV = 0.50;
    double S_kB = -9.0;
    double sigma_c_Pa = 14.5e9;
    double floor_fraction = 0.20;
    double a = 6.65607;
    double n = 2.15276;
    double eta0_s = 1.0e12;

    void validate() const {
        if (sigma_c_Pa <= 0.0) throw std::runtime_error("sigma_c_Pa must be positive");
        if (floor_fraction < 0.0 || floor_fraction > 1.0) throw std::runtime_error("floor_fraction must be in [0,1]");
        if (a <= 0.0 || n <= 0.0) throw std::runtime_error("EXP-floor shape parameters must be positive");
        if (eta0_s <= 0.0) throw std::runtime_error("eta0_s must be positive");
    }
};

inline double G0_eV(const ExpFloorBarrierParams& p, double T_K) {
    return p.H_eV - kB_eV_per_K * T_K * p.S_kB;
}

inline double exp_floor_barrier_eV(const ExpFloorBarrierParams& p, double tau_eff_Pa, double T_K) {
    p.validate();
    const double G0 = G0_eV(p, T_K);
    const double floor = p.floor_fraction * G0;
    const double x = std::max(tau_eff_Pa, 0.0) / p.sigma_c_Pa;
    const double shape = p.floor_fraction + (1.0 - p.floor_fraction) * std::exp(-p.a * std::pow(x, p.n));
    return std::max({0.0, floor, G0 * shape});
}

inline double arrhenius_rate_s(const ExpFloorBarrierParams& p, double tau_eff_Pa, double T_K) {
    const double G = exp_floor_barrier_eV(p, tau_eff_Pa, T_K);
    return p.eta0_s * std::exp(-G / (kB_eV_per_K * T_K));
}

inline double event_probability(const ExpFloorBarrierParams& p, double tau_eff_Pa, double T_K, double dt_s) {
    const double rdt = std::clamp(arrhenius_rate_s(p, tau_eff_Pa, T_K) * dt_s, 0.0, 50.0);
    return -std::expm1(-rdt);
}

inline double force_work_tau_eff_Pa(double F_event_N, double x_dagger_m, double v_star_m3) {
    if (v_star_m3 <= 0.0) throw std::runtime_error("v_star_m3 must be positive");
    return F_event_N * x_dagger_m / v_star_m3;
}

inline double signed_forward_minus_reverse_rate_s(const ExpFloorBarrierParams& p, double tau_eff_Pa, double T_K) {
    return arrhenius_rate_s(p, std::max(tau_eff_Pa, 0.0), T_K)
         - arrhenius_rate_s(p, std::max(-tau_eff_Pa, 0.0), T_K);
}

enum class MechanismKind {
    PeierlsGlide,
    ForestDepinning,
    JunctionZip,
    JunctionUnzip,
    CrossSlip,
    ActivatedCollision
};

struct MechanismAuditRecord {
    std::uint64_t step = 0;
    double time_s = 0.0;
    MechanismKind mechanism = MechanismKind::ForestDepinning;
    double tau_local_Pa = 0.0;
    double tau_eff_Pa = 0.0;
    double force_N = 0.0;
    double x_dagger_m = 0.0;
    double v_star_m3 = 0.0;
    double barrier_eV = 0.0;
    double rate_s = 0.0;
    double probability_dt = 0.0;
    bool selected = false;
    bool deterministic_geometry_only = false;
    std::string note;
};

// ArrheniusMobilityLaw
// --------------------
// Native target for replacing segment glide mobility with a signed
// forward-minus-reverse Peierls hazard.  The ExaDiS mobility law should call
// signed_velocity_increment() using the resolved glide stress.  It should not
// inject a deterministic Peierls threshold.
class ArrheniusMobilityLaw {
public:
    explicit ArrheniusMobilityLaw(ExpFloorBarrierParams peierls_params, double glide_jump_m)
        : peierls_params_(peierls_params), glide_jump_m_(glide_jump_m) {}

    double signed_rate_s(double tau_resolved_Pa, double T_K) const {
        return signed_forward_minus_reverse_rate_s(peierls_params_, tau_resolved_Pa, T_K);
    }

    double signed_velocity_m_s(double tau_resolved_Pa, double T_K) const {
        return glide_jump_m_ * signed_rate_s(tau_resolved_Pa, T_K);
    }

    MechanismAuditRecord audit(std::uint64_t step, double time_s, double tau_resolved_Pa, double T_K, double dt_s) const {
        MechanismAuditRecord r;
        r.step = step;
        r.time_s = time_s;
        r.mechanism = MechanismKind::PeierlsGlide;
        r.tau_local_Pa = tau_resolved_Pa;
        r.tau_eff_Pa = tau_resolved_Pa;
        r.barrier_eV = exp_floor_barrier_eV(peierls_params_, std::abs(tau_resolved_Pa), T_K);
        r.rate_s = std::abs(signed_rate_s(tau_resolved_Pa, T_K));
        r.probability_dt = std::clamp(r.rate_s * dt_s, 0.0, 50.0);
        return r;
    }

private:
    ExpFloorBarrierParams peierls_params_;
    double glide_jump_m_ = 0.0;
};

// ArrheniusTopology
// -----------------
// Native target for topology operations that are currently deterministic once
// geometry/force criteria are met.  Geometry should generate candidates;
// kinetics must be hazard based.
class ArrheniusTopology {
public:
    ArrheniusTopology(ExpFloorBarrierParams depinning_params,
                      ExpFloorBarrierParams junction_zip_params,
                      ExpFloorBarrierParams junction_unzip_params)
        : depinning_params_(depinning_params),
          junction_zip_params_(junction_zip_params),
          junction_unzip_params_(junction_unzip_params) {}

    double depinning_probability(double F_PK_N, double x_dagger_m, double v_star_m3,
                                 double T_K, double dt_s) const {
        const double tau_eff = force_work_tau_eff_Pa(F_PK_N, x_dagger_m, v_star_m3);
        return event_probability(depinning_params_, tau_eff, T_K, dt_s);
    }

    double junction_reaction_probability(double delta_work_J, double v_star_m3,
                                         double T_K, double dt_s, bool unzip) const {
        const double tau_eff = std::max(delta_work_J, 0.0) / v_star_m3;
        return event_probability(unzip ? junction_unzip_params_ : junction_zip_params_, tau_eff, T_K, dt_s);
    }

private:
    ExpFloorBarrierParams depinning_params_;
    ExpFloorBarrierParams junction_zip_params_;
    ExpFloorBarrierParams junction_unzip_params_;
};

// ArrheniusCrossSlip
// ------------------
// Native target for force-based cross slip.  Candidate planes remain generated
// by ExaDiS crystallography.  Selection is a competing hazard using the local
// driving difference between primary and cross-slip planes.
class ArrheniusCrossSlip {
public:
    explicit ArrheniusCrossSlip(ExpFloorBarrierParams cross_slip_params)
        : cross_slip_params_(cross_slip_params) {}

    double probability(double tau_primary_Pa, double tau_cross_Pa, double T_K, double dt_s) const {
        const double tau_eff = tau_cross_Pa - tau_primary_Pa;
        return event_probability(cross_slip_params_, tau_eff, T_K, dt_s);
    }

private:
    ExpFloorBarrierParams cross_slip_params_;
};

// ArrheniusCollision
// ------------------
// Native target for collisions only when the reaction is activated.  Pure core
// overlap or numerical collision cleanup should remain deterministic and should
// be audited as deterministic_geometry_only.
class ArrheniusCollision {
public:
    explicit ArrheniusCollision(ExpFloorBarrierParams activated_collision_params)
        : activated_collision_params_(activated_collision_params) {}

    double activated_probability(double reaction_force_N, double x_dagger_m,
                                 double v_star_m3, double T_K, double dt_s) const {
        const double tau_eff = force_work_tau_eff_Pa(reaction_force_N, x_dagger_m, v_star_m3);
        return event_probability(activated_collision_params_, tau_eff, T_K, dt_s);
    }

private:
    ExpFloorBarrierParams activated_collision_params_;
};

} // namespace arrhenius_exadis
