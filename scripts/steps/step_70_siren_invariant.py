#!/usr/bin/env python3
"""Standard-siren invariant check: is the source binary locally standard?

The standard-siren mapping of step_21 (M_det = M_tilde/(1+z)) assumes the
source binary radiates with Einstein-frame mass M_e = A_e M_tilde. A
natural objection: if the source's local physics is entirely standard,
the standard (1+z) mapping should return. This step settles it with the
dimensionless invariant

    I = G_loc M_chirp f / c^3      (a fixed point of the inspiral phase)

computed in each observer's own locally measured constants.

Source-local inference. A comoving observer at the source measures the
chirp on matter clocks (f_tilde_e), with a locally-measured Cavendish
constant G_loc,e = G_* A_e^2 (1 + alpha_e^2) — the bare coupling dressed
by the conformal factor and the ambient self-response. Matching the
standard template to the true dynamics (driven by G_dyn = G_*(1 +
alpha_A alpha_B) on Einstein-frame mass M_e = A_e M_tilde) returns

    M_loc = M_tilde (1 + alpha_A alpha_B) / (1 + alpha_e^2).

For screened systems (alpha_i, alpha_e -> 0), M_loc = M_tilde: the
source IS a standard binary in its own units — no varying-G anomaly is
needed at emission. The (1+z)^{-1} mass mapping is therefore not a
statement about nonstandard local physics; it is a transport statement:
on the static background the coordinate frequency is conserved
(f_o = f_e) rather than redshifted, so the standard
M_det = M_src(1+z) never applies.

Detector side. I_det = G_loc,o M_det f_tilde_o with
G_loc,o = G_* A_o^2 (1 + alpha_o^2). The full ratio

    I_det/I_src = A_o^2 (1 + alpha_o^2) / (1 + alpha_A alpha_B)

equals 1 today (A_o = 1, screened alphas): the invariant is preserved
end-to-end — the mass mapping is internally consistent.

Exposed constraint surface. The same formula implies the locally
measured gravitational 'constant' varies with the ambient field,
G_loc = G_* A^2(phi) (1 + alpha_amb^2), giving

    Gdot/G = 2 (A-dot/A) + d ln(1 + alpha_amb^2)/dt
           ~= 2 (A-dot/A) evaluated at the system.

The evolution is dominated by the conformal factor itself — the
clock-amplitude channel, through which an ambient drift reaches
interior conformal factors essentially fully on the leaking branch
(transmission T ~ 1, step_76). The screened self-response
alpha_amb ~ alpha_0 S_Sigma,^sun ~ 1e-7 suppresses only the
shear-squared correction, not the A^2 drift. Published LLR analyses
bound |Gdot/G| <~ 4e-13 /yr (Williams et al.), hence |A-dot/A| <~
2e-13 /yr at the Earth--Moon location today: the local ambient drift
is pinned to ~0.3% of the Hubble rate. The siren mapping thus stands
on a verified invariant, with the epoch-dependence of G_loc at the
source end identified as the physical content of the assumption — the
item to confront in the LVK mass-spectrum channel rather than a hidden
defect.
"""
import json
import os
import sympy as sp


def run():
    G, Ao, Ae = sp.symbols('G_* A_o A_e', positive=True)
    a_o, a_e, aA, aB = sp.symbols('a_o a_e a_A a_B', positive=True)
    Mt, ft_e = sp.symbols('M_t f_te', positive=True)

    G_dyn = G * (1 + aA * aB)
    M_e = Ae * Mt
    G_loc_e = G * Ae**2 * (1 + a_e**2)
    # template match at source (see docstring): M_loc
    M_loc = sp.simplify(G_dyn * M_e * Ae / G_loc_e)
    I_src = sp.simplify(G_loc_e * M_loc * ft_e)

    G_loc_o = G * Ao**2 * (1 + a_o**2)
    M_det = Ao * Ae * Mt                     # step_21 result
    I_det = sp.simplify(G_loc_o * M_det * (Ae / Ao) * ft_e)

    ratio = sp.simplify(I_det / I_src)
    screened = sp.simplify(ratio.subs({aA: 0, aB: 0, a_e: 0, a_o: 0}))
    mloc_ratio = sp.simplify(M_loc / Mt)

    # --- numeric Ġ/G implication on the drift branch -------------------
    # Gdot/G ~ 2 (A-dot/A): the A^2 factor dominates; the (1+alpha^2)
    # shear correction carries the screened charge a_amb = alpha_0
    # S_Sigma_sun ~ sqrt(2)*9.5e-8 only in its own term.  Ambient
    # field-value drift transmits to interior clocks at T ~ 1 on the
    # leaking branch (step_76), so the LLR bound applies to the ambient
    # drift essentially directly.
    H0 = 70e3 / 3.085677581e22          # s^-1
    alpha_amb_sun = 2**0.5 * 9.5e-8     # screened solar self-response
    yr = 365.25 * 86400
    T_transmission = 0.99685            # step_76, LLR excursion scale
    gdot_over_G_at_hubble_drift = 2 * H0 * yr   # drift at the Hubble rate
    llr_bound = 4e-13                    # /yr (Williams et al.)
    adot_bound = llr_bound / 2.0         # |A-dot/A| bound at Earth-Moon
    ambient_drift_bound = adot_bound / T_transmission
    local_drift_over_H0 = ambient_drift_bound / (H0 * yr)

    out = {
        "invariant": "I = G_loc M_chirp f / c^3",
        "source_local_mass_ratio": str(mloc_ratio),
        "I_det_over_I_src": str(ratio),
        "screened_limit": str(screened),
        "source_is_locally_standard":
            "M_loc = M_tilde (1+a_A a_B)/(1+a_e^2) -> M_tilde for "
            "screened systems: the source IS standard in local units",
        "mapping_status":
            "M_det = M_tilde/(1+z) survives the invariant check — the "
            "nonstandard link is transport (f conserved on the static "
            "background), not local physics",
        "exposed_constraint": {
            "G_local": "G_loc = G_* A^2(1+alpha_amb^2) -> "
                       "Gdot/G ~= 2 (A-dot/A); the A^2 factor dominates, "
                       "the screened shear charge enters only its own "
                       "term",
            "excursion_transmission": T_transmission,
            "gdot_over_G_per_yr_at_hubble_drift":
                gdot_over_G_at_hubble_drift,
            "llr_bound_per_yr": llr_bound,
            "adot_bound_per_yr": adot_bound,
            "ambient_drift_bound_per_yr": ambient_drift_bound,
            "local_drift_fraction_of_H0": local_drift_over_H0,
            "note": "the LLR bound pins the ambient drift at the "
                    "Earth--Moon location to ~0.3% of the Hubble rate; "
                    "the epoch-dependence of G_loc at the source end is "
                    "the physical content of the siren mapping — "
                    "confronting it in the LVK mass-spectrum channel is "
                    "the remaining item",
        },
    }
    return out


if __name__ == "__main__":
    res = run()
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_70_siren_invariant.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    for k, v in res.items():
        print(f"{k}: {v}")
    print(f"Saved {path}")
