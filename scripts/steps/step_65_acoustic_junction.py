#!/usr/bin/env python3
"""AUD-4 / T4.2 -- acoustic-side sensitivity of the junction z_j.

Question (plan AUD-4 condition b): does the acoustic sector bound z_j
from below?  Under the piecewise ambient A(eta) (step_64) the power-law
tail covers z > z_j, so a naive FLRW-image reading would place
recombination physics inside the tail and catastrophically shift r_s /
theta_s.  This step resolves which reading is correct inside the
framework's own rules:

1. Deep-ambient well nesting.  Re-run the step_58 well solve with
   u(r_max) = u_amb > 0 and measure the local depth u_loc = u_c - u_amb.
   The e^{-u} starvation factor makes the interior source shut off as
   the ambient deepens.

2. Depth decomposition at fixed observed redshift.  Rule 9: the
   observed 1+z is the endpoint clock ratio, u_emit = ln(1+z).  A source
   observed at z = 1100 sits at total depth u ~ 7 regardless of the
   ambient map; only the ambient/local split changes.  Under the tail,
   u_amb(eta ~ z_rec) is deep, so the split is ~99% ambient.

3. Acoustic insensitivity.  r_s in matter units is relational inside
   the emitter's environment: an overall uniform clock shift rescales
   every rate equally.  The only ambient-level dependence enters the
   well STRUCTURE (via starvation of the internal gradient structure),
   which is the transport/decoration level, not the sound-horizon
   level.  The FLRW-image reading -- where eta_lb(z) maps to the
   observed distance-redshift relation directly -- is exactly the
   pre-derivation ansatz the Projection Dictionary (plan Part VI)
   retires: distances come from the constraint slice / spatial
   geometry, not from the ambient bookkeeping map.

Outputs results/step_65_acoustic_junction.json
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step_58_well_family as w  # noqa: E402


def solve_well_amb(R, sector, u_amb, nr=4000, rmax=40.0):
    """step_58 solve_well with nonzero ambient boundary u(rmax)=u_amb."""
    r = np.exp(np.linspace(np.log(R * 1e-4), np.log(rmax), nr))
    dr = np.diff(r)
    shape = np.where(r < R, 3.0 / R ** 3, 0.0)
    u = np.full(nr, u_amb)
    it = 0
    for it in range(400):
        src = shape * np.exp(-u)
        F = np.concatenate([[0.0], np.cumsum(
            0.5 * (src[:-1] * r[:-1] ** 2 + src[1:] * r[1:] ** 2) * dr)])
        g = np.maximum(F, 0.0) / np.maximum(r ** 2, 1e-30)
        q = np.array([w.invert_flux(gi, sector) for gi in g])
        integ = 0.5 * (q[:-1] + q[1:]) * dr
        u_new = u_amb + np.concatenate((np.cumsum(integ[::-1])[::-1],
                                        [0.0]))
        du = np.max(np.abs(u_new - u))
        u = u_new
        if du < 1e-10:
            break
    return float(u[0]), it


def main():
    out = {"question": "does the acoustic sector bound z_j from below?",
           "sector": "two_branch", "k": w.K}

    # 1. starvation law: u_loc vs u_amb for a mid-compactness well
    scan = {}
    for R in (2.0, 0.5, 0.05):
        rows = []
        for ua in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0):
            uc, it = solve_well_amb(R, "two_branch", ua)
            rows.append({"u_amb": ua, "u_c": uc,
                         "u_loc": uc - ua, "iters": it})
        arr = np.array([r["u_loc"] for r in rows])
        uav = np.array([r["u_amb"] for r in rows])
        m = uav > 0
        slope = float(np.polyfit(uav[m], np.log(arr[m]), 1)[0])
        scan["R=%.2f" % R] = {"rows": rows,
                              "starvation_exponent": slope}
        out["deep_ambient_scan_" + ("R=%.2f" % R)] = scan["R=%.2f" % R]

    # 2. depth decomposition for emitters at observed redshift z.
    #    u_emit = ln(1+z) is fixed by the endpoint ratio; solve for the
    #    ambient level the source would sit on if u_loc follows the
    #    measured starvation law u_loc(u_amb).
    decomp = {}
    R_ref = "R=0.50"
    rows = scan["R=0.50"]["rows"]
    uav = np.array([r["u_amb"] for r in rows])
    ulv = np.array([r["u_loc"] for r in rows])
    m = uav > 0
    b, a = np.polyfit(uav[m], np.log(ulv[m]), 1)

    def u_loc(ua):
        return float(np.exp(a + b * ua))

    for z in (1100.0, 100.0, 10.0, 1.0):
        u_tot = np.log(1.0 + z)
        # solve u_amb + u_loc(u_amb) = u_tot
        lo, hi = 0.0, u_tot
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if mid + u_loc(mid) < u_tot:
                lo = mid
            else:
                hi = mid
        ua = 0.5 * (lo + hi)
        decomp["z=%g" % z] = {
            "u_total": float(u_tot), "u_amb_inferred": float(ua),
            "u_loc_inferred": float(u_loc(ua)),
            "ambient_fraction": float(ua / u_tot),
            "note": "well reference R=0.5 r* (galactic-scale); "
                    "starvation law u_loc = %.3f exp(%.3f u_amb)"
                    % (np.exp(a), b)}
    out["depth_decomposition"] = decomp
    out["starvation_law_fit"] = {"u_loc0": float(np.exp(a)),
                                 "decay_per_u_amb": float(b),
                                 "reference": R_ref}

    # 3. verdict
    out["acoustic_bound"] = {
        "verdict": "no lower bound on z_j from the acoustic sector",
        "reasoning": [
            "observed redshift = endpoint clock ratio (Rule 9): "
            "u_emit = ln(1+z) independent of the ambient map shape",
            "r_s in matter units is relational inside the emitter "
            "environment; a uniform clock shift rescales all rates "
            "equally and leaves r_s invariant",
            "the tail's eta_lb(z) divergence (step_64: 3.6e5x at "
            "z=1100) is a bookkeeping map, not the observed "
            "distance-redshift relation; distances are set by the "
            "constraint-slice geometry, not the ambient image",
            "ambient level DOES modulate well interiors via "
            "starvation (u_loc ~ e^{%.2f u_amb}) -- sources on deep "
            "ambient are ~99%% ambient-depth -- but that is a "
            "transport/structure effect, not an r_s shift" % b],
        "flrw_image_reading": "would fail catastrophically "
            "(eta_lb ratio 1.7x at z=2 growing to 3.6e5x at z=1100) "
            "-- that reading is the pre-derivation ansatz the "
            "Projection Dictionary retires",
        "binding_constraint": "regularity only: z_j <= 0.687 "
            "(step_64); physical z_j is where the well network hands "
            "the ambient its drift",
        "remaining_dependency": "well-frame r_s and peak morphology "
            "under deep-ambient nesting = AUD-6 / T6.3 "
            "(constraint-slice machinery); this step establishes only "
            "that the ambient MAP does not supply the bound"}

    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_65_acoustic_junction.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["acoustic_bound"], indent=1))
    print("starvation law: u_loc ~ %.3f exp(%.3f u_amb)"
          % (np.exp(a), b))
    for k, v in decomp.items():
        print("  %s: u_tot=%.3f u_amb=%.3f u_loc=%.4f (ambient %.1f%%)"
              % (k, v["u_total"], v["u_amb_inferred"],
                 v["u_loc_inferred"], 100 * v["ambient_fraction"]))
    print("wrote", dest)


if __name__ == "__main__":
    main()
