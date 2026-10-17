#!/usr/bin/env python3
"""AUD-3 / T3.1 (second pass) — in-well X_env volume weighting of the
growth vertex.

Step_62 established that inside collapsed regions the nonlinear pair
map x_env = (1-y)/y runs large (S_Sigma -> 0, G_eff -> 1): collapsed
matter grows Newtonianly.  The residual f sigma8 excess is therefore
carried only by the mode-weighted AMBIENT fraction of the perturbation
support.  This step computes a collapsed-fraction growth benchmark and a separate
environment-density response diagnostic.

The weighting is derived, not fit.  The mode support of a linear
perturbation at scale k is the full density field; the collapsed
fraction f_coll(z) removes the fraction of that support that has
virialised.  Standard Press-Schechter bookkeeping gives

    f_coll(z) = erfc( delta_c / (sqrt(2) sigma_M(z)) )

with sigma_M(z) = sigma_M(0) D_LCDM(z)/D_LCDM(0) the rms overdensity
on the threshold mass scale.  Two weightings are reported:

  (a) mass-weighted: collapsed fraction at the M ~ 10^12 M_sun
      nonlinearity boundary (sigma ~ 1 at z=0) — the conservative
      floor, since underdense regions also partially screen;
  (b) a lower threshold M ~ 10^10 M_sun (legacy output label "field").
      This is a collapsed-mass-threshold sensitivity calculation. The
      separate lognormal environment average is reported as a diagnostic,
      and is not substituted into the growth ODE.

In both weightings the ambient-mode response retains the two-branch
G_eff(z) while the collapsed fraction contributes G_eff = 1.  The
convolution G_eff,mode = (1 - w_coll) [1 + 2 b_A^2 S_amb] + w_coll
then re-enters the step_62 growth ODE.

Outputs results/step_74_fs8_inwell_convolution.json
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import lognorm, norm

OM, OL = 0.3, 0.7
S8_PLANCK = 0.834
BETA_A = -1.0
K = 16.03
DELTA_C = 1.686


def PX(xi, sector):
    xi = np.maximum(np.asarray(xi, dtype=float), 0.0)
    if sector == "lcdm":
        return np.inf * np.ones_like(xi)
    if sector == "baseline":
        return 1.0 + 2.0 * xi
    if sector == "two_branch":
        return K * np.sqrt(xi) + 2.0 * xi
    if sector == "exp_interp":
        return 1.0 - np.exp(-K * np.sqrt(xi)) + 2.0 * xi
    raise ValueError(sector)


def Ez(z):
    return np.sqrt(OM * (1.0 + z) ** 3 + OL)


# --- LCDM growth for the variance normalisation ---------------------
def D_lcdm(z):
    """Unnormalised LCDM growth factor via integral form."""
    def rhs(lna, y):
        a = np.exp(lna)
        zz = 1.0 / a - 1.0
        dlnH = -(1.5 * OM * (1.0 + zz) ** 3) / Ez(zz) ** 2
        src = 1.5 * OM * (1.0 + zz) ** 3 / Ez(zz) ** 2
        return [y[1], src * y[0] - (2.0 + dlnH) * y[1]]
    ai = 1.0 / 201.0
    sol = solve_ivp(rhs, [np.log(ai), 0.0], [ai, ai],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    a = 1.0 / (1.0 + np.asarray(z, dtype=float))
    D = sol.sol(np.log(a))[0]
    return D / sol.sol([0.0])[0][0]


# --- collapsed fraction ----------------------------------------------
def sigma_m_z0(m_over_msun):
    """Approximate sigma(M,0) on galactic-to-cluster scales: power-law
    anchored at sigma_8 (R8 -> M8 ~ 4 pi rho_m R8^3/3 ~ 5.9e14 Msun
    for Om=0.3, h=0.7).  Slope n_eff ~ -2.0 -> sigma ~ M^{-(3+n_eff)/6}
    ~ M^{-1/6}...  standard CDM local index gives sigma ~ M^{-0.27}
    near 1e12 Msun; we use the fitted form sigma_8 * (M/M8)^{-0.28}."""
    M8 = 5.9e14
    return S8_PLANCK * (m_over_msun / M8) ** (-0.28)


def f_collapsed(z, m_threshold):
    """Press-Schechter fraction of mass in structures above
    m_threshold at redshift z."""
    Dz = D_lcdm(np.atleast_1d(z))
    sig0 = sigma_m_z0(m_threshold)
    sig_z = sig0 * Dz
    nu = DELTA_C / sig_z
    return 2.0 * norm.sf(nu)  # erfc(nu/sqrt(2))


# --- field-weighted environment averaging ---------------------------
def xi_spatial_from_g(g):
    """xi = (grad phi)^2/(2 (H0/c)^2), grad phi ~ g/c^2 (beta_A=-1)."""
    H0_SI, C_L = 2.27e-18, 3.0e8
    return (g / C_L**2 / (H0_SI / C_L)) ** 2 / 2.0


def density_pdf(sigma_linear):
    """Unit-mean density ratio with variance sigma_linear**2."""
    variance_log = np.log1p(float(sigma_linear)**2)
    return lognorm(np.sqrt(variance_log), scale=np.exp(-0.5 * variance_log))


def env_weighted_S(z, sector):
    """Mode-weighted suppression: the perturbation support is the full
    density field.  Approximate the local scalar acceleration by
    delta -> g ~ delta * g_ambient (proportional environment depth),
    so in-well regions (delta >> 1, g > g_t) carry S_Sigma -> 0.
    Lognormal PDF matched to sigma_lin(z) = sigma_8 D(z)."""
    sig = 0.834 * float(D_lcdm([z])[0])
    if sig < 1e-3:
        sig = 1e-3
    # lognormal with <delta>=0, variance sig^2 -> ln delta field
    dist = density_pdf(sig)  # rho/mean = 1 + delta
    # sample the PDF
    rv = dist.rvs(size=20000, random_state=7)
    delta = rv - 1.0
    # environment depth ~ delta sets the local screening kinetic
    # variable; g_delta ~ delta * g_t at disk-mean (u ~ 0.18)
    g_t = 3.4e-10
    g_local = np.maximum(delta, 0.0) * 0.18 * g_t
    xi_loc = xi_spatial_from_g(g_local)
    xi_tot = Ez(z)**2 / 2.0 + np.where(xi_loc > 0, xi_loc, 0.0)
    S = 1.0 / PX(xi_tot, sector)
    return float(np.mean(S))


def growth_diluted(sector, weighting):
    """Growth ODE with G_eff = 1 + 2 b_A^2 S_amb * (1 - w_coll)."""
    ai = 1.0 / 201.0

    def w_coll(z):
        if weighting == "mass":
            return float(f_collapsed(z, 1e12)[0])
        if weighting == "field":
            return float(f_collapsed(z, 1e10)[0])
        raise ValueError(weighting)

    def rhs(lna, y):
        a = np.exp(lna)
        z = 1.0 / a - 1.0
        dlnH = -(1.5 * OM * (1.0 + z) ** 3) / Ez(z) ** 2
        S_amb = 1.0 / PX(Ez(z)**2 / 2.0, sector)
        G_eff = ((1.0 - w_coll(z)) * (1.0 + 2.0 * BETA_A**2 * S_amb)
                 + w_coll(z))
        src = 1.5 * OM * (1.0 + z) ** 3 / Ez(z) ** 2 * G_eff
        return [y[1], src * y[0] - (2.0 + dlnH) * y[1]]

    sol = solve_ivp(rhs, [np.log(ai), 0.0], [ai, ai],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    return sol


def main():
    out = {"step": "step_74_fs8_inwell_convolution",
           "purpose": ("evaluate the residual identified in step_62: "
                       "convolve the ambient X_env response over the "
                       "mode-weighted collapsed fraction; in-well "
                       "regions revert to G_eff = 1")}

    zout = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])

    # collapsed fractions
    out["collapsed_fraction"] = {
        "model": "Press-Schechter, delta_c = 1.686, "
                 "sigma_M(0) = sigma8 (M/M8)^{-0.28}",
        "z": zout.tolist(),
        "f_mass_above_1e12": f_collapsed(zout, 1e12).tolist(),
        "f_mass_above_1e10": f_collapsed(zout, 1e10).tolist(),
        "reading": ("The arrays are Press-Schechter collapsed mass fractions "
                    "at two threshold masses; they approximate the "
                    "screened fraction used by the growth benchmark.")
    }

    rd = json.load(open(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "..", "..", "results",
        "step_62_growth_fsigma8.json")))["reference_data"]

    out["weighting_labels"] = {"mass": "collapsed mass above 1e12 Msun",
                               "field": "legacy label: collapsed mass above 1e10 Msun"}
    out["normalization"] = "Present-day sigma8 fixed to 0.834 for every sector; S8_0 is assigned, not independently predicted."
    for weighting in ("mass", "field"):
        res = {}
        for sector in ("baseline", "two_branch"):
            sol = growth_diluted(sector, weighting)
            a = 1.0 / (1.0 + zout)
            D = sol.sol(np.log(a))[0]
            f = sol.sol(np.log(a))[1] / D
            D0 = sol.sol([0.0])[0][0]
            s8 = S8_PLANCK * D / D0
            fs8 = f * s8
            fs8_pts = np.interp([p["z"] for p in rd["fs8_points"]],
                                zout, fs8)
            chi2_fs8 = float(np.sum(
                [(fi - p["fs8"])**2 / p["err"]**2
                 for fi, p in zip(fs8_pts, rd["fs8_points"])]))
            S8v = s8 * np.sqrt(OM / 0.3)
            chi2_s8 = float(np.sum(
                [(S8v[0] - p["S8"])**2 / p["err"]**2
                 for p in rd["S8_points"]]))
            res[sector] = {
                "fs8": fs8.tolist(),
                "S8_0": float(S8v[0]),
                "fs8_at_0.61": float(np.interp(0.61, zout, fs8)),
                "chi2_total": chi2_fs8 + chi2_s8,
            }
        out[f"diluted_{weighting}"] = res

    # environment-resolved ambient suppression at the data redshifts
    out["env_resolved_S_amb"] = {
        "note": ("Separate diagnostic, not an ODE input: unit-mean lognormal density ratio with variance sigma_lin(z)^2; "
                 "in-well deltas screen through xi_spatial, ambient "
                 "modes keep the temporal xi_cosmo"),
        "z": [0.0, 0.5, 1.0],
        "two_branch": [env_weighted_S(z, "two_branch") for z in
                       (0.0, 0.5, 1.0)],
        "naive_S": [1.0 / PX(Ez(z)**2 / 2.0, "two_branch")
                    for z in (0.0, 0.5, 1.0)],
    }

    out["verdict"] = {
        "reading": ("The growth benchmark recomputes the response with two "
                    "collapsed-mass thresholds. The environment-density "
                    "average is reported separately. Compare computed "
                    "chi2 values using the declared common normalization."),
        "status": "collapsed-fraction growth benchmark computed",
    }

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "results",
                        "step_74_fs8_inwell_convolution.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    for w in ("mass", "field"):
        for s in ("baseline", "two_branch"):
            r = out[f"diluted_{w}"][s]
            print(f"{w} {s}: S8={r['S8_0']:.3f} "
                  f"fs8(0.61)={r['fs8_at_0.61']:.3f} "
                  f"chi2={r['chi2_total']:.1f}")
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
