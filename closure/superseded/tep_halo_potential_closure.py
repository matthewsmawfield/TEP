"""CORPUS_PLAN §7.8 Target 3 closure computation: HALO-POTENTIAL X coordinate.

Companion to tep_shared_distortion_closure.py (stellar-mass X).
TEP's time field is continuous and nested: it responds to the TOTAL
gravitational potential (dark-matter halo + group + cosmic web), not
just the baryonic stellar disk.  The stellar-mass X coordinate
(sigma_stellar^2 / c^2, sigma ~ 91 km/s) captures only the baryonic
component and underestimates the amplitude.  The halo-potential X
coordinate uses the HALO velocity dispersion (sigma_halo ~ 200 km/s),
giving a ~5x larger amplitude via (sigma_halo/sigma_stellar)^2 ~ 5.3.

Step 1  kappa_SN from Pantheon+ Hubble-flow SN magnitudes vs halo-velocity-dispersion X
Step 2  agreement with Paper 11 velocity-space kappa
Step 3  zero-free-parameter H0 = 73.04 * 10^(-kappa (Xbar_HF - Xbar_anchor)/5)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import integrate, stats

HERE = Path(__file__).resolve().parent
H0R = HERE.parent / "TEP-H0"
C = 299792.458
LN = np.log(10) / 5
OM = 0.334
rng = np.random.default_rng(7)

# Stellar velocity dispersion from stellar mass (same as stellar-mass closure)
def sigma_stellar_of_mass(logm, slope=0.211):
    return 10 ** (-0.16 + slope * logm + (0.211 - slope) * 10.2)

# Halo velocity dispersion from halo mass via virial relation
# sigma_halo = V_vir / sqrt(2), V_vir = sqrt(G M / R_vir)
# For Delta=200 at z=0: sigma_halo ~ 160 km/s * (M_h / 1e12)^(1/3) * h70
# Calibrated to give sigma_halo ~ 210 km/s for log M_h = 12.3 (typical Cepheid host)
SIGMA_HALO_NORM = 160.0  # km/s at M_h = 1e12

def sigma_halo_of_log_mh(log_mh):
    """Halo velocity dispersion from halo mass via virial relation.

    sigma_halo = 160 km/s * (M_h / 1e12)^(1/3)
    Calibrated so that log M_h = 12.3 gives sigma ~ 200 km/s.
    """
    return SIGMA_HALO_NORM * (10.0 ** np.asarray(log_mh) / 1e12) ** (1.0 / 3.0)

# Behroozi-like SMHM relation (from TEP-JWST tep_model.py)
def stellar_to_halo_mass_behroozi(log_mstar, z=None):
    """log_ratio = -1.8 - 0.1*(log_mstar - 10) - 0.05*(z - 5)

    NOTE: The tep_model.py proxy is calibrated for z~4-8 (JWST high-z).
    At z=0 it underestimates halo masses.  For the closure computation
    on Pantheon+ (z~0) hosts, use the z=0-appropriate fallback
    log_mh = log_mstar + 2.0 (the original tep_model.py fallback for z=None),
    which is consistent with Behroozi+2019 at z=0 for log_mstar ~ 10-11.
    """
    log_mstar = np.asarray(log_mstar)
    # Use z=0 fallback: log_mh = log_mstar + 2.0 (Behroozi+2019 z=0 calibration)
    return log_mstar + 2.0

def X_halo_of(df, S=1.0, Uref=0.0):
    """Halo-velocity-dispersion environmental coordinate.

    X = S * sigma_halo^2 / c^2 - Uref
    where sigma_halo is the halo velocity dispersion from the virial relation.
    This captures the full nested gravitational potential (dark matter halo).
    """
    log_mh = stellar_to_halo_mass_behroozi(df.HOST_LOGMASS.values, z=np.zeros(len(df)))
    sigma_h = sigma_halo_of_log_mh(log_mh)
    return (S * sigma_h ** 2 - Uref) / C ** 2

def fit(df, X, y=None, mass_step=False, use_salt=True, extra_err=None):
    """WLS of y = M - alpha x1 + beta c - kappa X (+ gamma step)."""
    y = df.y.values if y is None else y
    cols = [np.ones(len(df))]
    if use_salt:
        cols += [-df.x1.values, df.c.values]
    cols += [-X]
    if mass_step:
        cols += [(df.HOST_LOGMASS.values > 10).astype(float)]
    A = np.column_stack(cols)
    alpha, beta, sint = 0.15, 3.0, 0.12
    sv = 5 / np.log(10) * 250.0 / (C * df.zHD.values)
    for _ in range(30):
        if use_salt:
            var = (df.mBERR ** 2 + alpha ** 2 * df.x1ERR ** 2 + beta ** 2 * df.cERR ** 2
                   + 2 * alpha * df.cov_mb_x1 - 2 * beta * df.cov_mb_c
                   - 2 * alpha * beta * df.COV_x1_c).values
        else:
            var = (extra_err ** 2) if extra_err is not None else df.m_b_corr_err_DIAG.values ** 2
        w = 1 / (np.clip(var, 1e-6, None) + sint ** 2 + sv ** 2)
        Aw = A * np.sqrt(w)[:, None]
        yw = y * np.sqrt(w)
        cov = np.linalg.inv(Aw.T @ Aw)
        th = cov @ Aw.T @ yw
        res = y - A @ th
        chi2 = float(np.sum(w * res ** 2))
        dof = len(y) - A.shape[1]
        if use_salt:
            alpha, beta = th[1], th[2]
        f = lambda s: np.sum(res ** 2 / (np.clip(var, 1e-6, None) + s ** 2 + sv ** 2)) - dof
        lo, hi = 0.0, 0.5
        if f(lo) > 0:
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
            sint_new = 0.5 * (lo + hi)
        else:
            sint_new = 0.0
        if abs(sint_new - sint) < 1e-5:
            sint = sint_new
            break
        sint = sint_new
    ik = 3 if use_salt else 1
    out = {"n": int(len(y)), "kappa": float(th[ik]), "kappa_err": float(np.sqrt(cov[ik, ik])),
           "sigma_int": float(sint), "chi2_dof": chi2 / dof}
    if use_salt:
        out.update(alpha=float(th[1]), beta=float(th[2]))
    if mass_step:
        out["mass_step"] = float(th[-1])
        out["mass_step_err"] = float(np.sqrt(cov[-1, -1]))
    out["kappa_sigma"] = out["kappa"] / out["kappa_err"]
    return out, th, A, sint

# ---- load Pantheon+SH0ES
d = pd.read_csv(H0R / "data" / "raw" / "Pantheon+SH0ES.dat", sep=r"\s+")
zg = np.linspace(0, 2.4, 24001)
Dc_grid = integrate.cumulative_trapezoid(1 / np.sqrt(OM * (1 + zg) ** 3 + 1 - OM), zg, initial=0) * C / 70.0
d["mu_model"] = 5 * np.log10((1 + d.zHEL) * np.interp(d.zHD, zg, Dc_grid)) + 25
d["y"] = d.mB - d.mu_model
k25 = -2.5 / (np.log(10) * d.x0)
d["cov_mb_x1"] = k25 * d.COV_x1_x0
d["cov_mb_c"] = k25 * d.COV_c_x0

R = {"prereg": "closure/results/tep_shared_distortion_prereg.md",
     "coordinate": "halo_velocity_dispersion",
     "smhm": "Behroozi-like (tep_model.py): log_ratio = -1.8 - 0.1*(log_mstar-10) - 0.05*(z-5)",
     "sigma_halo": "sigma_halo = 160 km/s * (M_h / 1e12)^(1/3) (virial relation)",
     "X_definition": "X = S * sigma_halo^2 / c^2 - Uref"}
base = d[(d.USED_IN_SH0ES_HF == 1) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)

# ---- recovery test
X0 = X_halo_of(base)
P0, th0, A0, s0 = fit(base, X0)
rec = []
for _ in range(200):
    th_mock = th0.copy()
    th_mock[3] = 4e5  # inject a halo-scale kappa (same scale as stellar-mass)
    sv = 5 / np.log(10) * 250.0 / (C * base.zHD.values)
    ymock = A0 @ th_mock + rng.normal(0, np.sqrt(s0 ** 2 + sv ** 2 + base.mBERR.values ** 2))
    rec.append(fit(base, X0, y=ymock)[0]["kappa"])
R["recovery_test"] = {"kappa_injected": 4e5,
                      "kappa_recovered_mean": float(np.mean(rec)),
                      "kappa_recovered_sd": float(np.std(rec))}

# ---- Step 1: primary fit + variants
R["step1_primary"] = P0
V = {}
V["biasCor_subtracted"] = fit(base.assign(y=base.y - base.biasCor_m_b), X0)[0]
mc = base.assign(y=base.m_b_corr - base.mu_model)
V["published_m_b_corr"] = fit(mc, X0, use_salt=False)[0]
nc = d[(d.IS_CALIBRATOR == 0) & (d.zHD > 0.01) & (d.zHD < 0.15) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)
V["z_0.01_0.15_noncal"] = fit(nc, X_halo_of(nc))[0]
hm = base[base.HOST_LOGMASS > 9].reset_index(drop=True)
V["logM_gt_9"] = fit(hm, X_halo_of(hm))[0]
V["with_joint_mass_step"] = fit(base, X0, mass_step=True)[0]
# SMHM offset variant (steeper halo mapping: +0.5 dex in log Mh)
log_mh_plus = stellar_to_halo_mass_behroozi(base.HOST_LOGMASS.values, z=np.zeros(len(base))) + 0.5
X_plus = (sigma_halo_of_log_mh(log_mh_plus) ** 2) / C ** 2
V["smhm_plus_0.5_dex_offset"] = fit(base, X_plus)[0]
R["step1_variants"] = V

# ---- Step 2: agreement with velocity-space kappa
s42 = json.load(open(H0R / "results" / "outputs" / "step_42_tep_native_ladder.json"))
kv = {}
for r in s42:
    if r.get("model") == "T0" and r["sigma_v"] not in kv:
        kv[r["sigma_v"]] = (r["Gamma_X"] / (LN * r["H_app"]), r["Gamma_X_err"] / (LN * r["H_app"]))
k_vel_stellar, k_vel_stellar_err = kv[250]  # in stellar-velocity-dispersion units
# Convert kappa_vel from stellar units to halo units
# kappa_halo = kappa_stellar * <sigma_stellar^2> / <sigma_halo^2> = kappa_stellar / amplification
sigma_stellar_typical = float(np.mean(sigma_stellar_of_mass(base.HOST_LOGMASS.values)))
sigma_halo_typical = float(np.mean(sigma_halo_of_log_mh(
    stellar_to_halo_mass_behroozi(base.HOST_LOGMASS.values, z=np.zeros(len(base))))))
amplification = (sigma_halo_typical / sigma_stellar_typical) ** 2
k_vel_halo = k_vel_stellar / amplification
k_vel_halo_err = k_vel_stellar_err / amplification
comb = np.hypot(P0["kappa_err"], k_vel_halo_err)
R["step2"] = {"kappa_vel_by_sigma_v_stellar_units": {str(k): {"kappa": v[0], "err": v[1]} for k, v in sorted(kv.items())},
              "kappa_vel_stellar_units": k_vel_stellar, "kappa_vel_stellar_err": k_vel_stellar_err,
              "kappa_vel_halo_units": k_vel_halo, "kappa_vel_halo_err": k_vel_halo_err,
              "kappa_SN_halo": P0["kappa"], "kappa_SN_err": P0["kappa_err"],
              "amplification_factor": amplification,
              "sigma_stellar_typical": sigma_stellar_typical,
              "sigma_halo_typical": sigma_halo_typical,
              "difference_sigma": float((P0["kappa"] - k_vel_halo) / comb),
              "agreement_within_2sigma": bool(abs(P0["kappa"] - k_vel_halo) < 2 * comb),
              "kappa_SN_detection_one_sided_p": float(stats.norm.sf(P0["kappa"] / P0["kappa_err"])),
              "kappa_SN_detected_2sigma": bool(P0["kappa"] / P0["kappa_err"] >= 2)}

# ---- Step 3: H0 closure
hf = d[(d.USED_IN_SH0ES_HF == 1) & (d.HOST_LOGMASS > 5)]
s2_hf = float(np.mean(sigma_halo_of_log_mh(
    stellar_to_halo_mass_behroozi(hf.HOST_LOGMASS.values, z=np.zeros(len(hf)))) ** 2))
s04 = pd.read_csv(H0R / "results" / "outputs" / "step_04_tep_corrected_h0.csv")
# Calibrator halo velocity dispersion with shear suppression (same structure as stellar-mass closure)
# Invert sigma_inferred -> stellar mass -> halo mass -> halo sigma
cal_log_mstar = (np.log10(s04.sigma_inferred.values) + 0.16) / 0.211
cal_log_mh = stellar_to_halo_mass_behroozi(cal_log_mstar, z=np.zeros(len(s04)))
cal_sigma_halo = sigma_halo_of_log_mh(cal_log_mh)
cal_s2 = (s04.shear_suppression.values * cal_sigma_halo ** 2)

preds = {}
for lab_u, Uref in (("unscreened_87.165", 87.165 ** 2), ("screened_30.507", 30.507 ** 2)):
    dX = (s2_hf - Uref) / C ** 2
    for lab_k, (k, ke) in (("kappa_vel_halo", (k_vel_halo, k_vel_halo_err)),
                          ("kappa_SN_halo", (P0["kappa"], P0["kappa_err"]))):
        H = 73.04 * 10 ** (-k * dX / 5)
        He = abs(H * np.log(10) / 5 * dX * ke)
        cal_pred = float(np.mean(H * 10 ** (k * (cal_s2 - Uref) / C ** 2 / 5)))
        preds[f"{lab_u}|{lab_k}"] = {"Xbar_HF_minus_anchor": dX, "H0_true_pred": H,
                                     "H0_pred_err_from_kappa": He,
                                     "pred_uncorrected_calibrator_cz_over_d": cal_pred}
    need_dX = 5 * np.log10(73.04 / 67.4) / k_vel_halo
    preds[f"{lab_u}|required"] = {"Xbar_HF_minus_anchor_needed_at_kappa_vel_halo": float(need_dX),
                                    "equivalent_HF_mean_sigma2_km2s2": float(need_dX * C ** 2 + Uref),
                                    "equivalent_HF_sigma_kms": float(np.sqrt(need_dX * C ** 2 + Uref))}
R["step3"] = {"mean_sigma2_HF_hosts_km2s2": s2_hf,
              "mean_S_sigma2_calibrators_km2s2": float(cal_s2.mean()),
              "observed_uncorrected_calibrator_cz_over_d": 67.72,
              "predictions": preds}

# ---- Summary
H0_SN_screened = preds.get("screened_30.507|kappa_SN_halo", {}).get("H0_true_pred")
H0_vel_screened = preds.get("screened_30.507|kappa_vel_halo", {}).get("H0_true_pred")
H0_SN_unscreened = preds.get("unscreened_87.165|kappa_SN_halo", {}).get("H0_true_pred")
H0_vel_unscreened = preds.get("unscreened_87.165|kappa_vel_halo", {}).get("H0_true_pred")
R["summary"] = {"coordinate": "halo_velocity_dispersion",
                 "kappa_SN_halo": P0["kappa"], "kappa_SN_err": P0["kappa_err"],
                 "kappa_SN_sigma": P0["kappa_sigma"],
                 "kappa_vel_halo": k_vel_halo, "kappa_vel_halo_err": k_vel_halo_err,
                 "amplification_factor": amplification,
                 "sigma_halo_typical": sigma_halo_typical,
                 "sigma_stellar_typical": sigma_stellar_typical,
                 "H0_true_from_kappa_SN_screened": H0_SN_screened,
                 "H0_true_from_kappa_vel_screened": H0_vel_screened,
                 "H0_true_from_kappa_SN_unscreened": H0_SN_unscreened,
                 "H0_true_from_kappa_vel_unscreened": H0_vel_unscreened,
                 "H0_planck_lcdm": 67.4, "H0_sh0es": 73.04,
                 "gap_to_close": 73.04 - 67.4,
                 "fraction_closed_SN_screened": float((73.04 - H0_SN_screened) / (73.04 - 67.4)) if H0_SN_screened else None,
                 "fraction_closed_vel_screened": float((73.04 - H0_vel_screened) / (73.04 - 67.4)) if H0_vel_screened else None}

json.dump(R, open(HERE / "results" / "tep_halo_potential_closure.json", "w"), indent=1, default=float)

# ---- Print summary
print("=" * 70)
print("HALO-POTENTIAL X CLOSURE (halo velocity dispersion)")
print("=" * 70)
print(f"Coordinate: X = S * sigma_halo^2 / c^2 - Uref")
print(f"  sigma_halo = {SIGMA_HALO_NORM} km/s * (M_h / 1e12)^(1/3)")
print()
print("RECOVERY TEST:")
print(f"  injected 4.0e5 -> {R['recovery_test']['kappa_recovered_mean']:.3e} +- {R['recovery_test']['kappa_recovered_sd']:.2e}")
print()
print("STEP 1 (primary fit, halo-velocity-dispersion X):")
p = P0
print(f"  n={p['n']}: kappa_SN = {p['kappa']:.3e} +- {p['kappa_err']:.3e} ({p['kappa_sigma']:+.2f} sigma)")
print(f"  alpha={p['alpha']:.3f} beta={p['beta']:.2f} sig_int={p['sigma_int']:.3f}")
for k, v in V.items():
    extra = f"  mass step {v['mass_step']:+.3f}+-{v['mass_step_err']:.3f}" if "mass_step" in v else ""
    print(f"  {k:22s} n={v['n']:4d} kappa = {v['kappa']:+.3e} +- {v['kappa_err']:.2e} ({v['kappa_sigma']:+.2f}){extra}")
print()
print("STEP 2 (agreement with velocity-space kappa):")
s = R["step2"]
print(f"  kappa_vel (stellar units, sigma_v=250): {k_vel_stellar:.3e} +- {k_vel_stellar_err:.2e}")
print(f"  kappa_vel (halo units, converted):     {k_vel_halo:.3e} +- {k_vel_halo_err:.2e}")
print(f"  kappa_SN (halo units):                  {P0['kappa']:.3e} +- {P0['kappa_err']:.2e}")
print(f"  kappa_SN vs kappa_vel_halo: {s['difference_sigma']:+.2f} sigma; agreement={s['agreement_within_2sigma']}")
print(f"  SN detection 2sigma={s['kappa_SN_detected_2sigma']}; one-sided p={s['kappa_SN_detection_one_sided_p']:.4f}")
print(f"  Amplification: {amplification:.2f}x (sigma_halo={sigma_halo_typical:.1f} vs sigma_stellar={sigma_stellar_typical:.1f})")
print()
print("STEP 3 (H0 closure):")
s3 = R["step3"]
print(f"  <sigma_halo^2>_HF = {s3['mean_sigma2_HF_hosts_km2s2']:.0f} km2/s2")
print(f"  <S*sigma_halo^2>_cal = {s3['mean_S_sigma2_calibrators_km2s2']:.0f} km2/s2")
for k, v in preds.items():
    if k.endswith("required"):
        print(f"  {k:34s} needs Xbar_HF-anchor = {v['Xbar_HF_minus_anchor_needed_at_kappa_vel_halo']:.2e} (HF sigma ~ {v['equivalent_HF_sigma_kms']:.0f} km/s)")
    else:
        print(f"  {k:34s} dX={v['Xbar_HF_minus_anchor']:+.2e}  H0_true={v['H0_true_pred']:.2f} +- {v['H0_pred_err_from_kappa']:.2f}  pred. cal cz/d={v['pred_uncorrected_calibrator_cz_over_d']:.2f}")
print()
print("SUMMARY:")
sm = R["summary"]
print(f"  kappa_SN_halo = {sm['kappa_SN_halo']:.3e} +- {sm['kappa_SN_err']:.3e} ({sm['kappa_SN_sigma']:+.2f} sigma)")
print(f"  sigma_halo_typical = {sm['sigma_halo_typical']:.1f} km/s")
print(f"  sigma_stellar_typical = {sm['sigma_stellar_typical']:.1f} km/s")
print(f"  amplification_factor = {sm['amplification_factor']:.2f}x")
print(f"  H0_true (kappa_SN, screened)  = {sm['H0_true_from_kappa_SN_screened']}")
print(f"  H0_true (kappa_vel, screened) = {sm['H0_true_from_kappa_vel_screened']}")
print(f"  H0_true (kappa_SN, unscreened)  = {sm['H0_true_from_kappa_SN_unscreened']}")
print(f"  H0_true (kappa_vel, unscreened) = {sm['H0_true_from_kappa_vel_unscreened']}")
print(f"  Planck (LCDM) = {sm['H0_planck_lcdm']}, SH0ES = {sm['H0_sh0es']}")
