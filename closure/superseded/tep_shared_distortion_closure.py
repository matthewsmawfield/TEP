"""CORPUS_PLAN §7.8 Target 3 closure computation: shared environmental clock distortion.

Pre-registered in results/tep_shared_distortion_prereg.md.
Step 1  kappa_SN from Pantheon+ Hubble-flow SN magnitudes vs host potential X (no host-mass step)
Step 2  agreement with Paper 11 velocity-space kappa = Gamma_X / ((ln10/5) H_app)
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

d = pd.read_csv(H0R / "data" / "raw" / "Pantheon+SH0ES.dat", sep=r"\s+")
zg = np.linspace(0, 2.4, 24001)
Dc_grid = integrate.cumulative_trapezoid(1 / np.sqrt(OM * (1 + zg) ** 3 + 1 - OM), zg, initial=0) * C / 70.0
d["mu_model"] = 5 * np.log10((1 + d.zHEL) * np.interp(d.zHD, zg, Dc_grid)) + 25
d["y"] = d.mB - d.mu_model
k25 = -2.5 / (np.log(10) * d.x0)
d["cov_mb_x1"] = k25 * d.COV_x1_x0
d["cov_mb_c"] = k25 * d.COV_c_x0

def sigma_of_mass(logm, slope=0.211):
    return 10 ** (-0.16 + slope * logm + (0.211 - slope) * 10.2)

def X_of(df, S=1.0, slope=0.211, Uref=0.0):
    return (S * sigma_of_mass(df.HOST_LOGMASS.values, slope) ** 2 - Uref) / C ** 2

def fit(df, X, y=None, mass_step=False, use_salt=True, extra_err=None):
    """WLS of y = M - alpha x1 + beta c - kappa X (+ gamma step); iterate alpha,beta-dependent errors and sigma_int."""
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
                   + 2 * alpha * df.cov_mb_x1 - 2 * beta * df.cov_mb_c - 2 * alpha * beta * df.COV_x1_c).values
        else:
            var = (extra_err ** 2) if extra_err is not None else df.m_b_corr_err_DIAG.values ** 2
        w = 1 / (np.clip(var, 1e-6, None) + sint ** 2 + sv ** 2)
        Aw = A * np.sqrt(w)[:, None]; yw = y * np.sqrt(w)
        cov = np.linalg.inv(Aw.T @ Aw); th = cov @ Aw.T @ yw
        res = y - A @ th
        chi2 = float(np.sum(w * res ** 2)); dof = len(y) - A.shape[1]
        if use_salt:
            alpha, beta = th[1], th[2]
        # tune sigma_int so chi2/dof = 1
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
            sint = sint_new; break
        sint = sint_new
    ik = 3 if use_salt else 1
    out = {"n": int(len(y)), "kappa": float(th[ik]), "kappa_err": float(np.sqrt(cov[ik, ik])),
           "sigma_int": float(sint), "chi2_dof": chi2 / dof}
    if use_salt:
        out.update(alpha=float(th[1]), beta=float(th[2]))
    if mass_step:
        out["mass_step"] = float(th[-1]); out["mass_step_err"] = float(np.sqrt(cov[-1, -1]))
    out["kappa_sigma"] = out["kappa"] / out["kappa_err"]
    return out, th, A, sint

R = {"prereg": "closure/results/tep_shared_distortion_prereg.md"}
base = d[(d.USED_IN_SH0ES_HF == 1) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)

# ---- recovery test (guards sign and coding errors)
X0 = X_of(base)
P0, th0, A0, s0 = fit(base, X0)
rec = []
for _ in range(200):
    th_mock = th0.copy(); th_mock[3] = 4e5
    sv = 5 / np.log(10) * 250.0 / (C * base.zHD.values)
    ymock = A0 @ th_mock + rng.normal(0, np.sqrt(s0 ** 2 + sv ** 2 + base.mBERR.values ** 2))
    rec.append(fit(base, X0, y=ymock)[0]["kappa"])
R["recovery_test"] = {"kappa_injected": 4e5, "kappa_recovered_mean": float(np.mean(rec)), "kappa_recovered_sd": float(np.std(rec))}

# ---- Step 1
R["step1_primary"] = P0
V = {}
V["biasCor_subtracted"] = fit(base.assign(y=base.y - base.biasCor_m_b), X0)[0]
mc = base.assign(y=base.m_b_corr - base.mu_model)
V["published_m_b_corr"] = fit(mc, X0, use_salt=False)[0]
nc = d[(d.IS_CALIBRATOR == 0) & (d.zHD > 0.01) & (d.zHD < 0.15) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)
V["z_0.01_0.15_noncal"] = fit(nc, X_of(nc))[0]
hm = base[base.HOST_LOGMASS > 9].reset_index(drop=True)
V["logM_gt_9"] = fit(hm, X_of(hm))[0]
V["sigma_slope_x0.5"] = fit(base, X_of(base, slope=0.211 * 0.5))[0]
V["sigma_slope_x1.5"] = fit(base, X_of(base, slope=0.211 * 1.5))[0]
V["S_0.867"] = fit(base, X_of(base, S=0.867))[0]
V["with_joint_mass_step"] = fit(base, X0, mass_step=True)[0]
R["step1_variants"] = V

# ---- Step 2
s42 = json.load(open(H0R / "results" / "outputs" / "step_42_tep_native_ladder.json"))
kv = {}
for r in s42:
    if r.get("model") == "T0" and r["sigma_v"] not in kv:
        kv[r["sigma_v"]] = (r["Gamma_X"] / (LN * r["H_app"]), r["Gamma_X_err"] / (LN * r["H_app"]))
k_vel, k_vel_err = kv[250]
comb = np.hypot(P0["kappa_err"], k_vel_err)
R["step2"] = {"kappa_vel_by_sigma_v": {str(k): {"kappa": v[0], "err": v[1]} for k, v in sorted(kv.items())},
              "kappa_SN": P0["kappa"], "kappa_SN_err": P0["kappa_err"],
              "difference_sigma_vs_canonical": float((P0["kappa"] - k_vel) / comb),
              "agreement_within_2sigma": bool(abs(P0["kappa"] - k_vel) < 2 * comb),
              "kappa_SN_detection_one_sided_p": float(stats.norm.sf(P0["kappa"] / P0["kappa_err"])),
              "kappa_SN_detected_2sigma": bool(P0["kappa"] / P0["kappa_err"] >= 2)}

# ---- Step 3
hf = d[(d.USED_IN_SH0ES_HF == 1) & (d.HOST_LOGMASS > 5)]
s2_hf = float(np.mean(sigma_of_mass(hf.HOST_LOGMASS.values) ** 2))
s04 = pd.read_csv(H0R / "results" / "outputs" / "step_04_tep_corrected_h0.csv")
cal_s2 = (s04.shear_suppression * s04.sigma_inferred ** 2).values
preds = {}
for lab_u, Uref in (("unscreened_87.165", 87.165 ** 2), ("screened_30.507", 30.507 ** 2)):
    dX = (s2_hf - Uref) / C ** 2
    for lab_k, (k, ke) in (("kappa_vel_250", (k_vel, k_vel_err)), ("kappa_SN", (P0["kappa"], P0["kappa_err"])),
                           ("kappa_0.96e6", (0.96e6, 0.0))):
        H = 73.04 * 10 ** (-k * dX / 5)
        He = abs(H * np.log(10) / 5 * dX * ke)
        cal_pred = float(np.mean(H * 10 ** (k * (cal_s2 - Uref) / C ** 2 / 5)))
        preds[f"{lab_u}|{lab_k}"] = {"Xbar_HF_minus_anchor": dX, "H0_true_pred": H, "H0_pred_err_from_kappa": He,
                                     "pred_uncorrected_calibrator_cz_over_d": cal_pred}
    need_dX = 5 * np.log10(73.04 / 67.4) / k_vel
    preds[f"{lab_u}|required"] = {"Xbar_HF_minus_anchor_needed_at_kappa_vel": float(need_dX),
                                  "equivalent_HF_mean_sigma2_km2s2": float(need_dX * C ** 2 + Uref),
                                  "equivalent_HF_sigma_kms": float(np.sqrt(need_dX * C ** 2 + Uref))}
R["step3"] = {"mean_sigma2_HF_hosts_km2s2": s2_hf, "mean_S_sigma2_calibrators_km2s2": float(cal_s2.mean()),
              "observed_uncorrected_calibrator_cz_over_d": 67.72, "predictions": preds}
json.dump(R, open(HERE / "results" / "tep_shared_distortion_closure.json", "w"), indent=1, default=float)

print(f"RECOVERY: injected 4.0e5 -> {R['recovery_test']['kappa_recovered_mean']:.3e} +- {R['recovery_test']['kappa_recovered_sd']:.2e}")
p = P0
print(f"STEP1 primary n={p['n']}: kappa_SN = {p['kappa']:.3e} +- {p['kappa_err']:.3e} ({p['kappa_sigma']:+.2f} sigma); alpha={p['alpha']:.3f} beta={p['beta']:.2f} sig_int={p['sigma_int']:.3f}")
for k, v in V.items():
    extra = f"  mass step {v['mass_step']:+.3f}+-{v['mass_step_err']:.3f}" if "mass_step" in v else ""
    print(f"  {k:22s} n={v['n']:4d} kappa = {v['kappa']:+.3e} +- {v['kappa_err']:.2e} ({v['kappa_sigma']:+.2f}){extra}")
s = R["step2"]
print("STEP2 kappa_vel:", {k: f"{v['kappa']:.2e}+-{v['err']:.2e}" for k, v in s["kappa_vel_by_sigma_v"].items()})
print(f"  kappa_SN vs kappa_vel(250): {s['difference_sigma_vs_canonical']:+.2f} sigma; agreement={s['agreement_within_2sigma']}; SN detection 2sigma={s['kappa_SN_detected_2sigma']}")
print(f"STEP3 <sigma^2>_HF = {s2_hf:.0f}; <S sigma^2>_cal = {cal_s2.mean():.0f} km2/s2")
for k, v in preds.items():
    if k.endswith("required"):
        print(f"  {k:34s} needs Xbar_HF-anchor = {v['Xbar_HF_minus_anchor_needed_at_kappa_vel']:.2e} (HF sigma ~ {v['equivalent_HF_sigma_kms']:.0f} km/s)")
    else:
        print(f"  {k:34s} dX={v['Xbar_HF_minus_anchor']:+.2e}  H0_true={v['H0_true_pred']:.2f} +- {v['H0_pred_err_from_kappa']:.2f}  pred. calibrator cz/d={v['pred_uncorrected_calibrator_cz_over_d']:.2f}")
