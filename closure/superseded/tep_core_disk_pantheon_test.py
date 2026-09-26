"""TEP core-disk observing-chain test on Pantheon+ (pre-registered: results/tep_core_disk_prereg.md).

Prediction: Hubble residual r_i = -k T_i with k > 0, where
  D_i = (sigma_i^2/c^2) ln(1 + R_i^2/R_s^2)   core-to-SN potential difference
  T_i = D_i (dmu/dz + 1.5)                   distance (1/z) + time-dilation response
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import integrate, stats

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "TEP-H0" / "data" / "raw"
C_KMS = 299792.458
OM = 0.334
rng = np.random.default_rng(2026)

d = pd.read_csv(RAW / "Pantheon+SH0ES.dat", sep=r"\s+")
N = len(d)
with open(RAW / "external" / "Pantheon+SH0ES_STAT+SYS.cov") as f:
    assert int(f.readline()) == N
COV = np.loadtxt(RAW / "external" / "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(N, N)

# --- cosmology (H0 = 70 placeholder; absorbed by intercept)
zg = np.linspace(0, 2.4, 24001)
Ez = lambda z: np.sqrt(OM * (1 + z) ** 3 + 1 - OM)
Dc_grid = integrate.cumulative_trapezoid(1 / Ez(zg), zg, initial=0) * C_KMS / 70.0  # Mpc
Dc = lambda z: np.interp(z, zg, Dc_grid)
d["Dc"] = Dc(d.zHD)
d["dL"] = (1 + d.zHEL) * d.Dc
d["mu_model"] = 5 * np.log10(d.dL) + 25
dz = 1e-4
d["dmu_dz"] = 5 / np.log(10) * (np.log(Dc(d.zHD + dz)) - np.log(Dc(d.zHD - dz).clip(1e-9))) / (2 * dz)
d["resid"] = d.m_b_corr - d.mu_model
d["dA_kpc_per_arcsec"] = d.Dc / (1 + d.zHD) * 1e3 / 206264.806

def build(df, slope=0.211, Rs_mode="scaled", Rs_fixed=3.0):
    sig = 10 ** (-0.16 + slope * df.HOST_LOGMASS + (0.211 - slope) * 10.2)  # pivot at logM=10.2 keeps median sigma
    R = df.HOST_ANGSEP * df.dA_kpc_per_arcsec
    Rs = 3.0 * 10 ** (0.25 * (df.HOST_LOGMASS - 10.5)) if Rs_mode == "scaled" else Rs_fixed
    D = (sig / C_KMS) ** 2 * np.log1p((R / Rs) ** 2)
    return D.values, np.log10((R / Rs).values), (D * (df.dmu_dz + 1.5)).values

def gls(idx, cols, diag=False):
    y = d.resid.values[idx]
    Cs = np.diag(d.m_b_corr_err_DIAG.values[idx] ** 2) if diag else COV[np.ix_(idx, idx)]
    L = np.linalg.cholesky(Cs)
    A = np.column_stack([np.ones(len(idx))] + cols)
    Aw, yw = np.linalg.solve(L, A), np.linalg.solve(L, y)
    cov = np.linalg.inv(Aw.T @ Aw)
    th = cov @ Aw.T @ yw
    chi2 = float(np.sum((yw - Aw @ th) ** 2))
    return th, cov, chi2, L, Aw, yw

def select(zmin=0.01, zmax=0.8, dedup=False, survey_out=None):
    m = (d.IS_CALIBRATOR == 0) & (d.HOST_ANGSEP > 0) & (d.HOST_LOGMASS > 5) & (d.zHD > zmin) & (d.zHD < zmax)
    if survey_out is not None:
        m &= d.IDSURVEY != survey_out
    idx = np.where(m)[0]
    if dedup:
        keep = ~d.iloc[idx].CID.duplicated().values
        idx = idx[keep]
    return idx

def fit(idx, controls=True, diag=False, **kw):
    df = d.iloc[idx]
    D, logRn, T = build(df, **kw)
    cols = ([(df.HOST_LOGMASS.values > 10).astype(float), logRn] if controls else []) + [T]
    th, cov, chi2, L, Aw, yw = gls(idx, cols, diag)
    k, ke = -th[-1], np.sqrt(cov[-1, -1])
    return {"n": int(len(idx)), "k": float(k), "k_err": float(ke), "k_sigma": float(k / ke),
            "p_one_sided_k_gt_0": float(stats.norm.sf(k / ke)), "chi2": chi2,
            "mass_step": float(th[1]) if controls else None, "offset_slope": float(th[2]) if controls else None,
            "_D": D, "_T": T, "_L": L, "_Aw": Aw, "_yw": yw, "_df": df}

def clean(r):
    return {k: v for k, v in r.items() if not k.startswith("_")}

out = {"prereg": "closure/results/tep_core_disk_prereg.md"}
P = fit(select())
out["primary"] = clean(P)

# permutation of D across SNe (keeps z, mass, offset controls fixed)
df = P["_df"]; D = P["_D"]; w = (df.dmu_dz + 1.5).values
Aw = P["_Aw"].copy(); L = P["_L"]; yw = P["_yw"]
ks = []
for _ in range(2000):
    Tp = rng.permutation(D) * w
    Aw[:, -1] = np.linalg.solve(L, Tp)
    th = np.linalg.lstsq(Aw, yw, rcond=None)[0]
    ks.append(-th[-1])
ks = np.array(ks)
out["permutation"] = {"n_draws": 2000, "p_k_ge_observed": float(np.mean(ks >= P["k"])), "null_k_sd": float(ks.std())}

V = {}
V["z_0.01_0.15"] = fit(select(0.01, 0.15))
V["z_0.023_0.15"] = fit(select(0.023, 0.15))
V["z_all"] = fit(select(0.01, 3.0))
for rs in (1.0, 3.0, 5.0):
    V[f"Rs_fixed_{rs:g}kpc"] = fit(select(), Rs_mode="fixed", Rs_fixed=rs)
V["sigma_slope_x0.5"] = fit(select(), slope=0.211 * 0.5)
V["sigma_slope_x1.5"] = fit(select(), slope=0.211 * 1.5)
V["diag_errors"] = fit(select(), diag=True)
V["no_controls"] = fit(select(), controls=False)
V["dedup"] = fit(select(dedup=True))
for s in sorted(d.iloc[select()].IDSURVEY.value_counts().index[:8]):
    V[f"leave_out_survey_{s}"] = fit(select(survey_out=s))
out["robustness"] = {k: clean(v) for k, v in V.items()}
signs = [v["k"] > 0 for v in V.values()]
out["robustness_same_sign_as_TEP_fraction"] = float(np.mean(signs))

# amplitude required to resolve the tension
hf = d.iloc[select()]
hfm = hf.USED_IN_SH0ES_HF.values == 1
_, _, T_all = build(hf)
cal = d[(d.IS_CALIBRATOR == 1) & (d.HOST_ANGSEP > 0) & (d.HOST_LOGMASS > 5)]
D_cal = build(cal)[0] if len(cal) else np.array([])
target = np.log(73.04 / 67.4)
mean_T_hf = float(T_all[hfm].mean()) if hfm.any() else float("nan")
cal_term = 1.5 * float(D_cal.mean()) if len(D_cal) else 0.0
k_req = target / ((np.log(10) / 5) * (mean_T_hf - cal_term))
k, ke = P["k"], P["k_err"]
out["amplitude"] = {
    "n_SH0ES_HF_in_sample": int(hfm.sum()), "mean_T_HF": mean_T_hf,
    "n_calibrators_with_offsets": int(len(cal)), "calibrator_time_dilation_term": cal_term,
    "k_required_for_H0_67p4": float(k_req),
    "k_fitted": k, "k_fitted_95": [k - 1.96 * ke, k + 1.96 * ke],
    "k_required_in_95_interval": bool(k - 1.96 * ke <= k_req <= k + 1.96 * ke),
    "k_required_excluded_sigma": float((k_req - k) / ke),
    "H0_shift_implied_by_fit_kms_Mpc": float(73.04 * (np.exp(-(np.log(10) / 5) * k * (mean_T_hf - cal_term)) - 1)),
    "median_D": float(np.median(D)), "median_redshift_offset_at_k_req_kms": float(k_req * np.median(D) * C_KMS),
}
pk = out["primary"]
out["verdict"] = (
    "SUPPORT" if (pk["k"] > 0 and pk["p_one_sided_k_gt_0"] < 0.01 and out["permutation"]["p_k_ge_observed"] < 0.01
                  and out["robustness_same_sign_as_TEP_fraction"] >= 0.8)
    else "EXCLUDED_AT_PRIMARY" if pk["k"] <= 0 else "NOT_SIGNIFICANT")
(HERE / "results").mkdir(exist_ok=True)
json.dump(out, open(HERE / "results" / "tep_core_disk_pantheon_test.json", "w"), indent=1)

print(f"PRIMARY n={pk['n']}  k = {pk['k']:.3e} +- {pk['k_err']:.3e}  ({pk['k_sigma']:+.2f} sigma, one-sided p={pk['p_one_sided_k_gt_0']:.3g})")
print(f"  mass step {pk['mass_step']:+.3f} mag, offset slope {pk['offset_slope']:+.3f} mag/dex")
print(f"PERMUTATION p(k>=obs) = {out['permutation']['p_k_ge_observed']:.3f}  null sd {out['permutation']['null_k_sd']:.2e}")
for n, v in out["robustness"].items():
    print(f"  {n:24s} n={v['n']:4d} k={v['k']:+.2e} +- {v['k_err']:.2e} ({v['k_sigma']:+.2f})")
print(f"same-sign fraction: {out['robustness_same_sign_as_TEP_fraction']:.2f}")
a = out["amplitude"]
print(f"AMPLITUDE: k_req = {a['k_required_for_H0_67p4']:.3e} (SH0ES HF rows in sample {a['n_SH0ES_HF_in_sample']}); "
      f"fitted 95% [{a['k_fitted_95'][0]:.2e}, {a['k_fitted_95'][1]:.2e}]; k_req excluded at {a['k_required_excluded_sigma']:.1f} sigma; "
      f"implied dH0 = {a['H0_shift_implied_by_fit_kms_Mpc']:+.2f}; redshift offset at k_req ~ {a['median_redshift_offset_at_k_req_kms']:.0f} km/s")
print("VERDICT:", out["verdict"])
