"""TEP (achromatic clock distortion) vs dust (chromatic) discriminator for the Step 1 kappa_SN signal.
Pre-registered: results/tep_shared_distortion_prereg.md (addendum).
"""
import io, json, contextlib, importlib.util
from pathlib import Path
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("sd", HERE / "tep_shared_distortion_closure.py")
m = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(m)

d = m.d
base = d[(d.USED_IN_SH0ES_HF == 1) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)
wide = d[(d.IS_CALIBRATOR == 0) & (d.zHD > 0.01) & (d.zHD < 0.15) & (d.HOST_LOGMASS > 5)].reset_index(drop=True)

def split(df, col, cut):
    lo = df[df[col] < cut].reset_index(drop=True); hi = df[df[col] >= cut].reset_index(drop=True)
    a = m.fit(lo, m.X_of(lo))[0]; b = m.fit(hi, m.X_of(hi))[0]
    diff = (b["kappa"] - a["kappa"]) / np.hypot(a["kappa_err"], b["kappa_err"])
    return {"low": a, "high": b, "high_minus_low_sigma": float(diff)}

def interaction(df):
    X = m.X_of(df); cc = df.c.values - df.c.values.mean()
    y = df.y.values
    A = np.column_stack([np.ones(len(df)), -df.x1.values, df.c.values, -X, -X * cc])
    sv = 5 / np.log(10) * 250.0 / (m.C * df.zHD.values)
    base_fit = m.fit(df, X)[0]
    al, be, si = base_fit["alpha"], base_fit["beta"], base_fit["sigma_int"]
    var = (df.mBERR ** 2 + al ** 2 * df.x1ERR ** 2 + be ** 2 * df.cERR ** 2 + 2 * al * df.cov_mb_x1
           - 2 * be * df.cov_mb_c - 2 * al * be * df.COV_x1_c).values
    w = 1 / (np.clip(var, 1e-6, None) + si ** 2 + sv ** 2)
    Aw = A * np.sqrt(w)[:, None]; cov = np.linalg.inv(Aw.T @ Aw); th = cov @ Aw.T @ (y * np.sqrt(w))
    return {"kappa": float(th[3]), "kappa_err": float(np.sqrt(cov[3, 3])),
            "kappa_c": float(th[4]), "kappa_c_err": float(np.sqrt(cov[4, 4])),
            "kappa_c_sigma": float(th[4] / np.sqrt(cov[4, 4]))}

out = {}
for lab, df in (("SH0ES_HF", base), ("wide_0.01_0.15", wide)):
    col = split(df, "c", 0.0); x1s = split(df, "x1", 0.0); inter = interaction(df)
    blue, red = col["low"], col["high"]
    tep = (abs(col["high_minus_low_sigma"]) < 2 and blue["kappa"] / blue["kappa_err"] >= 1.5 and abs(inter["kappa_c_sigma"]) < 2)
    dust = ((abs(blue["kappa"] / blue["kappa_err"]) < 1 and col["high_minus_low_sigma"] > 2) or inter["kappa_c_sigma"] >= 2)
    out[lab] = {"colour_split": col, "stretch_split": x1s, "colour_interaction": inter,
                "verdict": "ACHROMATIC_SUPPORTS_TEP" if tep else ("CHROMATIC_SUPPORTS_DUST" if dust else "UNDETERMINED")}
json.dump(out, open(HERE / "results" / "tep_shared_distortion_discriminator.json", "w"), indent=1)
for lab, o in out.items():
    c, s, i = o["colour_split"], o["stretch_split"], o["colour_interaction"]
    print(f"== {lab}")
    print(f"  blue c<0  n={c['low']['n']:3d} kappa={c['low']['kappa']:+.2e}+-{c['low']['kappa_err']:.2e} ({c['low']['kappa_sigma']:+.2f})")
    print(f"  red  c>=0 n={c['high']['n']:3d} kappa={c['high']['kappa']:+.2e}+-{c['high']['kappa_err']:.2e} ({c['high']['kappa_sigma']:+.2f})   red-blue {c['high_minus_low_sigma']:+.2f} sigma")
    print(f"  x1<0      n={s['low']['n']:3d} kappa={s['low']['kappa']:+.2e}+-{s['low']['kappa_err']:.2e} ({s['low']['kappa_sigma']:+.2f})")
    print(f"  x1>=0     n={s['high']['n']:3d} kappa={s['high']['kappa']:+.2e}+-{s['high']['kappa_err']:.2e} ({s['high']['kappa_sigma']:+.2f})   diff {s['high_minus_low_sigma']:+.2f} sigma")
    print(f"  interaction kappa_c = {i['kappa_c']:+.2e} +- {i['kappa_c_err']:.2e} ({i['kappa_c_sigma']:+.2f} sigma)")
    print(f"  VERDICT: {o['verdict']}")
