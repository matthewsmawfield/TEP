"""Cepheid environmental term through Paper 11's own full SH0ES matrix machinery.

Uses TEP-H0/scripts/steps/step_45 fit_sh0es_with_tep_correction (fixed kappa, all
other ladder parameters refit) to profile chi^2(kappa) under every reference origin,
three anchor-coordinate conventions and diagonal vs full covariance, and records the
provenance of the H0 = 66.65 and 25.24 sigma figures. Output:
results/tep_cepheid_matrix_profile.json
"""
import sys, io, json, math, contextlib, importlib.util
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
H0R = HERE.parent / "TEP-H0"
sys.path.insert(0, str(H0R)); sys.path.insert(0, str(H0R / "scripts" / "steps"))
spec = importlib.util.spec_from_file_location("s45", H0R / "scripts/steps/step_45_full_ladder_h0_propagation.py")
m = importlib.util.module_from_spec(spec)
quiet = lambda: contextlib.redirect_stdout(io.StringIO())
with quiet(): spec.loader.exec_module(m)
with quiet():
    L, y, C, q, src = m.load_sh0es_data(include_sources=True)
    hs, hz, hS = m.load_host_metadata()
    s_std, s_scr = m.compute_sigma_ref(screened=False), m.compute_sigma_ref(screened=True)

def fit(k, sr, hs_=hs, CC=C):
    with quiet():
        return m.fit_sh0es_with_tep_correction(L, y, CC, q, src, hs_, hS, k, sr)

def profile(hs_=hs, CC=C, sr=s_std):
    ks = np.array([-0.4, -0.2, 0.0, 0.2, 0.4]) * 1e6
    ch = [fit(k, sr, hs_, CC)["chi2"] for k in ks]
    a, b, _ = np.polyfit(ks / 1e6, ch, 2)
    kb, s = -b / (2 * a), 1 / math.sqrt(a)
    return {"kappa_best_1e6": kb, "kappa_err_1e6": s, "canonical_0p96_excluded_sigma": (0.96 - kb) / s,
            "endpoint_0p365_offset_sigma": (0.365 - kb) / s, "kappa_95_upper_1e6": kb + 1.96 * s}

R = {}
R["grid_standard_origin"] = [{"kappa": k, "chi2": (r := fit(k, s_std))["chi2"], "H0": r["H0"], "H0_err": r["H0_err"],
                              "n_rows": r["n_corrected_rows"]} for k in np.array([-0.6, -0.2, 0, 0.365, 0.96, 1.4]) * 1e6]
R["origin_invariance"] = {lab: fit(0.96e6, sr)["chi2"] for lab, sr in
                          (("87.165", s_std), ("30.507", s_scr), ("0", 1e-6), ("150", 150.0))}
R["profile_published_anchors"] = profile()
h2 = dict(hs); h2.update({"MW": 230 / math.sqrt(2), "LMC": 69 / math.sqrt(2), "N4258": 208 / math.sqrt(2)})
R["profile_anchors_Vrot_over_sqrt2"] = profile(h2)
h3 = {k: v for k, v in hs.items() if k not in ("MW", "LMC", "N4258", "SMC", "M31")}
R["profile_anchors_uncorrected"] = profile(h3)
R["profile_diagonal_covariance"] = profile(CC=np.diag(np.diag(C)))
out = H0R / "results" / "outputs"
s04 = json.load(open(out / "step_04_tep_correction_results.json"))
s34 = json.load(open(out / "step_34_full_ladder_likelihood_results.json"))
R["provenance"] = {
  "H0_66p65": {"source": "step_04 unified_h0 = unweighted mean of per-host cz/d after in-sample kappa that zeroes an unweighted slope",
               "unified_h0": s04["unified_h0"], "unified_h0_screened_origin": s04["unified_h0_screened"],
               "origin_dependence_kms_Mpc": s04["unified_h0"] - s04["unified_h0_screened"],
               "optimal_kappa": s04["optimal_kappa_cep"], "wls_kappa": s04["wls_kappa"], "wls_kappa_err_scaled": s04["wls_kappa_err_scaled"],
               "bootstrap_kappa_mean": s04["bootstrap_kappa_mean"], "bootstrap_kappa_std": s04["bootstrap_kappa_std"]},
  "kappa_25p24_sigma": {"source": "step_34 stage1 summary likelihood", **s34["stage1"]["stage1"],
                        "pipeline_label": s34["comparison"]["notes"]["1"]},
  "injection_recovery_matrix": s34["injection_test"],
}
R["statement"] = ("The full SH0ES matrix is not blind to a host-level Cepheid term: an injected kappa=0.96e6 is recovered "
                  "at 100% with +-0.21e6. On the real data kappa is consistent with zero and the canonical value is excluded. "
                  "H0=66.65 is origin-dependent and not a full-ladder estimate.")
(HERE / "results").mkdir(exist_ok=True)
json.dump(R, open(HERE / "results" / "tep_cepheid_matrix_profile.json", "w"), indent=1, default=float)
for k in ("profile_published_anchors", "profile_anchors_Vrot_over_sqrt2", "profile_anchors_uncorrected", "profile_diagonal_covariance"):
    p = R[k]; print(f"{k:34s} kappa = {p['kappa_best_1e6']:+.3f} +- {p['kappa_err_1e6']:.3f} e6 ; 0.96e6 excluded at {p['canonical_0p96_excluded_sigma']:.1f} sigma")
print("origin invariance chi2 at 0.96e6:", {k: round(v, 3) for k, v in R["origin_invariance"].items()})
print("66.65 vs screened origin:", round(s04["unified_h0"], 2), round(s04["unified_h0_screened"], 2))
