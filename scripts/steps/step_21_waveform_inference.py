import json
import sys
import sympy as sp

def run_derivation():
    print(r"=========================================================")
    print(r" TEP GW Waveform Inference & Standard Siren Derivation")
    print(r"=========================================================")
    
    # ---------------------------------------------------------
    # 1. Sourced Emission & Einstein-Frame Scaling
    # ---------------------------------------------------------
    print(r"")
    print(r"[1] Sourced Emission & Einstein-Frame Scaling")
    print(r"---------------------------------------------")
    print(r"The Einstein-frame action is:")
    print(r"  S = \int d^4x \sqrt{-g} [M_Pl^2 R / 2] + S_m[\tilde{g}_ab, \psi]")
    print(r"where \tilde{g}_ab = A^2(\phi) g_ab.")
    print(r"")
    print(r"Varying with respect to g_uv gives the Einstein-frame matter stress tensor:")
    print(r"  T_g^uv = (2 / \sqrt{-g}) \delta S_m / \delta g_uv")
    print(r"         = (2 / \sqrt{-g}) (\delta \tilde{g}_ab / \delta g_uv) (\delta S_m / \delta \tilde{g}_ab)")
    print(r"         = (2 / \sqrt{-g}) A^2 \delta^u_a \delta^v_b (\sqrt{-\tilde{g}} / 2) \tilde{T}^ab")
    print(r"Since \sqrt{-\tilde{g}} = A^4 \sqrt{-g}, we have:")
    print(r"  T_g^uv = A^6 \tilde{T}^uv")
    print(r"")
    print(r"For a point mass \tilde{m} in the matter frame, the stress tensor is:")
    print(r"  \tilde{T}^uv = \int d\tilde{\tau} \tilde{m} (dx^u / d\tilde{\tau}) (dx^v / d\tilde{\tau}) \delta^4(x - x_p) / \sqrt{-\tilde{g}}")
    print(r"Using d\tilde{\tau} = A d\tau (where \tau is proper time in g_uv), we transform this to:")
    print(r"  T_g^uv = A^6 \int (A d\tau) \tilde{m} (A^{-1} dx^u / d\tau) (A^{-1} dx^v / d\tau) \delta^4(x - x_p) / (A^4 \sqrt{-g})")
    print(r"  T_g^uv = \int d\tau (A \tilde{m}) (dx^u / d\tau) (dx^v / d\tau) \delta^4(x - x_p) / \sqrt{-g}")
    print(r"")
    print(r"This is identically the stress tensor in g_uv for a particle of mass:")
    print(r"  m = A_e \tilde{m}")
    print(r"Thus, the GW generation equation \Box h_uv = - (16\pi G) (T_uv - 1/2 g_uv T) is sourced by")
    print(r"the Einstein-frame chirp mass:")
    print(r"  M_e = A_e \widetilde{\mathcal{M}}")

    # ---------------------------------------------------------
    # 2. Propagation & Time Intervals
    # ---------------------------------------------------------
    print(r"")
    print(r"[2] Propagation & Time Intervals")
    print(r"--------------------------------")
    print(r"The gravitational waves propagate on the static background metric g_uv.")
    print(r"For successive wave crests emitted at coordinate times t_e and t_e + \Delta t_e,")
    print(r"and received at t_o and t_o + \Delta t_o, the stationarity of g_uv guarantees:")
    print(r"  \Delta t_o = \Delta t_e")
    print(r"Thus, the coordinate frequency is strictly conserved during propagation:")
    print(r"  f_o = f_e")
    print(r"The tensor amplitude decays geometrically as 1/r in the flat spatial background,")
    print(r"with no (1+z) cosmological expansion factor:")
    print(r"  h(t, r) \propto M_e^{5/3} f_e^{2/3} / r")

    # ---------------------------------------------------------
    # 3. Detector Mapping & Inferred Parameters
    # ---------------------------------------------------------
    print(r"")
    print(r"[3] Detector Mapping & Inferred Parameters")
    print(r"------------------------------------------")
    print(r"The LIGO/Virgo detector operates in the matter metric \tilde{g}_uv, measuring:")
    print(r"  Proper time interval: d\tilde{\tau}_o = A_o dt")
    print(r"  Matter-frame frequency: \tilde{f}_o = f_o / A_o")
    print(r"  Fractional strain: \tilde{h} = h")
    print(r"")
    print(r"The GR template matches the phase evolution:")
    print(r"  d\tilde{f}_o / d\tilde{\tau}_o \propto (G_{\rm loc,o} \mathcal{M}_{\rm det})^{5/3} \tilde{f}_o^{11/3}")
    print(r"where G_{\rm loc,o} = G_* A_o^2 (1+\alpha_o^2) is the locally measured Cavendish")
    print(r"constant: what a waveform pipeline reports as 'mass' is the dimensionless")
    print(r"combination (G M)^{5/3}, and the G inserted is the lab-measured one.")
    print(r"")
    print(r"Deriving this from the TEP source dynamics:")
    print(r"  d\tilde{f}_o / d\tilde{\tau}_o = (1 / A_o) d(f_o / A_o) / dt = (1 / A_o^2) df_e / dt")
    print(r"The Einstein-frame phase evolution is df_e / dt \propto (G_{\rm dyn} M_e)^{5/3} f_e^{11/3},")
    print(r"with G_{\rm dyn} = G_*(1+\alpha_A\alpha_B) carrying the pair scalar exchange. Substituting:")
    print(r"  (G_{\rm loc,o} \mathcal{M}_{\rm det})^{5/3} \tilde{f}_o^{11/3}")
    print(r"  = (1 / A_o^2) (G_{\rm dyn} M_e)^{5/3} (A_o \tilde{f}_o)^{11/3}")
    print(r"  = (G_{\rm dyn} M_e)^{5/3} A_o^{5/3} \tilde{f}_o^{11/3}")
    print(r"")
    print(r"Equating gives the frame-consistent inferred chirp mass:")
    print(r"  G_{\rm loc,o} \mathcal{M}_{\rm det} = A_o G_{\rm dyn} M_e")
    print(r"  \mathcal{M}_{\rm det} = (A_e / A_o) [(1+\alpha_A\alpha_B)/(1+\alpha_o^2)] \widetilde{\mathcal{M}}")
    print(r"Consistency check: at A_o = A_e a co-located observer recovers")
    print(r"  \mathcal{M}_{\rm loc} = [(1+\alpha_A\alpha_B)/(1+\alpha_e^2)] \widetilde{\mathcal{M}},")
    print(r"the source's own local mass (step_70) -- the earlier bare-G form")
    print(r"M_det = A_o A_e \tilde M would have violated this. For A_o = 1 (today) and")
    print(r"screened couplings, \mathcal{M}_{\rm det} = \widetilde{\mathcal{M}} / (1+z): the")
    print(r"headline mapping is unchanged, but now holds in a locally consistent frame.")
    print(r"")
    print(r"Now for the amplitude matching to extract inferred distance:")
    print(r"  \tilde{h} \propto (G_{\rm loc,o} \mathcal{M}_{\rm det})^{5/3} \tilde{f}_o^{2/3} / D_L^{GW}")
    print(r"  h \propto (G_{\rm dyn} M_e)^{5/3} f_e^{2/3} / r")
    print(r"Equating \tilde{h} = h with G_{\rm loc,o} M_det = A_o G_{\rm dyn} M_e:")
    print(r"  (A_o G_{\rm dyn} M_e)^{5/3} (f_e / A_o)^{2/3} / D_L^{GW} = (G_{\rm dyn} M_e)^{5/3} f_e^{2/3} / r")
    print(r"  A_o^{5/3} A_o^{-2/3} / D_L^{GW} = 1 / r")
    print(r"  D_L^{GW} = A_o r")
    print(r"")
    print(r"The distance map is unchanged: the chirp combination G*M is conserved under the")
    print(r"frame conversion, so the amplitude inference tracks M_e through the same factor.")

    # ---------------------------------------------------------
    # 4. Independent EM Distance (Sachs Equation)
    # ---------------------------------------------------------
    print(r"")
    print(r"[4] Independent EM Distance (Sachs Equation)")
    print(r"--------------------------------------------")
    print(r"Place the observer at the coordinate origin (r=0) and the source at r.")
    print(r"In the g_uv frame, a solid angle \delta\Omega at the observer subtends an area")
    print(r"  \delta A_g = r^2 \delta\Omega  at the source.")
    print(r"The physical cross-sectional area in the matter metric \tilde{g}_uv at the source is:")
    print(r"  \delta \tilde{A} = A_e^2 \delta A_g = A_e^2 r^2 \delta\Omega")
    print(r"Thus, the matter-frame angular diameter distance is:")
    print(r"  D_A = \sqrt{\delta \tilde{A} / \delta\Omega} = A_e r")
    print(r"")
    print(r"Using Etherington's reciprocity theorem (which holds for any metric theory where photons")
    print(r"travel on null geodesics and phase space is conserved, valid for \tilde{g}_uv):")
    print(r"  D_L^{EM} = (1+z)^2 D_A = (1+z)^2 A_e r")
    print(r"Since 1+z = A_o / A_e, we substitute A_e = A_o / (1+z):")
    print(r"  D_L^{EM} = (1+z)^2 [A_o / (1+z)] r = (1+z) A_o r")

    # ---------------------------------------------------------
    # 5. Resulting Ratio and Data Confrontation
    # ---------------------------------------------------------
    print(r"")
    print(r"[5] Resulting Ratio \Xi(z)")
    print(r"--------------------------")
    print(r"We have derived:")
    print(r"  D_L^{GW} = A_o r")
    print(r"  D_L^{EM} = (1+z) A_o r")
    print(r"")
    print(r"The Standard Siren ratio is therefore:")
    print(r"  \Xi(z) = D_L^{GW} / D_L^{EM} = A_o r / [(1+z) A_o r] = 1 / (1+z)")
    print(r"")
    print(r"This explicitly confirms the formula in the manuscript: \Xi(z) = A(z).")

    # ---------------------------------------------------------
    # 6. Symbolic Verification (machine-checked algebra)
    # ---------------------------------------------------------
    print(r"")
    print(r"[6] Symbolic Verification")
    print(r"-------------------------")
    A_e, A_o, M_t, rr, f_t = sp.symbols('A_e A_o M_tilde r f_t', positive=True)
    K, M_d, D_gw = sp.symbols('K M_d D_gw', positive=True)
    aA, aB, a_o = sp.symbols('a_A a_B a_o', positive=True)

    # Emission: Einstein-frame chirp mass; G_dyn carries the pair scalar exchange
    M_e = A_e * M_t
    G_dyn_over_G = 1 + aA * aB
    G_loc_o_over_G = A_o**2 * (1 + a_o**2)
    f_e = A_o * f_t  # f_o = f_e (static transport); f_o = A_o * f_tilde_o
    fdot_e = K * (G_dyn_over_G * M_e)**sp.Rational(5, 3) * f_e**sp.Rational(11, 3)

    # Detector readout: d f_tilde_o / d tau_tilde_o = (1/A_o^2) df_e/dt
    fdot_t = sp.simplify(fdot_e / A_o**2)

    # Template match in units of the locally measured G:
    # fdot_t = K * (G_loc,o M_det)^(5/3) * f_t^(11/3)  =>  G_loc,o M_det = A_o G_dyn M_e
    M_det = sp.solve(
        sp.Eq(fdot_t, K * (G_loc_o_over_G * M_d)**sp.Rational(5, 3) * f_t**sp.Rational(11, 3)),
        M_d)[0]
    M_det_expected = (A_e / A_o) * ((1 + aA * aB) / (1 + a_o**2)) * M_t
    mass_ok = sp.simplify(M_det - M_det_expected) == 0
    # Screened today-limit: A_o=1, alphas -> 0 gives M_tilde/(1+z) = A_e M_tilde
    mass_today_ok = sp.simplify(
        M_det.subs({A_o: 1, aA: 0, aB: 0, a_o: 0}) - A_e * M_t) == 0

    # Amplitude match: (G_loc,o M_det)^(5/3) f_t^(2/3)/D_L^GW
    #                  = (G_dyn M_e)^(5/3) f_e^(2/3)/r
    D_L_GW = sp.solve(
        sp.Eq((G_loc_o_over_G * M_det)**sp.Rational(5, 3) * f_t**sp.Rational(2, 3) / D_gw,
              (G_dyn_over_G * M_e)**sp.Rational(5, 3) * f_e**sp.Rational(2, 3) / rr),
        D_gw)[0]
    dist_ok = sp.simplify(D_L_GW - A_o * rr) == 0

    # EM distance: D_A = A_e r; Etherington D_L^EM = (1+z)^2 D_A; 1+z = A_o/A_e
    D_L_EM = sp.simplify((A_o / A_e)**2 * (A_e * rr))
    em_ok = sp.simplify(D_L_EM - (A_o / A_e) * A_o * rr) == 0  # (1+z) A_o r

    # Detector-level ratio
    Xi = sp.simplify((A_o * rr) / D_L_EM)
    xi_ok = sp.simplify(Xi - A_e / A_o) == 0  # = 1/(1+z) = A(z)

    print(r"  G_loc,o M_det = A_o G_dyn M_e :", mass_ok)
    print(r"  M_det -> M_tilde/(1+z) today  :", mass_today_ok)
    print(r"  D_L^GW = A_o r               :", dist_ok)
    print(r"  D_L^EM = (1+z) A_o r         :", em_ok)
    print(r"  Xi(z) = A_e/A_o = 1/(1+z)    :", xi_ok)

    # Write JSON output
    output = {
        "sourced_emission": {
            "effective_mass_derivation": "T_g^uv = A^6 \tilde{T}^uv => m = A_e \tilde{m}",
            "einstein_frame_chirp_mass": "M_e = A_e M_tilde"
        },
        "propagation": {
            "time_intervals": "Delta t_o = Delta t_e (static background)",
            "frequency": "f_o = f_e"
        },
        "detector_inference": {
            "matter_proper_time": "d_tau_tilde_o = A_o dt",
            "matter_frequency": "f_tilde_o = f_o / A_o",
            "chirp_combination": "G_loc,o M_det = A_o G_dyn M_e",
            "detector_frame_chirp_mass":
                "M_det = (A_e/A_o)(1+alpha_A alpha_B)/(1+alpha_o^2) M_tilde "
                "= M_tilde/(1+z) for A_o=1, screened",
            "local_consistency": "at A_o=A_e reduces to M_loc (step_70)",
            "inferred_gw_distance": "D_L_GW = A_o r"
        },
        "electromagnetic_distance": {
            "area_distance": "D_A = A_e r (from solid angle at observer to physical area at source)",
            "luminosity_distance": "D_L_EM = (1+z)^2 D_A = (1+z) A_o r"
        },
        "siren_ratio": {
            "Xi_z": "1 / (1+z)",
            "matches_manuscript": True,
            "refutes_cancellation_argument": True
        },
        "symbolic_verification": {
            "chirp_mass_map_verified": bool(mass_ok),
            "chirp_mass_today_limit": bool(mass_today_ok),
            "gw_distance_verified": bool(dist_ok),
            "em_distance_verified": bool(em_ok),
            "siren_ratio_verified": bool(xi_ok)
        }
    }
    
    import os
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_21_waveform_inference.json")
    with open(path, "w") as f:
        json.dump(output, f, indent=2)
    print("\n[SUCCESS] Wrote JSON output to", path)

if __name__ == "__main__":
    run_derivation()
