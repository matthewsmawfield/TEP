#!/usr/bin/env python3
"""
step_22_tilted_congruence.py
==========================
Re-derive the zero-volume-expansion condition for a TILTED matter
congruence in the inhomogeneous temporal landscape.  Every prior bound
(steps 16-21) assumed the volume-normal congruence (v=0), forcing all
compensation onto the disformal Q.  Here the congruence itself may flow.

Matter metric:  g~_ab = A^2(phi) * g^_ab,  g^_ab = eta_ab + b~ phi_a phi_b
Congruence:     u~^a = Gamma/A * (u^_n^a + v^i s^i_a)   (tilt v in the
                g^-normal frame)

Identity used (exact, conformal):
    theta~ = A^{-1} [ theta^ + 3 u^.grad ln A ]

For a phi(t,x) field (plane symmetry), theta^ is computed by direct
divergence of the tilted congruence.  Tests:

  T1  v=0 reproduces the volume-normal identity: theta~=0 demands the
      disformal spatial structure (the A^6(1+Q)=C law).
  T2  b~=0 (pure conformal) + converging coordinate flow v = -Pi_A x:
      theta~=0 achievable with NO disformal structure -- but the observed
      z(L) between congruence members cancels at first order (the
      <H>=theta~/3 statement).  Confirms the tilted route is not a free
      lunch in the homogeneous limit.
  T3  Inhomogeneous landscape phi(t,x) = drift + structured lumps +
      ensemble of moving matter: solve for the flow that keeps <theta~>=0
      statistically, then measure the accumulated photon transport z(L)
      -- the 'local+local+local' channel on the true congruence.

Units c=1.  All fields prescribed (kinematic test); Einstein/KG residual
gates deferred to the coupled solve.
"""
import json
import numpy as np
import os

EPS = 1e-5


def ghat(phi_t, phi_x, b):
    """g^_ab = eta + b phi_a phi_b (2D, t-x)."""
    return np.array([[-1.0 + b*phi_t*phi_t, b*phi_t*phi_x],
                     [b*phi_t*phi_x,        1.0 + b*phi_x*phi_x]])


def normal_frame(gh):
    """Unit timelike normal u_n^a (normalise -gh(u,u)=1) and unit
    spatial vector s^a orthogonal to it."""
    # normal one-form n_a ~ delta_at (constant-t slices)
    gi = np.linalg.inv(gh)
    nt2 = -gi[0, 0]                      # g^{ab} n_a n_b, n_a=(1,0)
    un = np.array([-gi[0, 0], -gi[0, 1]]) / np.sqrt(nt2)   # u^a = g^{ab} n_b
    # spatial unit vector orthogonal: s^a prop to (u_1, -u_0)
    u_low = gh @ un
    s = np.array([u_low[1], -u_low[0]])
    s /= np.sqrt(abs(gh @ s @ s))
    if gh @ s @ s < 0:                   # keep spacelike
        s = -s
    return un, s


def divergence(F_t, F_x, gh_vol, t, x, i, j):
    """(1/sqrt|g|) d_a (sqrt|g| F^a) with central differences."""
    dt = t[1]-t[0]; dx = x[1]-x[0]
    dT = (gh_vol[i+1, j]*F_t[i+1, j] - gh_vol[i-1, j]*F_t[i-1, j])/(2*dt)
    dX = (gh_vol[i, j+1]*F_x[i, j+1] - gh_vol[i, j-1]*F_x[i, j-1])/(2*dx)
    return (dT + dX)/gh_vol[i, j]


def theta_tilde(phi_fn, dphi_fn, b, v_field, t, x):
    """theta~ on a grid for congruence u~ = Gamma/A (u_n + v s)."""
    NT, NX = len(t), len(x)
    th = np.zeros((NT, NX))
    # build congruence and volume factor on the grid
    UT = np.zeros((NT, NX)); UX = np.zeros((NT, NX))
    vol = np.zeros((NT, NX))
    lnA_t = np.zeros((NT, NX)); lnA_x = np.zeros((NT, NX))
    th_hat = np.zeros((NT, NX))
    udot_lnA = np.zeros((NT, NX))
    for i, ti in enumerate(t):
        for j, xj in enumerate(x):
            pt, px = dphi_fn(ti, xj)
            gh = ghat(pt, px, b)
            un, s = normal_frame(gh)
            v = v_field(ti, xj)
            Gam = 1.0/np.sqrt(1.0 - v*v)
            ua = Gam*(un + v*s)
            UT[i, j], UX[i, j] = ua
            vol[i, j] = np.sqrt(abs(np.linalg.det(gh)))
            # u^.grad ln A = u^a d_a ln A = -u^a d_a phi  (A = e^{-phi})
            udot_lnA[i, j] = -(ua[0]*pt + ua[1]*px)
    # divergence of u^ (in g^)
    for i in range(1, NT-1):
        for j in range(1, NX-1):
            th_hat[i, j] = divergence(UT, UX, vol, t, x, i, j)
    A = np.exp(-np.array([[phi_fn(ti, xj) for xj in x] for ti in t]))
    # 2D toy: one spatial dim -> conformal coefficient is (d-1)=1, not 3
    th = (th_hat + 1.0*udot_lnA)/A        # theta~ = A^-1(theta^ + n u.d lnA)
    return th, th_hat, udot_lnA


# ---------- fields -------------------------------------------------------
def drift_field(Pi):
    return (lambda t, x: -Pi*t,
            lambda t, x: (-Pi, 0.0))


def wave_field(Phi0, Pi, b_target, center):
    """zero-volume converging wave (step-18 type) + drift."""
    def f(u):
        # f' = sqrt((1 - e^{-6 f})/b)  -- pre-integrated numerically
        pass
    return None


def lumpy_field(Pi, lumps):
    """phi = -Pi t + sum_k A_k exp(-(x-x_k)^2/2 l_k^2) * (1 + e_k t)"""
    def phi(t, x):
        v = -Pi*t
        for (A_k, x_k, l_k, e_k) in lumps:
            v += A_k*(1+e_k*t)*np.exp(-(x-x_k)**2/(2*l_k*l_k))
        return v
    def dphi(t, x):
        pt = -Pi; px = 0.0
        for (A_k, x_k, l_k, e_k) in lumps:
            g = np.exp(-(x-x_k)**2/(2*l_k*l_k))
            pt += A_k*e_k*g
            px += A_k*(1+e_k*t)*g*(-(x-x_k)/l_k**2)
        return pt, px
    return phi, dphi


def main():
    t = np.linspace(-10, 110, 400)
    x = np.linspace(-100, 100, 400)
    res = {}

    # ---- T1: volume-normal check -------------------------------------
    # pure drift, b=0, v=0: theta~ should equal 3*Pi/A != 0 (expansion)
    Pi = 0.02
    phi, dphi = drift_field(Pi)
    th, thh, ud = theta_tilde(phi, dphi, 0.0, lambda tt, xx: 0.0, t, x)
    res["T1_drift_B0_v0"] = {
        "theta_tilde_rms": float(np.sqrt(np.mean(th[5:-5, 5:-5]**2))),
        "expected": "3*Pi/A ~ nonzero -> drift alone IS expansion"}

    # pure drift, b>0 CANNOT fix (needs spatial gradient); verify
    th, _, _ = theta_tilde(phi, dphi, 5.0, lambda tt, xx: 0.0, t, x)
    res["T1_drift_B5_v0"] = {
        "theta_tilde_rms": float(np.sqrt(np.mean(th[5:-5, 5:-5]**2))),
        "expected": "still nonzero -- disformal needs SPATIAL gradient"}

    # ---- T2: conformal + converging flow ------------------------------
    # In the 2D toy: theta~ = A^-1(theta^ + u.dlnA); pure drift gives
    # u.dlnA = -Pi (lnA = phi_hat... A=e^{-phi}, phi=-Pi t -> lnA = Pi t)
    # need div v = -Pi  ->  v = -Pi x  (converging coordinate flow)
    for Hfac in [1.0]:
        vfield = lambda tt, xx: np.clip(-Pi*Hfac*xx, -0.9, 0.9)
        th, thh, ud = theta_tilde(phi, dphi, 0.0, vfield, t, x)
        res["T2_converging_flow"] = {
            "H_over_Pi": Hfac,
            "theta_tilde_rms": float(np.sqrt(np.mean(th[5:-5, 5:-5]**2))),
            "theta_tilde_mean": float(np.mean(th[5:-5, 5:-5])),
            "expected": "~0 if flow divergence compensates drift"}

    # ---- T3: landscape-correlated flow ---------------------------------
    # TEP picture: galaxies drift toward slow-time wells (phi maxima):
    # v = kappa * d_x phi.  Volume conservation then requires
    # <v d_x phi> = Pi  (flow-gradient correlation carries the budget,
    # not a coherent converging pattern).
    lumps = [(0.4, -30.0, 8.0, 0.001), (0.3, 20.0, 6.0, -0.002),
             (0.5, 60.0, 10.0, 0.0015)]
    phi3, dphi3 = lumpy_field(Pi, lumps)
    # measure <(d_x phi)^2> to pick kappa = Pi/<px^2>
    px2 = np.mean([dphi3(50.0, xj)[1]**2 for xj in x])
    for kap in [0.0, Pi/px2, 2*Pi/px2]:
        def vf(tt, xx, kap=kap):
            _, px = dphi3(tt, xx)
            return np.clip(kap*px, -0.9, 0.9)
        th, thh, ud = theta_tilde(phi3, dphi3, 0.0, vf, t, x)
        res[f"T3_kappa_{kap:.3e}"] = {
            "theta_tilde_rms": float(np.sqrt(np.mean(th[10:-10, 10:-10]**2))),
            "theta_tilde_mean": float(np.mean(th[10:-10, 10:-10])),
            "note": "v = kappa d_x phi: flow toward slow-time wells"}

    # ---- T4: observed z(L) on the zero-expansion congruence ------------
    # photon geodesic in g~ = A^2 eta (B=0 branch) between congruence
    # members: emitter at x=L with v_E = -Pi L, receiver at 0.
    # conserved covariant k_t (eta-static):  w~ = -k.u~ = Gamma(1+v)k_c/A
    # (signs resolved in code by direct evaluation)
    A = lambda ph: np.exp(-ph)
    Ls = np.linspace(0, 60, 13)
    zvals = []
    for L in Ls:
        PiL = Pi*L
        vE = -PiL                      # converging flow of T2
        GamE = 1.0/np.sqrt(1-vE*vE) if abs(vE) < 1 else np.nan
        # photon emitted at t_E = -L (travel time), moving -x
        # w_E = k_c Gamma (1 - vE_photon) / A_E  with photon k_x = -k_c
        # source moves -x (vE<0): 1 - vE*(n)=1 - (-PiL)*(-1)...
        # direct: -k.u = -(k_t u^t + k_x u^x); k_t=-kc, k_x=-kc,
        # u^t = Gam/A, u^x = Gam v/A  ->  w = kc Gam (1+v)/A
        phiE = -Pi*(-L)                # phi at emission event (x=L,t=-L)
        phiR = 0.0                     # at reception (x=0,t=0)
        wE = GamE*(1.0+vE)/A(phiE)
        wR = 1.0/A(phiR)
        zvals.append(wE/wR - 1.0)
    res["T4_z_of_L_B0_converging"] = {
        "L": Ls.tolist(), "z": [float(z) for z in zvals],
        "note": "conformal drift redshift vs converging-flow Doppler: "
                "first-order cancellation is the <H>=theta/3 statement"}

    # ---- T5: same but ensemble inhomogeneous ---------------------------
    # redshift through the lumpy landscape with the T3 congruence:
    # integrate photon transport numerically (RK2 on g~ = A^2 eta +
    # b phi_a phi_b) -- reuse of step-21 machinery inline.
    def gtil(t_, x_, b_):
        pt, px = dphi3(t_, x_)
        A2 = np.exp(-2*phi3(t_, x_))
        return A2*np.array([[-1+b_*pt*pt, b_*pt*px],
                            [b_*pt*px, 1+b_*px*px]])

    def chris(t_, x_, b_):
        gi = np.linalg.inv(gtil(t_, x_, b_))
        dg = np.array([(gtil(t_+EPS, x_, b_)-gtil(t_-EPS, x_, b_))/(2*EPS),
                       (gtil(t_, x_+EPS, b_)-gtil(t_, x_-EPS, b_))/(2*EPS)])
        return 0.5*np.einsum("ad,dbc->abc", gi,
               dg.transpose(1, 0, 2)+dg.transpose(2, 0, 1)-dg)

    def cross(b_, vflow, L=80.0, dl=0.005):
        """photon from x=-L emitter (moving with flow) to x=0 receiver."""
        gg = gtil(0.0, -L, b_)
        aa_, bb_, cc_ = gg[0, 0], 2*gg[0, 1], gg[1, 1]
        disc = bb_*bb_-4*aa_*cc_
        if disc <= 0: return None
        rts = [(-bb_+np.sqrt(disc))/(2*aa_), (-bb_-np.sqrt(disc))/(2*aa_)]
        kt = max(rts)                    # future-directed root (k^t > 0)
        K = np.array([kt, 1.0]); X = np.array([0.0, -L])
        while X[1] < 0 and X[0] < 4*L and X[0] > -4*L:
            G = chris(*X, b_)
            Km = K-0.5*dl*np.einsum("abc,b,c->a", G, K, K)
            Xm = X+0.5*dl*K
            G2 = chris(*Xm, b_)
            K = K-dl*np.einsum("abc,b,c->a", G2, Km, Km)
            X = X+dl*Km
        ge = gtil(0.0, -L, b_); gr = gtil(X[0], 0.0, b_)
        # emitter/receiver on the flow congruence
        vE = vflow(0.0, -L); vR = vflow(X[0], 0.0)
        # matter-frame freq: -k.u with u^a = Gam/A (u_n + v s) in g~ frame
        def omega(g_, K_, v_):
            # normal frame of g~ (which includes A): reuse normal_frame
            un, s = normal_frame(g_)
            Gam = 1.0/np.sqrt(1-v_*v_)
            ua = Gam*(un + v_*s)
            return -(g_ @ K_) @ ua
        wE = omega(ge, np.array([kt, 1.0]), vE)
        wR = omega(gr, K, vR)
        return {"z": wE/wR-1.0, "t_arr": X[0]}

    out = {}
    kap = Pi/px2
    for b_ in [0.0, 5.0, 20.0]:
        vf = lambda tt, xx: np.clip(kap*dphi3(tt, xx)[1], -0.9, 0.9)
        out[f"b_{b_}"] = cross(b_, vf)
    res["T5_z_inhomogeneous"] = out

    res["interpretation"] = {
        "T1": "v=0: drift needs disformal structure (old identity)",
        "T2": "homogeneous converging flow gives theta~=0 but kills "
              "linear z(L) -- pure relabeling, no Hubble law",
        "T3": "inhomogeneous landscape: flow toward/through gradients "
              "shares the burden statistically",
        "T5": "the question that matters: transport z on the true "
              "congruence in a lumpy landscape"}
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_22_tilted_congruence.json"),
              "w") as f:
        json.dump(res, f, indent=2)
    for k, v in res.items():
        print(k, ":", json.dumps(v, indent=1)[:900], "\n")


if __name__ == "__main__":
    main()
