#!/usr/bin/env python3
"""Coupled interior: the lapse-tracking scalar is allowed to move the metric.

Exterior seed, source-free, on a fixed background:
    U = -ln(1-2M/r),    flux P = r^2 F U' = -2M.
The invariant F (U')^2 diverges as F -> 0. Here m(r) and the redshift
factor are integrated with that stress tensor, inward from a mild
exterior point. Independent variable is the areal radius while F stays
positive. The integration stops if F reaches a floor or if a turning
point F' = 0 appears, where the accumulated gradient energy starts
driving the mass aspect down and F back up (the well reopens).

Canonical sector first. A second pass multiplies the radial kinetic
term by P_X = 1 + I/I_0, the leading piece of the corpus completion
X|X|/Lambda^4, and asks whether the invariant then stays finite while
the matter-frame clock A * sqrt(F) still falls.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

# G = 1. Einstein equation uses m' = 4 pi r^2 rho.
FOUR_PI = 4.0 * np.pi
M0 = 1.0


def pack(r, m, delta, u, up, px):
    F = 1.0 - 2.0 * m / r
    I = F * up * up  # g^{rr} (u')^2
    A = float(np.exp(-u))
    # Matter-frame clock for a static observer: the gravitational lapse is
    # e^delta * sqrt(F) for ds^2 = -e^{2d}F dt^2 + F^{-1}dr^2 + r^2 dO^2.
    clock = A * float(np.exp(delta)) * float(np.sqrt(max(F, 0.0)))
    # Kretschmann of the spherical metric is dominated near a horizon
    # by the mass aspect; report the scalar invariant and 48 m^2/r^6.
    K_schw = 48.0 * m * m / r ** 6
    return {
        "r": float(r),
        "m": float(m),
        "F": float(F),
        "u": float(u),
        "A": A,
        "matter_clock": float(clock),
        "lapse2": float(np.exp(2.0 * delta) * F),
        "invariant": float(I),
        "P_X": float(px),
        "K48": float(K_schw),
    }


def make_rhs(I0):
    """I0 = None -> canonical P_X = 1. Else P_X = 1 + I/I0."""

    def rhs(r, y):
        # Inward: the solver will be given a decreasing interval, so the
        # returned derivatives are with respect to r (not to -r).
        m, delta, u, up = y
        F = 1.0 - 2.0 * m / r
        if F < 1.0e-8:
            return [0.0, 0.0, 0.0, 0.0]
        I = F * up * up
        px = 1.0 if I0 is None else 1.0 + I / I0
        # rho = 1/2 P_X * I   for V = 0 and the radial kinetic term.
        # p_r = rho for a purely radial canonical gradient (V=0).
        rho = 0.5 * px * I
        dm = FOUR_PI * r * r * rho
        # δ' = 4 pi r (rho + p_r) / F * F / I * up^2 = 4 pi r px up^2
        # since rho + p_r = px * I = px F up^2, divided by F.
        ddelta = FOUR_PI * r * px * up * up
        # (r^2 e^δ F px u')' = 0  =>  flux constant only if we do not
        # differentiate px. Full product rule:
        # d/dr (F px up) = - (2/r + δ') F px up
        # F' = -2 (dm r - m) / r^2
        Fp = -2.0 * (dm * r - m) / r ** 2
        # Let w = F * px * up
        # w' = -(2/r + δ') w
        # F' px up + F px' up + F px up' = -(2/r+δ') w
        # px = 1 + F up^2 / I0, px' depends on up'.
        if I0 is None:
            # F up' + Fp up = -(2/r+δ') F up
            up_p = (-(2.0 / r + ddelta) * F * up - Fp * up) / F
        else:
            # px = 1 + (F up^2)/I0
            # d(px)/d(up) = 2 F up / I0
            # w = F px up
            # dw/dup excluding up' : treat F' separately
            # w' = Fp*px*up + F*(dpx/dr)*up + F*px*up'
            # dpx/dr = (Fp up^2 + F * 2 up up') / I0
            # Collect up':
            # F * (2 F up^2 / I0) * up / wait
            # F * (F * 2 up * up' / I0) * up + F px up'
            # = F^2 * 2 up^2 / I0 * up' + F px up'
            coeff = F * px + (2.0 * F * F * up * up / I0)
            # known pieces of w' that must equal -(2/r+δ') w
            # Fp*px*up + F*(Fp up^2 / I0)*up + coeff*up' = -(2/r+δ') w
            known = Fp * px * up + F * (Fp * up * up / I0) * up
            target = -(2.0 / r + ddelta) * (F * px * up)
            up_p = (target - known) / coeff
        return [dm, ddelta, up, up_p]

    return rhs


def integrate(I0, r_start=80.0, r_floor=1.8, flux=-2.0):
    # flux = -2M is the full lapse-tracking charge. Its gradient energy
    # outside a moderate radius is comparable to M, so the interior mass
    # must be reduced by that energy before the inward march begins.
    P = flux * M0
    E_out = 2.0 * np.pi * P * P / r_start
    m0 = max(M0 - E_out, 0.2)
    F0 = 1.0 - 2.0 * m0 / r_start
    u0 = -np.log(max(F0, 1e-8))
    up0 = P / (r_start ** 2 * F0)
    y0 = [m0, 0.0, u0, up0]

    def hit_floor(r, y):
        return (1.0 - 2.0 * y[0] / r) - 1.0e-4

    hit_floor.terminal = True
    hit_floor.direction = -1

    def turned(r, y):
        # F' = 0: the mass aspect turns around, i.e. the radius where the
        # scalar's accumulated gradient energy starts driving F back up.
        m = y[0]
        dm = make_rhs(I0)(r, y)[0]
        return -2.0 * (dm * r - m) / r ** 2

    turned.terminal = True
    turned.direction = -1

    sol = solve_ivp(
        lambda r, y: make_rhs(I0)(r, y),
        (r_start, r_floor),
        y0,
        method="Radau",
        rtol=1.0e-6,
        atol=1.0e-8,
        max_step=0.25,
        events=(hit_floor, turned),
        dense_output=True,
    )
    samples = []
    if sol.sol is not None:
        for rf in (40.0, 10.0, 6.0, 4.0, 3.0, 2.5, 2.2, 2.1, 2.05, 2.0, 1.95, 1.9):
            if rf <= sol.t[-1] or rf > r_start:
                continue
            yy = sol.sol(rf)
            px = 1.0 if I0 is None else 1.0 + (1 - 2 * yy[0] / rf) * yy[3] ** 2 / I0
            samples.append(pack(rf, yy[0], yy[1], yy[2], yy[3], px))
    end = pack(sol.t[-1], *sol.y[:, -1], 1.0 if I0 is None else np.nan)
    return {
        "I0": I0,
        "flux": flux,
        "success": bool(sol.success),
        "message": sol.message,
        "r_end": float(sol.t[-1]),
        "end": end,
        "samples": samples,
    }


def main():
    runs = [integrate(None, flux=f) for f in (-2.0, -0.3, -0.05)]
    runs.append(integrate(1.0, flux=-0.05))  # kinetic-cap reference run
    for run in runs:
        e = run["end"]
        print(
            f"flux={run.get('flux')} r_end={run['r_end']:.4f} F={e['F']:.3e} "
            f"u={e['u']:.3f} A={e['A']:.3e} clock={e['matter_clock']:.3e} "
            f"I={e['invariant']:.3e} m={e['m']:.3f}"
        )
    out = {
        "step": "step_48_coupled_interior",
        "seed": "U=-ln(1-2M/r) profile initialized at r_start=80, then m and u integrated inward together",
        "runs": runs,
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "step_48_coupled_interior.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
