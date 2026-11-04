#!/usr/bin/env python3
"""Offline order-of-magnitude formalization for issue #16 (stdlib only).

Part A: "black hole = folded Zero sheet" against GR horizon formulas.
Part B: "spiral cut of a drop neck" as a delay-cascade criterion.
Every model choice is labelled; literature inputs are cited in README.md.
"""
import argparse, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
C, G, HBAR = 299792458.0, 6.67430e-11, 1.054571817e-34  # CODATA 2018/2022
KB = 1.380649e-23
M_SUN = 1.98841e30
L_P = math.sqrt(HBAR * G / C**3)
F_P = C**4 / G            # "Planck force", N
LAMBDA = 1.1056e-52       # Planck 2018 cosmological constant, m^-2
# Numerical relativity, equal-mass non-spinning (Scheel et al. 2009, PRD 79 024003)
NR_EQUAL = {'mf': 0.95162, 'chi': 0.68646}
# GW150914 source-frame medians (Abbott et al. 2016, PRL 116 061102)
GW150914 = {'m1': 36.0, 'm2': 29.0, 'mf': 62.0, 'chi': 0.67}


def rs(m):
    if not (math.isfinite(m) and m > 0):
        raise ValueError('mass must be finite positive')
    return 2 * G * m / C**2


def kerr_area(m, chi):
    """Horizon area 8 pi (GM/c^2)^2 (1+sqrt(1-chi^2))."""
    if not 0 <= chi <= 1:
        raise ValueError('0<=chi<=1')
    return 8 * math.pi * (G * m / C**2) ** 2 * (1 + math.sqrt(1 - chi**2))


def sphere_merge_area_loss():
    """Two equal 3D spheres keep total volume: A_f/A_i = 2^(2/3)/2."""
    return 1 - 2 ** (2 / 3) / 2


def horizon_check():
    out = {}
    for name, m in [('sun', M_SUN), ('sgrA', 4.297e6 * M_SUN)]:
        r = rs(m)
        a, v = 4 * math.pi * r**2, 4 / 3 * math.pi * r**3
        out[name] = {
            'r_s_m': r, 'area_m2': a, 'euclid_volume_m3': v,
            'A_over_V_per_m': a / v,           # = 3/r_s, depends on units
            'entropy_bits': a / (4 * L_P**2) / math.log(2),
            # Mc^2 = F_P * r_s / 2: "mass = connectivity" restated exactly
            'Mc2_over_FP_rs_half': m * C**2 / (F_P * r / 2),
            # de Sitter "pressure" scale vs horizon curvature
            'Lambda_rs2': LAMBDA * r**2,
        }
    # area = volume only in units where r_s = 3 length units: not invariant
    out['area_equals_volume_radius_m'] = 3.0
    return out


def merger_check():
    m = 1.0  # each, arbitrary units; ratios are scale-free
    ai = 2 * kerr_area(m * M_SUN, 0)
    af = kerr_area(2 * m * NR_EQUAL['mf'] * M_SUN, NR_EQUAL['chi'])
    g = GW150914
    ai_g = kerr_area(g['m1'] * M_SUN, 0) + kerr_area(g['m2'] * M_SUN, 0)
    af_g = kerr_area(g['mf'] * M_SUN, g['chi'])
    return {
        'sphere_volume_conserving_loss': sphere_merge_area_loss(),
        'gr_equal_mass_area_ratio': af / ai,
        'gr_equal_mass_radiated_fraction': 1 - NR_EQUAL['mf'],
        'gw150914_area_ratio': af_g / ai_g,
        'hawking_area_theorem_holds': af > ai and af_g > ai_g,
        # max area loss allowed by mass-conservation-free bound? none in GR:
        'max_radiated_equal_area_bound': 1 - 1 / math.sqrt(2),
    }


# ---- Part B: drop neck --------------------------------------------------
def neck_radius(tau, sigma=0.072, rho=1000.0):
    """Inviscid pinch r = 0.7 (sigma/rho)^(1/3) tau^(2/3) (Day, Hinch, Lister 1998)."""
    if tau <= 0:
        raise ValueError('tau>0')
    return 0.7 * (sigma / rho) ** (1 / 3) * tau ** (2 / 3)


def cut_length(step, dslow):
    """Delay accumulates dt(L)=L*dslow (dslow = 1/v1-1/v2); cut when dt=step."""
    if step <= 0 or dslow <= 0:
        raise ValueError('positive inputs')
    return step / dslow


def cascade_check(sigma=0.072, rho=1000.0, v=1500.0, dv_rel=1e-3):
    """Model ansatz: two channels differing by dv_rel in speed v (sound in water).

    Levels: step = characteristic reorganization period of that floor.
    Atomic ~ 1 fs (electronic, ~4 eV), molecular ~ 0.2 ps (H-bond libration),
    capillary-hydrodynamic ~ the pinch time itself. Values are illustrative.
    """
    dslow = dv_rel / v
    levels = {'atomic': 1e-15, 'molecular': 2e-13, 'hbond_network': 1e-12}
    res = {}
    for name, step in levels.items():
        lc = cut_length(step, dslow)
        # time-to-pinch at which neck radius equals lc
        tau = (lc / (0.7 * (sigma / rho) ** (1 / 3))) ** 1.5
        res[name] = {'step_s': step, 'L_cut_m': lc, 'tau_match_s': tau,
                     'resolvable': lc > 2.75e-10}  # > water molecule size
    return {'dslow_s_per_m': dslow, 'levels': res}


def angular_momentum(chi_p, grad):
    """Linear ansatz L = chi_p * grad. grad is a polar vector, L is axial, so
    parity allows chi_p != 0 only for a pseudoscalar (chiral medium, rotation,
    magnetic field). Achiral medium: chi_p = 0 and <L> = 0 for either sign."""
    return tuple(chi_p * x for x in grad)


def chirality_check():
    g = (0.0, 0.0, 1.0)
    chiral = angular_momentum(0.3, g)[2], angular_momentum(0.3, (0, 0, -1.0))[2]
    return {'chiral_Lz_for_+g_and_-g': list(chiral),
            'achiral_Lz': angular_momentum(0.0, g)[2],
            'sign_flips_with_gradient': chiral[0] == -chiral[1] != 0}


def close(a, b, rel=1e-9):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(map(close, a, b))
    if isinstance(a, float) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=rel, abs_tol=1e-300)
    return a == b


def results():
    return {'horizons': horizon_check(), 'merger': merger_check(),
            'cascade': cascade_check(), 'chirality': chirality_check()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true', help='compare with results/results.json')
    a = p.parse_args()
    r = json.loads(json.dumps(results()))
    path = ROOT / 'results' / 'results.json'
    if a.check:
        old = json.loads(path.read_text())
        if not close(r, old):
            raise SystemExit('results.json out of date')
        print('results.json matches')
    else:
        path.write_text(json.dumps(r, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(r, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
