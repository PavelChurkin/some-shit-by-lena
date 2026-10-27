#!/usr/bin/env python3
"""Bounded, offline checks for issue #14; standard physics and labelled ansatzes.

No function implements an unspecified holographic operator. See SOURCES.md for
provenance, units, independent inputs and the limits of every comparison.
"""

import argparse
import cmath
import csv
import io
import json
import logging
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
C = 299792458.0  # exact SI, m/s
H = 6.62607015e-34  # exact SI, J s
KB = 1.380649e-23  # exact SI, J/K
E = 1.602176634e-19  # exact SI, J/eV
HBAR = H / (2 * math.pi)


def codata():
    """Read NIST's fixed-width ASCII, retaining source names and SI units."""
    values = {}
    for line in (DATA / 'codata-2022.txt').read_text().splitlines()[11:]:
        if len(line) < 85:
            continue
        name, raw = line[:60].strip(), line[60:85].replace(' ', '')
        try:
            values[name] = float(raw.replace('...', ''))
        except ValueError:
            continue
    return values


CONSTANTS = codata()


def positive(*values):
    if any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('expected finite positive values')


def planck(nu, temperature):
    """B_nu in W m^-2 Hz^-1 sr^-1; frequency, NOT wavelength radiance."""
    positive(nu, temperature)
    x = H * nu / (KB * temperature)
    if x > 700:
        return 0.0
    return 2 * H * nu**3 / C**2 / math.expm1(x)


def spectral_fraction(nu, tau):
    """Illustrative suppression exp[-(nu*tau)^2]-1, not derived from attachment."""
    positive(nu)
    if not math.isfinite(tau) or tau < 0:
        raise ValueError('tau must be nonnegative')
    return math.expm1(-(nu * tau)**2)


def firas_rows():
    return [list(map(float, line.split()))
            for line in (DATA / 'firas_monopole_spec_v1.txt').read_text().splitlines()
            if line.strip() and not line.startswith('#')]


def blackbody_check():
    rows = firas_rows()
    # Column 2 was constructed from a historical 2.725 K blackbody + residual.
    # Recover THAT baseline, avoiding spurious residuals from newer SI constants.
    frequencies = [r[0] * 100 * C for r in rows]
    baseline = [r[1] - r[2] / 1000 for r in rows]  # MJy/sr
    shape = [-1000 * b * (nu / 1e12)**2 for nu, b in zip(frequencies, baseline)]
    denominator = sum((s / r[3])**2 for s, r in zip(shape, rows))
    amplitude = sum(s * r[2] / r[3]**2 for s, r in zip(shape, rows)) / denominator
    chi2 = sum((r[2] / r[3])**2 for r in rows)
    logging.debug('FIRAS diagonal linear ansatz A=%g, chi2=%g', amplitude, chi2)
    return {
        'kind': 'published_residuals_and_illustrative_ansatz',
        'channels': len(rows),
        'frequency_range_GHz': [frequencies[0] / 1e9, frequencies[-1] / 1e9],
        'diagonal_chi2_null': chi2,
        'residual_rms_kJy_sr': math.sqrt(sum(r[2]**2 for r in rows) / len(rows)),
        'linear_amplitude_A_at_1THz': amplitude,
        'diagonal_standard_error_A': 1 / math.sqrt(denominator),
        'formal_confidence_bound': False,
        'reason': 'No covariance matrix, instrument response, or nuisance refit in this table.',
        'illustrative_tau_sweep': [
            {'tau_s': tau, 'diagonal_chi2': sum(
                ((r[2] - 1000 * b * spectral_fraction(nu, tau)) / r[3])**2
                for r, nu, b in zip(rows, frequencies, baseline))}
            for tau in [0, 1e-14, 1e-13, 1e-12]
        ],
        'planck_time_fraction_at_1THz': spectral_fraction(
            1e12, math.sqrt(HBAR * CONSTANTS['Newtonian constant of gravitation'] / C**5)),
        'planck_Bnu_160GHz_2_725K_MJy_sr': planck(160e9, 2.725) / 1e-20,
    }


def photon_probabilities(kind, mean, cutoff=200):
    positive(mean)
    if kind == 'coherent':
        probabilities = [math.exp(-mean)]
        for n in range(1, cutoff + 1):
            probabilities.append(probabilities[-1] * mean / n)
        return probabilities
    if kind == 'thermal':
        return [(mean / (1 + mean))**n / (1 + mean) for n in range(cutoff + 1)]
    raise ValueError('unknown photon state')


def g2(probabilities):
    if not probabilities or any(p < 0 or not math.isfinite(p) for p in probabilities):
        raise ValueError('invalid photon probabilities')
    norm = sum(probabilities)
    positive(norm)
    mean = sum(n * p for n, p in enumerate(probabilities)) / norm
    positive(mean)
    return sum(n * (n - 1) * p for n, p in enumerate(probabilities)) / norm / mean**2


def hom(delay, sigma, visibility):
    """Balanced beam splitter, one photon per input; Gaussian overlap model."""
    positive(sigma)
    if not 0 <= visibility <= 1:
        raise ValueError('visibility must lie in [0, 1]')
    return (1 - visibility * math.exp(-delay**2 / (2 * sigma**2))) / 2


def photon_check():
    return {
        'kind': 'quantum_optics_reference_not_new_correlations',
        'g2_coherent': g2(photon_probabilities('coherent', 2)),
        'g2_thermal': g2(photon_probabilities('thermal', 2)),
        'g2_one_photon': g2([0, 1]),
        'gaussian_HOM_example_sigma_fs': 100,
        'HOM_probabilities': [{'delay_fs': t, 'coincidence_probability': hom(t, 100, 1)}
                              for t in [-300, -100, 0, 100, 300]],
        'experimental_residual_test_performed': False,
        'reason': 'Attachment contains no time tags/counts, detector model or predicted residual.',
    }


def pressure_work_eV(pressure_GPa, delta_volume_A3):
    return pressure_GPa * 1e9 * delta_volume_A3 * 1e-30 / E


def sodium_helium_check():
    return {
        'kind': 'published_counterexample_and_dimensional_check',
        'compound': 'Na2He', 'reported_stability_above_GPa': 113,
        'structure': 'fluorite-type electride',
        'illustrative_pressure_GPa': 113, 'illustrative_delta_volume_A3': -1,
        'illustrative_P_deltaV_eV': pressure_work_eV(113, -1),
        'predicted_phase_boundary': False,
        'reason': 'Symmetry alone gives no Gibbs energies or temperature dependence.',
    }


def ice(temperature, pressure):
    """IAPWS R10-06(2009), Tables 2–4: rho=1/g_p, alpha_V=g_Tp/g_p.

    Only pressure derivatives are needed, so r1, s0 and g00 cancel. Valid box
    is checked here; phase stability must additionally be checked by the caller.
    """
    if not 0 < temperature <= 273.16 or not 0 <= pressure <= 210e6:
        raise ValueError('outside the ice Ih formulation box')
    tt, pt, p0 = 273.16, 611.657, 101325
    tau, dp = temperature / tt, (pressure - p0) / pt
    g0 = [0, 0.655022213658955, -0.189369929326131e-7,
          0.339746123271053e-14, -0.556464869058991e-21]
    t2 = complex(0.337315741065416, 0.335449415919309)
    r21 = complex(-0.557107698030123e-4, 0.464578634580806e-4)
    r22 = complex(0.234801409215913e-10, -0.285651142904972e-10)
    gp0 = sum(k * g0[k] * dp**(k - 1) for k in range(1, 5)) / pt
    r2p = (r21 + 2 * r22 * dp) / pt
    f = ((t2 - tau) * cmath.log(t2 - tau) + (t2 + tau) * cmath.log(t2 + tau)
         - 2 * t2 * cmath.log(t2) - tau**2 / t2)
    ft = -cmath.log(t2 - tau) + cmath.log(t2 + tau) - 2 * tau / t2
    gp = gp0 + tt * (r2p * f).real
    gtp = (r2p * ft).real
    return {'T_K': temperature, 'p_Pa': pressure, 'rho_kg_m3': 1 / gp,
            'alpha_volume_K_inverse': gtp / gp,
            'alpha_linear_average_K_inverse': gtp / gp / 3}


def helium_check():
    mass, density, lambda_temperature = 6.646476989e-27, 146.15, 2.1768
    # Noninteracting Bose gas; not a model of strongly interacting liquid He.
    ideal = 2 * math.pi * HBAR**2 / (mass * KB) * (density / mass / 2.612375348685)**(2 / 3)
    return {'kind': 'ideal_gas_baseline_compared_with_liquid', 'density_kg_m3': density,
            'measured_lambda_K': lambda_temperature, 'ideal_bec_K': ideal,
            'relative_error': (ideal - lambda_temperature) / lambda_temperature,
            'holographic_prediction': False,
            'reason': 'Liquid interactions and superfluid density are not specified.'}


def london_field(depth, penetration):
    positive(penetration)
    if depth < 0:
        raise ValueError('depth must be nonnegative')
    return math.exp(-depth / penetration)


def superconductivity_check():
    return {'kind': 'published_same_element_phase_comparison_and_London_reference',
            'alpha_Ga': {'Tc_K': 1.08, 'space_group': 'Cmca', 'atoms_per_cell': 8},
            'beta_Ga': {'Tc_K': 6.04, 'Tc_uncertainty_K': 0.05,
                        'space_group': 'C2/c', 'atoms_per_cell': 4},
            'Tc_ratio': 6.04 / 1.08,
            'illustrative_penetration_nm': 40,
            'field_ratio_at_5_penetrations': london_field(200e-9, 40e-9),
            'holographic_prediction': False,
            'reason': 'Two phases show structure matters, but do not isolate symmetry from bands/phonons.'}


def clock_rate(speed):
    if not math.isfinite(speed) or not 0 <= speed < C:
        raise ValueError('speed must lie in [0, c)')
    return math.sqrt((1 - speed / C) * (1 + speed / C))


def clock_lag(speed, seconds):
    """Stable t(1-sqrt(1-beta²)) even for laboratory speeds."""
    return seconds * (speed / C)**2 / (1 + clock_rate(speed))


def centrifuge_check():
    radius, rpm = 1.0, 60000
    speed = 2 * math.pi * radius * rpm / 60
    gamma = 100000 * 365.25  # external 100000 yr versus internal one day
    return {'kind': 'special_relativity_reference', 'radius_m': radius, 'rpm': rpm,
            'speed_m_s': speed, 'lag_per_day_s': clock_lag(speed, 86400),
            'proper_time_over_lab': clock_rate(speed),
            'internal_speedup_possible_by_rotation': False,
            'gamma_for_external_100000yr_internal_1day': gamma,
            'one_minus_beta_for_that_gamma_approx': 1 / (2 * gamma**2),
            'acceleration_m_s2': speed**2 / radius}


def bands(diffusivity, spacing, first_position, count):
    """Empirical Liesegang spacing law + x²=2Dt, not simulated agate growth."""
    positive(diffusivity, spacing, first_position)
    if not isinstance(count, int) or not 1 <= count <= 100:
        raise ValueError('use 1..100 bands')
    return [{'band': n + 1, 'position_m': first_position * (1 + spacing)**n,
             'time_s': (first_position * (1 + spacing)**n)**2 / (2 * diffusivity)}
            for n in range(count)]


def agate_check():
    return {'kind': 'empirical_pattern_illustration_not_agate_synthesis',
            'assumed_D_m2_s': 1e-9, 'assumed_spacing_coefficient': 0.1,
            'bands': bands(1e-9, 0.1, 1e-3, 8),
            'experimental_agate_growth_performed': False,
            'reason': 'No solution recipe, measured pattern or kinetic operator in supplied material.'}


def ionization_rows():
    text = (DATA / 'nist-ionization-raw.csv').read_text().split('Notes:')[0]
    result = {}
    for raw in csv.DictReader(io.StringIO(text)):
        row = {k: v.removeprefix('="').removesuffix('"')
               for k, v in raw.items() if k and v}
        if not row.get('At. num', '').isdigit():
            continue
        result[int(row['At. num'])] = row
    return result


def slater_first20(z, configuration):
    """Outer ns/np Slater estimate, only Z<=20; no fitting to radii or energies."""
    if z > 20:
        return None
    cores = {'[He]': '1s2', '[Ne]': '1s2.2s2.2p6', '[Ar]': '1s2.2s2.2p6.3s2.3p6'}
    for core, expanded in cores.items():
        configuration = configuration.replace(core, expanded)
    shells = [(int(n), int(occupancy or '1'))
              for n, orbital, occupancy in re.findall(r'(\d)([spdf])(\d*)', configuration)]
    outer = max(n for n, occupancy in shells)
    same = sum(occupancy for n, occupancy in shells if n == outer) - 1
    previous = sum(occupancy for n, occupancy in shells if n == outer - 1)
    inner = sum(occupancy for n, occupancy in shells if n < outer - 1)
    return z - same * (0.3 if outer == 1 else 0.35) - previous * 0.85 - inner


def atomic_check():
    a0 = CONSTANTS['Bohr radius'] * 1e12
    ionizations = ionization_rows()
    radii = []
    for row in csv.DictReader((DATA / 'covalent-radii.csv').read_text().splitlines()):
        z, n, radius = int(row['Z']), int(row['period']), float(row['single_bond_radius_pm'])
        reference = ionizations.get(z, {})
        zeff = slater_first20(z, reference['Ground Shells (a)']) if z <= 20 else None
        radii.append({'Z': z, 'symbol': row['symbol'], 'period_n_heuristic': n,
                      'reference_radius_pm': radius, 'unscreened_radius_pm': a0 * n**2 / z,
                      'inferred_zeff_from_reference_not_prediction': a0 * n**2 / radius,
                      'slater_zeff': zeff,
                      'slater_radius_pm': a0 * n**2 / zeff if zeff else None,
                      'ionization_eV': float(reference['Ionization Energy (eV)']) if reference else None,
                      'ionization_source_marker': reference.get('Prefix'),
                      'ionization_uncertainty_eV': reference.get('Uncertainty (b) (eV)')})
    alpha = CONSTANTS['fine-structure constant']
    me, mp = CONSTANTS['electron mass'], CONSTANTS['proton mass']
    rydberg = alpha**2 * me * C / (2 * H)
    rh = rydberg / (1 + me / mp)
    return {'kind': 'atomic_reference_table_and_explicit_heuristics',
            'radius_definition': 'Pyykko-Atsumi single-bond covalent fitted radii, not free-atom boundaries',
            'radii': radii, 'radii_count': len(radii), 'ionization_rows_available': len(ionizations),
            'R_infinity_m_inverse': rydberg, 'R_hydrogen_m_inverse': rh,
            'hydrogen_ionization_nonrel_eV': H * C * rh / E,
            'hydrogen_ionization_reference_eV': float(ionizations[1]['Ionization Energy (eV)']),
            'balmer_alpha_vacuum_nm': 1e9 / (rh * (1 / 4 - 1 / 9)),
            'H_levels_eV': {str(n): -H * C * rh / E / n**2 for n in range(1, 6)},
            'reason': 'alpha, mass, h and atomic shell rules are inputs, not operator outputs.'}


def constants_check():
    alpha = CONSTANTS['fine-structure constant']
    me = CONSTANTS['electron mass energy equivalent in MeV']
    mp = CONSTANTS['proton mass energy equivalent in MeV']
    mn = CONSTANTS['neutron mass energy equivalent in MeV']
    return {'kind': 'measured_targets_and_QED_one_loop_reference',
            'alpha_reference': alpha, 'inverse_alpha_reference': 1 / alpha,
            'inverse_alpha_rounded_issue_relative_error': (137.036 * alpha - 1),
            'electron_MeV': me, 'proton_MeV': mp, 'neutron_MeV': mn,
            'proton_over_neutron': mp / mn,
            'electron_over_proton': me / mp,
            'proton_over_electron_reference': CONSTANTS['proton-electron mass ratio'],
            'anomaly_reference': CONSTANTS['electron mag. mom. anomaly'],
            'g_abs_reference': abs(CONSTANTS['electron g factor']),
            'anomaly_one_loop': alpha / (2 * math.pi),
            'g_abs_one_loop': 2 + alpha / math.pi,
            'alpha_derived_from_geometry': False,
            'reason': 'No dimensionless coupling or mass-generating eigenproblem is supplied.'}


def beta_shape(kinetic_MeV, endpoint_MeV):
    if not 0 <= kinetic_MeV <= endpoint_MeV:
        raise ValueError('kinetic energy must lie within endpoint')
    mass = CONSTANTS['electron mass energy equivalent in MeV']
    momentum = math.sqrt(kinetic_MeV * (kinetic_MeV + 2 * mass))
    return momentum * (kinetic_MeV + mass) * (endpoint_MeV - kinetic_MeV)**2


def decay_check():
    q = (CONSTANTS['neutron mass energy equivalent in MeV']
         - CONSTANTS['proton mass energy equivalent in MeV']
         - CONSTANTS['electron mass energy equivalent in MeV'])
    alpha_q = 5.40745  # Po-210, DDEP total Q, NOT alpha kinetic energy
    return {'kind': 'energy_conservation_and_kinematics_reference',
            'Po210_Q_MeV_input': alpha_q,
            'Po210_alpha_MeV_nonrel_mass_number_approx': alpha_q * 206 / 210,
            'Po210_alpha_MeV_reference': 5.30433,
            'neutron_Q_MeV': q,
            'beta_spectrum_no_Coulomb_recoil_radiative_corrections': [
                {'electron_kinetic_MeV': q * i / 20, 'arbitrary_density': beta_shape(q * i / 20, q)}
                for i in range(21)],
            'Co60_gamma_keV_inputs': [1173.228, 1332.492],
            'Co60_gamma_wavelength_pm': [H * C / (energy * 1000 * E) * 1e12
                                        for energy in [1173.228, 1332.492]],
            'operator_predicts_nuclear_energies': False,
            'reason': 'Q values/level gaps and interactions are inputs; kinematics is not a nuclear model.'}


def gravity_check():
    bodies = []
    for name, gm, radius in [('Earth', 3.986004418e14, 6.371e6),
                             ('Sun', 1.3271244e20, 6.957e8)]:
        compactness = 2 * gm / (radius * C**2)
        bodies.append({'body': name, 'GM_m3_s2_input': gm, 'radius_m_input': radius,
                       'potential_m2_s2': -gm / radius, 'surface_g_m_s2': gm / radius**2,
                       'schwarzschild_radius_m': 2 * gm / C**2,
                       'compactness_Rs_over_R': compactness,
                       'Kretschmann_m_inverse4': 48 * gm**2 / (C**4 * radius**6)})
    gm, earth, orbit = 3.986004418e14, 6.371e6, 26.56e6
    grav = gm * (1 / earth - 1 / orbit) / C**2 * 86400 * 1e6
    kinetic = -gm / (2 * orbit * C**2) * 86400 * 1e6
    return {'kind': 'Newton_Schwarzschild_weak_field_reference', 'bodies': bodies,
            'gps_gravity_us_day': grav, 'gps_kinetic_us_day': kinetic,
            'gps_net_us_day': grav + kinetic,
            'gps_approximations': 'circular orbit, spherical nonrotating Earth, no multipoles',
            'reason': 'GM and field equations are external inputs, not derived closure parameters.'}


def kerr_area_factor(spin):
    if not 0 <= spin <= 1:
        raise ValueError('dimensionless Kerr spin must lie in [0, 1]')
    return (1 + math.sqrt(1 - spin**2)) / 2


def horizon_check():
    gm_sun = 1.3271244e20
    gm, spin = 62 * gm_sun, 0.67
    radius = 2 * gm / C**2
    area = 4 * math.pi * radius**2 * kerr_area_factor(spin)
    planck_area = HBAR * CONSTANTS['Newtonian constant of gravitation'] / C**3
    return {'kind': 'classical_horizon_countercheck_not_invariant_interior_volume',
            'input_masses_solar': [36, 29, 62], 'final_spin_input': spin,
            'area_ratio_nonrotating': 62**2 / (36**2 + 29**2),
            'area_ratio_final_spin_0_67_initial_spins_zero':
                62**2 * kerr_area_factor(spin) / (36**2 + 29**2),
            'euclidean_volume_proxy_ratio': 62**3 / (36**3 + 29**3),
            'final_entropy_over_kB': area / (4 * planck_area),
            'Schwarzschild_Hawking_temperature_62_solar_K': HBAR * C**3 / (8 * math.pi * KB * gm),
            'Kerr_horizon_angular_velocity_rad_s':
                spin * C**3 / (2 * gm * (1 + math.sqrt(1 - spin**2))),
            'reason': 'Area increases in this example; Euclidean sphere volume is not a GR interior invariant.'}


def galaxy_check():
    # A declared illustrative flat velocity, not a fitted observational data set.
    v, g = 200e3, CONSTANTS['Newtonian constant of gravitation']
    return {'kind': 'flat_curve_required_mass_profile_illustration',
            'assumed_velocity_km_s': v / 1000,
            'enclosed_mass_ratio_20kpc_over_10kpc': 2,
            'required_density_at_10kpc_kg_m3': v**2 / (4 * math.pi * g * (10 * 3.085677581491367e19)**2),
            'reason': 'Calling it connectivity supplies neither a stress tensor nor lensing predictions.'}


def operator_check():
    # Countermodels to identifiability: every example obeys a/tau=c and relaxation
    # of mismatch. Free dimensionless couplings remain unrestricted by that rule.
    return {'kind': 'explicit_nonuniqueness_witness_not_a_proposed_physical_theory',
            'admissible_examples': [
                {'length_m': a, 'step_s': a / C, 'signal_speed_m_s': C,
                 'illustrative_alpha': coupling, 'dimensionless_relaxation_rate': rate,
                 'mismatch_at_dimensionless_time_1': math.exp(-rate)}
                for a, coupling, rate in [(1e-20, 0.01, 1), (1e-30, 0.02, 2), (1e-35, 0.03, 3)]
            ],
            'unique_alpha_or_lambda_temperature': False,
            'reason': 'A finite speed and convergence of mismatch leave scale and coupling free.'}


def run_all():
    return {
        'blackbody': blackbody_check(), 'photons': photon_check(),
        'sodium_helium': sodium_helium_check(),
        'ice': {'kind': 'IAPWS_reference_EOS',
                'rows_at_101325Pa': [ice(t, 101325) for t in [100, 150, 200, 250, 270]],
                'holographic_prediction': False},
        'constants': constants_check(), 'helium': helium_check(),
        'superconductivity': superconductivity_check(), 'centrifuge': centrifuge_check(),
        'agate': agate_check(), 'atoms': atomic_check(), 'decays': decay_check(),
        'gravity': gravity_check(), 'horizons': horizon_check(),
        'galaxies': galaxy_check(), 'operator': operator_check(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--write', action='store_true', help='regenerate tracked results')
    action.add_argument('--check', action='store_true', help='fail if tracked results are stale')
    parser.add_argument('--verbose', action='store_true', help='enable diagnostic logging')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING)
    result = run_all()
    path = ROOT / 'results/results.json'
    if args.write:
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
        print('Wrote', path)
    elif args.check:
        saved = json.loads(path.read_text())
        # Compare floating-point values at a precision portable across libm builds.
        def equivalent(a, b):
            if isinstance(a, float) and isinstance(b, (float, int)):
                return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-25)
            if isinstance(a, dict) and isinstance(b, dict):
                return a.keys() == b.keys() and all(equivalent(a[k], b[k]) for k in a)
            if isinstance(a, list) and isinstance(b, list):
                return len(a) == len(b) and all(equivalent(x, y) for x, y in zip(a, b))
            return a == b
        if not equivalent(result, saved):
            parser.exit(1, 'Stored results differ; run checks.py --write and review.\n')
        print('All 15 stored calculation groups match.')
    else:
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
