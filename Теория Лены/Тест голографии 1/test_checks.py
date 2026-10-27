"""Independent physical limits and reference values; python3 -m unittest discover."""

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import checks


class RadiationTests(unittest.TestCase):
    def test_rayleigh_jeans_limit_and_wien_tail(self):
        nu, temperature = 1e6, 300.0
        expected = 2 * checks.KB * temperature * nu**2 / checks.C**2
        self.assertAlmostEqual(checks.planck(nu, temperature) / expected, 1, places=6)
        self.assertEqual(checks.planck(1e18, 1), 0)

    def test_stefan_boltzmann_integral(self):
        # Integrate the dimensionless spectrum independently, including its tail.
        step = 0.005
        integral = sum((i * step)**3 / math.expm1(i * step)
                       for i in range(1, 10000)) * step
        self.assertAlmostEqual(integral, math.pi**4 / 15, places=7)

    def test_firas_column_units(self):
        rows = checks.firas_rows()
        self.assertEqual(len(rows), 43)
        self.assertEqual(rows[0][2:4], [5.0, 14.0])  # kJy/sr, not MJy/sr
        self.assertAlmostEqual(rows[0][0] * 100 * checks.C / 1e9, 68.052887966)
        result = checks.blackbody_check()
        self.assertGreater(result['diagonal_chi2_null'], 1)
        self.assertLess(result['diagonal_chi2_null'], 100)
        self.assertFalse(result['formal_confidence_bound'])

    def test_correction_zero_and_small_argument(self):
        self.assertEqual(checks.spectral_fraction(1e12, 0), 0)
        self.assertAlmostEqual(checks.spectral_fraction(1e12, 1e-20) / -1e-16, 1)

    def test_photon_factorial_moments(self):
        self.assertEqual(checks.g2([0, 1]), 0)
        self.assertAlmostEqual(checks.g2(checks.photon_probabilities('coherent', 2)), 1)
        self.assertAlmostEqual(checks.g2(checks.photon_probabilities('thermal', 2)), 2)

    def test_hom_distinguishability_and_probability(self):
        self.assertEqual(checks.hom(0, 1, 1), 0)
        self.assertEqual(checks.hom(0, 1, 0), 0.5)
        self.assertEqual(checks.hom(100, 1, 1), 0.5)
        self.assertAlmostEqual(checks.hom(-1, 2, 0.9), checks.hom(1, 2, 0.9))


class MatterTests(unittest.TestCase):
    def test_ice_official_verification_points(self):
        for temperature, pressure, rho, expansion in [
            (273.16, 611.657, 916.709492200, 1.59863102566e-4),
            (273.152519, 101325, 916.721463419, 1.59841589458e-4),
            (100, 1e8, 941.678203297, 2.58495528207e-5),
        ]:
            with self.subTest(temperature=temperature):
                got = checks.ice(temperature, pressure)
                self.assertAlmostEqual(got['rho_kg_m3'], rho, places=8)
                self.assertAlmostEqual(got['alpha_volume_K_inverse'], expansion, places=14)

    def test_ice_derivative_and_linear_volume_conversion(self):
        temperature, pressure, step = 250, 101325, 0.01
        v_plus = 1 / checks.ice(temperature + step, pressure)['rho_kg_m3']
        v_minus = 1 / checks.ice(temperature - step, pressure)['rho_kg_m3']
        mid = checks.ice(temperature, pressure)
        numerical = (v_plus - v_minus) / (2 * step) * mid['rho_kg_m3']
        self.assertAlmostEqual(numerical, mid['alpha_volume_K_inverse'], places=12)
        self.assertEqual(mid['alpha_volume_K_inverse'], 3 * mid['alpha_linear_average_K_inverse'])

    def test_helium_ideal_gas_is_not_lambda_prediction(self):
        result = checks.helium_check()
        self.assertGreater(result['ideal_bec_K'], result['measured_lambda_K'])
        self.assertGreater(result['relative_error'], 0.3)
        self.assertFalse(result['holographic_prediction'])

    def test_enthalpy_pressure_conversion(self):
        # One GPa times one cubic angstrom is 0.0062415 eV.
        self.assertAlmostEqual(checks.pressure_work_eV(1, 1), 0.006241509074460763)

    def test_london_profile(self):
        self.assertEqual(checks.london_field(0, 40e-9), 1)
        self.assertAlmostEqual(checks.london_field(40e-9, 40e-9), 1 / math.e)

    def test_full_radius_table_and_screening(self):
        rows = checks.atomic_check()['radii']
        self.assertEqual([r['Z'] for r in rows], list(range(1, 119)))
        self.assertEqual(rows[0]['reference_radius_pm'], 32)
        self.assertEqual(rows[-1]['symbol'], 'Og')
        self.assertIsNone(rows[-1]['ionization_eV'])
        self.assertAlmostEqual(rows[10]['slater_zeff'], 2.2)
        self.assertAlmostEqual(rows[1]['slater_zeff'], 1.7)
        self.assertNotEqual(rows[0]['unscreened_radius_pm'], rows[0]['reference_radius_pm'])

    def test_hydrogen_reduced_mass(self):
        result = checks.atomic_check()
        self.assertGreater(result['R_infinity_m_inverse'], result['R_hydrogen_m_inverse'])
        self.assertAlmostEqual(result['hydrogen_ionization_nonrel_eV'], 13.5982873, places=6)
        self.assertAlmostEqual(result['balmer_alpha_vacuum_nm'], 656.4696, places=3)

    def test_anomaly_is_g_over_two_minus_one(self):
        result = checks.constants_check()
        self.assertAlmostEqual((result['g_abs_reference'] - 2) / 2,
                               result['anomaly_reference'], places=15)
        self.assertGreater(abs(result['g_abs_one_loop'] - result['g_abs_reference']), 1e-6)

    def test_beta_continuum_endpoints(self):
        self.assertEqual(checks.beta_shape(0, 0.782), 0)
        self.assertEqual(checks.beta_shape(0.782, 0.782), 0)
        self.assertGreater(checks.beta_shape(0.3, 0.782), 0)
        result = checks.decay_check()
        self.assertAlmostEqual(result['neutron_Q_MeV'], 0.78233356, places=7)

    def test_liesegang_is_parameter_dependent(self):
        first = checks.bands(1e-9, 0.1, 1e-3, 8)
        faster = checks.bands(4e-9, 0.1, 1e-3, 8)
        for a, b in zip(first, faster):
            self.assertEqual(a['position_m'], b['position_m'])
            self.assertAlmostEqual(a['time_s'] / b['time_s'], 4)


class GravityAndMethodTests(unittest.TestCase):
    def test_archived_source_integrity(self):
        expected = {
            'data/codata-2022.txt': '77fb90e66c40db3e6eb16630bc9c88e4c7c8beddbe5e71be406f2f26e3f67e67',
            'data/covalent-radii.csv': '973baf5e9cdd6d18e551758464a1bbcdb12facfed85a03f89fc6d53cc6381081',
            'data/firas_monopole_spec_v1.txt': 'df793c3dca09ebfa7dbc5aa0ec1951daa8884431bc30eff28a710d7516cf50fa',
            'data/nist-ionization-raw.csv': '716e2a37da3316f6f544ace8607ab26b1e71b6e62ab0c132c005f87e231665a4',
            '../Голография/Рефлексия.txt': 'e9f95dd81b60d8ddd9de910c161f4beefdf8f08b148389bd10e99f53710c3eff',
        }
        for filename, digest in expected.items():
            with self.subTest(filename=filename):
                self.assertEqual(hashlib.sha256((checks.ROOT / filename).read_bytes()).hexdigest(), digest)

    def test_clock_limits_and_stable_small_speed(self):
        self.assertEqual(checks.clock_rate(0), 1)
        self.assertAlmostEqual(checks.clock_rate(0.6 * checks.C), 0.8)
        self.assertGreater(checks.clock_lag(1, 1), 0)
        self.assertLess(checks.clock_rate(0.99 * checks.C), 1)

    def test_gravity_scales_and_gps_signs(self):
        result = checks.gravity_check()
        earth, sun = result['bodies']
        self.assertAlmostEqual(earth['surface_g_m_s2'], 9.820, places=2)
        self.assertAlmostEqual(sun['schwarzschild_radius_m'], 2953.25, places=1)
        self.assertGreater(result['gps_net_us_day'], 38)
        self.assertLess(result['gps_kinetic_us_day'], 0)

    def test_horizon_area_increases_in_example(self):
        result = checks.horizon_check()
        self.assertGreater(result['area_ratio_nonrotating'], 1)
        self.assertNotAlmostEqual(result['euclidean_volume_proxy_ratio'], 1)
        self.assertEqual(checks.kerr_area_factor(0), 1)
        self.assertEqual(checks.kerr_area_factor(1), 0.5)

    def test_identifiability_rescaling(self):
        result = checks.operator_check()
        self.assertEqual(len(result['admissible_examples']), 3)
        self.assertEqual(len({x['signal_speed_m_s'] for x in result['admissible_examples']}), 1)
        self.assertEqual(len({x['illustrative_alpha'] for x in result['admissible_examples']}), 3)

    def test_invalid_physical_inputs(self):
        for call in [lambda: checks.planck(-1, 1), lambda: checks.planck(1, 0),
                     lambda: checks.clock_rate(checks.C), lambda: checks.clock_rate(-1),
                     lambda: checks.hom(0, 0, 1), lambda: checks.hom(0, 1, 2),
                     lambda: checks.ice(300, 101325), lambda: checks.kerr_area_factor(1.1),
                     lambda: checks.london_field(1, 0), lambda: checks.g2([1])]:
            with self.subTest(call=call), self.assertRaises(ValueError):
                call()

    def test_snapshot_and_external_working_directory(self):
        script = Path(checks.__file__).resolve()
        with tempfile.TemporaryDirectory() as directory:
            output = subprocess.run([sys.executable, str(script), '--check'], cwd=directory,
                                    text=True, capture_output=True, timeout=30)
        self.assertEqual(output.returncode, 0, output.stdout + output.stderr)
        saved = json.loads((script.parent / 'results/results.json').read_text())
        self.assertEqual(set(saved), set(checks.run_all()))


if __name__ == '__main__':
    unittest.main()
