import math, unittest
import checks as c


class Horizon(unittest.TestCase):
    def test_sun_schwarzschild_radius(self):
        self.assertAlmostEqual(c.rs(c.M_SUN), 2953.25, delta=0.1)

    def test_mass_is_planck_force_times_half_radius(self):
        m = 10 * c.M_SUN
        self.assertAlmostEqual(m * c.C**2 / (c.F_P * c.rs(m) / 2), 1.0, places=12)

    def test_area_volume_ratio_not_invariant(self):
        r = c.rs(c.M_SUN)
        self.assertAlmostEqual(4 * math.pi * r**2 / (4 / 3 * math.pi * r**3), 3 / r)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            c.rs(0)
        with self.assertRaises(ValueError):
            c.kerr_area(1, 1.5)


class Merger(unittest.TestCase):
    def test_sphere_loss_is_21_percent(self):
        self.assertAlmostEqual(c.sphere_merge_area_loss(), 0.2063, places=4)

    def test_gr_area_grows(self):
        m = c.merger_check()
        self.assertTrue(m['hawking_area_theorem_holds'])
        self.assertGreater(m['gr_equal_mass_area_ratio'], 1.5)

    def test_extremal_kerr_half_area(self):
        self.assertAlmostEqual(c.kerr_area(1e30, 1) / c.kerr_area(1e30, 0), 0.5)


class Neck(unittest.TestCase):
    def test_cut_length(self):
        self.assertAlmostEqual(c.cut_length(1e-12, 1e-6), 1e-6)

    def test_neck_shrinks(self):
        self.assertLess(c.neck_radius(1e-9), c.neck_radius(1e-6))

    def test_levels_ordered(self):
        lv = c.cascade_check()['levels']
        self.assertLess(lv['atomic']['L_cut_m'], lv['molecular']['L_cut_m'])
        for v in lv.values():
            self.assertAlmostEqual(c.neck_radius(v['tau_match_s']), v['L_cut_m'])

    def test_chirality(self):
        r = c.chirality_check()
        self.assertTrue(r['sign_flips_with_gradient'])
        self.assertEqual(r['achiral_Lz'], 0.0)


if __name__ == '__main__':
    unittest.main()
