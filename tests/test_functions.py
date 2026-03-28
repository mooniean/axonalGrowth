"""
tests/general.py – unit tests for pyaxon/functions.py.

Each public function has its own TestCase class so failures are easy to locate.
A module-level RNG with a fixed seed is used wherever randomness is required
(NPY002 – modern Generator API only).
"""

import unittest

import numpy as np
from scipy.ndimage import convolve

from pyaxon.functions import (
    build_transport_field,
    coeff_beta_,
    com,
    conservative_upwind_advection,
    container,
    h,
    heaviside,
    init_cell,
    init_growth_factor,
    interface_pos,
    rasterized_line,
    smoothing,
    unit_vector,
)

# Fixed-seed RNG shared across all tests that need random numbers.
RNG = np.random.default_rng(0)

# Minimal 5×5 4th-order isotropic Laplacian stencil (same as parameters.py).
_dL = 1.0
STENCIL = (1.0 / (12.0 * _dL**2)) * np.array(
    [
        [0, 0, -1, 0, 0],
        [0, 0, 16, 0, 0],
        [-1, 16, -60, 16, -1],
        [0, 0, 16, 0, 0],
        [0, 0, -1, 0, 0],
    ]
)


# =============================================================================
# heaviside
# =============================================================================
class TestHeaviside(unittest.TestCase):
    def test_positive_input_returns_one(self):
        self.assertEqual(heaviside(1.0), 1)

    def test_zero_input_returns_one(self):
        """H(0) = 1 by convention."""
        self.assertEqual(heaviside(0.0), 1)

    def test_negative_input_returns_zero(self):
        self.assertEqual(heaviside(-0.001), 0)

    def test_vectorised_on_array(self):
        x = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
        expected = np.array([0, 0, 1, 1, 1])
        np.testing.assert_array_equal(heaviside(x), expected)


# =============================================================================
# h  (smoothed-step interpolant)
# =============================================================================
class TestH(unittest.TestCase):
    def test_h_at_zero(self):
        self.assertAlmostEqual(h(0.0), 0.0)

    def test_h_at_one(self):
        self.assertAlmostEqual(h(1.0), 1.0)

    def test_h_at_half(self):
        self.assertAlmostEqual(h(0.5), 0.5)

    def test_h_monotone_on_unit_interval(self):
        x = np.linspace(0.0, 1.0, 200)
        diffs = np.diff(h(x))
        self.assertTrue(np.all(diffs >= 0), "h should be monotonically non-decreasing on [0, 1]")

    def test_h_on_array(self):
        x = np.array([0.0, 0.5, 1.0])
        result = h(x)
        np.testing.assert_allclose(result, [0.0, 0.5, 1.0], atol=1e-12)


# =============================================================================
# smoothing
# =============================================================================
class TestSmoothing(unittest.TestCase):
    def _sharp_disc(self, shape=(30, 30), centre=(15, 15), radius=5):
        phi = np.zeros(shape)
        return init_cell(phi, centre, radius)

    def test_returns_array_same_shape(self):
        phi = self._sharp_disc()
        result = smoothing(phi.copy(), tstep=10, stencil=STENCIL)
        self.assertEqual(result.shape, phi.shape)

    def test_volume_is_conserved(self):
        phi = self._sharp_disc()
        V0 = h(phi).sum()
        result = smoothing(phi.copy(), tstep=50, stencil=STENCIL)
        V1 = h(result).sum()
        # Allow 0.5 % drift – 50 explicit steps on a coarse grid have a small
        # but nonzero Lagrange-multiplier residual that is physically acceptable.
        self.assertAlmostEqual(V0, V1, delta=V0 * 5e-3)

    def test_smoothing_reduces_gradient_magnitude(self):
        """After smoothing the interface should be more diffuse."""
        phi = self._sharp_disc()
        result = smoothing(phi.copy(), tstep=100, stencil=STENCIL)
        before = np.linalg.norm(np.gradient(phi))
        after = np.linalg.norm(np.gradient(result))
        # A diffuse interface has smaller peak gradients.
        self.assertLess(after, before)

    def test_values_stay_bounded(self):
        phi = self._sharp_disc()
        result = smoothing(phi.copy(), tstep=50, stencil=STENCIL)
        self.assertTrue(np.all(result >= -0.05))
        self.assertTrue(np.all(result <= 1.05))


# =============================================================================
# com  (centre of mass)
# =============================================================================
class TestCom(unittest.TestCase):
    def test_single_point_mass(self):
        phi = np.zeros((10, 10))
        phi[3, 7] = 1.0
        r = com(phi, None, ndim=2)
        np.testing.assert_allclose(r, [3.0, 7.0])

    def test_symmetric_disc_centred_on_grid(self):
        phi = np.zeros((20, 20))
        centre = np.array([10, 10])
        phi = init_cell(phi, centre, cell_radius=4)
        r = com(phi, None, ndim=2)
        np.testing.assert_allclose(r, centre, atol=0.5)

    def test_uniform_field_returns_geometric_centre(self):
        phi = np.ones((8, 12))
        r = com(phi, None, ndim=2)
        np.testing.assert_allclose(r, [3.5, 5.5], atol=1e-10)


# =============================================================================
# interface_pos  (1-D centre-of-mass)
# =============================================================================
class TestInterfacePos(unittest.TestCase):
    def test_single_spike(self):
        phi = np.zeros((1, 10))
        phi[0, 4] = 1.0
        self.assertAlmostEqual(interface_pos(phi), 4.0)

    def test_uniform_returns_middle(self):
        phi = np.ones((1, 10))
        # weighted mean of indices 0..9 under uniform weights = 4.5
        self.assertAlmostEqual(interface_pos(phi), 4.5)


# =============================================================================
# init_cell
# =============================================================================
class TestInitCell(unittest.TestCase):
    def test_centre_is_one(self):
        phi = np.zeros((20, 20))
        phi = init_cell(phi, cell_position=[10, 10], cell_radius=3)
        self.assertAlmostEqual(phi[10, 10], 1.0)

    def test_far_corner_is_zero(self):
        phi = np.zeros((20, 20))
        phi = init_cell(phi, cell_position=[10, 10], cell_radius=3)
        self.assertAlmostEqual(phi[0, 0], 0.0)

    def test_only_zeros_and_ones(self):
        phi = np.zeros((20, 20))
        phi = init_cell(phi, cell_position=[10, 10], cell_radius=5)
        unique = set(np.unique(phi))
        self.assertTrue(unique.issubset({0.0, 1.0}))

    def test_disc_radius_respected(self):
        """All active cells must lie within the specified radius."""
        shape = (30, 30)
        centre = np.array([15, 15])
        radius = 6
        phi = np.zeros(shape)
        phi = init_cell(phi, centre, radius)
        active = np.argwhere(phi == 1.0)
        distances = np.linalg.norm(active - centre, axis=1)
        self.assertTrue(np.all(distances <= radius + 1e-9))

    def test_returns_modified_array(self):
        phi = np.zeros((10, 10))
        result = init_cell(phi, [5, 5], 2)
        self.assertIs(result, phi)


# =============================================================================
# init_growth_factor
# =============================================================================
class TestInitGrowthFactor(unittest.TestCase):
    def setUp(self):
        self.shape = (20, 20)
        self.pos = np.array([10, 10])
        self.radius = 5

    def test_inside_exclusion_zone_stays_zero(self):
        ngf = np.zeros(self.shape)
        ngf = init_growth_factor(ngf, self.pos, self.radius, rng=RNG)
        for i in range(ngf.size):
            s = np.asarray(np.unravel_index(i, ngf.shape))
            d = np.linalg.norm(s - self.pos)
            if d <= self.radius + 2:
                self.assertAlmostEqual(ngf.flat[i], 0.0,
                    msg=f"Cell at {s} (d={d:.1f}) inside exclusion zone should be 0")

    def test_outside_exclusion_zone_is_nonnegative(self):
        ngf = np.zeros(self.shape)
        ngf = init_growth_factor(ngf, self.pos, self.radius, rng=RNG)
        self.assertTrue(np.all(ngf >= 0.0))

    def test_outside_exclusion_zone_bounded_by_ten(self):
        ngf = np.zeros(self.shape)
        ngf = init_growth_factor(ngf, self.pos, self.radius, rng=RNG)
        self.assertTrue(np.all(ngf <= 10.0))

    def test_reproducible_with_seeded_rng(self):
        ngf1 = init_growth_factor(np.zeros(self.shape), self.pos, self.radius,
                                   rng=np.random.default_rng(42))
        ngf2 = init_growth_factor(np.zeros(self.shape), self.pos, self.radius,
                                   rng=np.random.default_rng(42))
        np.testing.assert_array_equal(ngf1, ngf2)

    def test_different_seeds_give_different_fields(self):
        ngf1 = init_growth_factor(np.zeros(self.shape), self.pos, self.radius,
                                   rng=np.random.default_rng(0))
        ngf2 = init_growth_factor(np.zeros(self.shape), self.pos, self.radius,
                                   rng=np.random.default_rng(99))
        self.assertFalse(np.array_equal(ngf1, ngf2))


# =============================================================================
# container
# =============================================================================
class TestContainer(unittest.TestCase):
    def test_output_shape_matches_input(self):
        f1 = RNG.random((5, 5))
        f2 = RNG.random((5, 5))
        out = container(f1, f2, alpha=1e-4)
        self.assertEqual(out.shape, f1.shape)

    def test_zero_field_gives_zero(self):
        f1 = np.zeros((5, 5))
        f2 = RNG.random((5, 5))
        out = container(f1, f2, alpha=1e-4)
        np.testing.assert_allclose(out, 0.0, atol=1e-12)

    def test_inside_carrier_weaker_than_outside(self):
        """Confinement should be stronger where field_2 ≈ 0 (outside)."""
        f1 = np.ones((2, 2)) * 0.5
        alpha = 1e-4
        inside = container(f1, np.ones((2, 2)) * 0.99, alpha)   # field_2 ≈ 1
        outside = container(f1, np.zeros((2, 2)), alpha)         # field_2 ≈ 0
        # The absolute confinement value is larger outside the carrier.
        self.assertTrue(np.all(np.abs(outside) > np.abs(inside)))

    def test_finite_for_typical_values(self):
        f1 = RNG.random((10, 10))
        f2 = RNG.random((10, 10))
        out = container(f1, f2, alpha=1e-4)
        self.assertTrue(np.all(np.isfinite(out)))


# =============================================================================
# unit_vector
# =============================================================================
class TestUnitVector(unittest.TestCase):
    def test_normalises_nonzero_scalar(self):
        self.assertAlmostEqual(unit_vector(3.0, 3.0), 1.0)
        self.assertAlmostEqual(unit_vector(6.0, 2.0), 3.0)

    def test_returns_zero_for_zero_norm(self):
        self.assertEqual(unit_vector(5.0, 0.0), 0)

    def test_vectorised_over_array(self):
        v = np.array([3.0, -4.0, 0.0])
        n = np.array([3.0, 4.0, 1.0])
        result = unit_vector(v, n)
        np.testing.assert_allclose(result, [1.0, -1.0, 0.0])


# =============================================================================
# coeff_beta_
# =============================================================================
class TestCoeffBeta(unittest.TestCase):
    def test_above_threshold_returns_reciprocal(self):
        self.assertAlmostEqual(coeff_beta_(0.5), 1.0 / 0.5)

    def test_at_saturation_returns_one(self):
        self.assertAlmostEqual(coeff_beta_(1.0), 1.0)

    def test_below_threshold_returns_zero(self):
        self.assertAlmostEqual(coeff_beta_(0.0), 0.0)
        self.assertAlmostEqual(coeff_beta_(0.05), 0.0)

    def test_vectorised_on_array(self):
        mtb = np.array([0.0, 0.05, 0.5, 1.0])
        result = coeff_beta_(mtb)
        expected = np.array([0.0, 0.0, 2.0, 1.0])
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_monotone_above_threshold(self):
        """1/mtb is strictly decreasing for mtb > threshold."""
        mtb = np.linspace(0.2, 1.0, 50)
        result = coeff_beta_(mtb)
        diffs = np.diff(result)
        self.assertTrue(np.all(diffs <= 0))


# =============================================================================
# rasterized_line
# =============================================================================
class TestRasterizedLine(unittest.TestCase):
    def test_horizontal_line_keeps_endpoints(self):
        line = rasterized_line((5, 2), (5, 8))
        self.assertTupleEqual(tuple(line[0]), (5, 2))
        self.assertTupleEqual(tuple(line[-1]), (5, 8))

    def test_vertical_line_keeps_endpoints(self):
        line = rasterized_line((2, 5), (8, 5))
        self.assertTupleEqual(tuple(line[0]), (2, 5))
        self.assertTupleEqual(tuple(line[-1]), (8, 5))

    def test_diagonal_line_keeps_endpoints(self):
        line = rasterized_line((0, 0), (5, 5))
        self.assertTupleEqual(tuple(line[0]), (0, 0))
        self.assertTupleEqual(tuple(line[-1]), (5, 5))

    def test_no_duplicate_points(self):
        line = rasterized_line((20, 20), (20, 35))
        # All consecutive pairs must differ in at least one coordinate.
        diffs = np.diff(line, axis=0)
        self.assertTrue(np.all(np.any(diffs != 0, axis=1)))

    def test_single_point_returns_one_row(self):
        line = rasterized_line((4, 7), (4, 7))
        self.assertEqual(line.shape[0], 1)
        self.assertTupleEqual(tuple(line[0]), (4, 7))

    def test_monotone_column_for_horizontal_segment(self):
        line = rasterized_line((20, 20), (20, 35))
        self.assertTrue(np.all(np.diff(line[:, 1]) >= 0))

    def test_integer_dtype(self):
        line = rasterized_line((1.4, 2.6), (5.5, 7.5))
        self.assertTrue(np.issubdtype(line.dtype, np.integer))

    def test_float_inputs_rounded_correctly(self):
        line = rasterized_line((1.4, 2.6), (1.4, 5.4))
        self.assertTupleEqual(tuple(line[0]), (1, 3))
        self.assertTupleEqual(tuple(line[-1]), (1, 5))


# =============================================================================
# build_transport_field
# =============================================================================
class TestBuildTransportField(unittest.TestCase):
    def test_zero_field_for_single_point(self):
        track = [(5, 5)]
        v = build_transport_field(track, (10, 10), support_radius=0)
        np.testing.assert_array_equal(v, np.zeros((2, 10, 10)))

    def test_zero_field_for_empty_track(self):
        v = build_transport_field([], (10, 10), support_radius=0)
        np.testing.assert_array_equal(v, np.zeros((2, 10, 10)))

    def test_horizontal_track_points_in_column_direction(self):
        """A track along row 10, columns 2→6, should give v[1] > 0, v[0] ≈ 0."""
        track = rasterized_line((10, 2), (10, 6))
        v = build_transport_field(track, (20, 20), support_radius=0)
        self.assertGreater(v[1, 10, 3], 0.0)
        self.assertAlmostEqual(v[0, 10, 3], 0.0)

    def test_vertical_track_points_in_row_direction(self):
        """A track along col 10, rows 2→6, should give v[0] > 0, v[1] ≈ 0."""
        track = rasterized_line((2, 10), (6, 10))
        v = build_transport_field(track, (20, 20), support_radius=0)
        self.assertGreater(v[0, 3, 10], 0.0)
        self.assertAlmostEqual(v[1, 3, 10], 0.0)

    def test_velocity_magnitude_leq_one_everywhere(self):
        track = rasterized_line((5, 5), (5, 15))
        v = build_transport_field(track, (20, 20), support_radius=1)
        mag = np.linalg.norm(v, axis=0)
        self.assertTrue(np.all(mag <= 1.0 + 1e-9))

    def test_output_shape(self):
        track = rasterized_line((3, 3), (3, 8))
        v = build_transport_field(track, (15, 20), support_radius=0)
        self.assertEqual(v.shape, (2, 15, 20))

    def test_support_radius_spreads_field(self):
        """Larger support_radius should fill more non-zero cells."""
        track = rasterized_line((10, 2), (10, 8))
        v0 = build_transport_field(track, (20, 20), support_radius=0)
        v1 = build_transport_field(track, (20, 20), support_radius=2)
        n0 = np.count_nonzero(np.linalg.norm(v0, axis=0))
        n1 = np.count_nonzero(np.linalg.norm(v1, axis=0))
        self.assertGreater(n1, n0)


# =============================================================================
# conservative_upwind_advection
# =============================================================================
class TestConservativeUpwindAdvection(unittest.TestCase):
    def test_mass_conservation_row_flow(self):
        """Sum of advection term must be zero (conservation)."""
        field = np.zeros((6, 6))
        field[3, 2] = 1.0
        v = np.zeros((2, 6, 6))
        v[0] = 1.0   # flow in row direction
        adv = conservative_upwind_advection(field, v, dL=1.0)
        self.assertAlmostEqual(float(adv.sum()), 0.0, places=12)

    def test_mass_conservation_column_flow(self):
        field = np.zeros((6, 6))
        field[2, 3] = 1.0
        v = np.zeros((2, 6, 6))
        v[1] = 1.0
        adv = conservative_upwind_advection(field, v, dL=1.0)
        self.assertAlmostEqual(float(adv.sum()), 0.0, places=12)

    def test_upwind_sign_column_direction(self):
        """With positive column velocity, mass should flow left→right."""
        field = np.zeros((5, 5))
        field[2, 1] = 1.0
        v = np.zeros((2, 5, 5))
        v[1] = 1.0
        adv = conservative_upwind_advection(field, v, dL=1.0)
        self.assertLess(adv[2, 1], 0.0)    # source cell loses mass
        self.assertGreater(adv[2, 2], 0.0)  # downstream cell gains mass

    def test_upwind_sign_row_direction(self):
        """With positive row velocity, mass should flow upward in rows."""
        field = np.zeros((5, 5))
        field[1, 2] = 1.0
        v = np.zeros((2, 5, 5))
        v[0] = 1.0
        adv = conservative_upwind_advection(field, v, dL=1.0)
        self.assertLess(adv[1, 2], 0.0)
        self.assertGreater(adv[2, 2], 0.0)

    def test_zero_velocity_gives_zero_advection(self):
        field = RNG.random((8, 8))
        v = np.zeros((2, 8, 8))
        adv = conservative_upwind_advection(field, v, dL=1.0)
        np.testing.assert_allclose(adv, 0.0, atol=1e-15)

    def test_dL_scales_flux_correctly(self):
        """Doubling dL should halve the advection term."""
        field = np.zeros((5, 5))
        field[2, 2] = 1.0
        v = np.zeros((2, 5, 5))
        v[1] = 1.0
        adv1 = conservative_upwind_advection(field, v, dL=1.0)
        adv2 = conservative_upwind_advection(field, v, dL=2.0)
        np.testing.assert_allclose(adv2, adv1 * 0.5, atol=1e-15)

    def test_output_shape_matches_field(self):
        field = RNG.random((7, 9))
        v = np.zeros((2, 7, 9))
        adv = conservative_upwind_advection(field, v)
        self.assertEqual(adv.shape, field.shape)

    def test_finite_output_for_typical_values(self):
        field = RNG.random((10, 10))
        v = (RNG.random((2, 10, 10)) - 0.5) * 2.0   # in (-1, 1)
        adv = conservative_upwind_advection(field, v, dL=1.0)
        self.assertTrue(np.all(np.isfinite(adv)))


if __name__ == "__main__":
    unittest.main()
