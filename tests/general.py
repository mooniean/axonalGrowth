import unittest

import numpy as np
from pyaxon.functions import build_transport_field, conservative_upwind_advection, rasterized_line


class TestGeneral(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Use the modern Generator API (NPY002) with a fixed seed for reproducibility.
        cls.rng = np.random.default_rng(0)

    def test_init(self):
        # tests should be written here
        self.assertAlmostEqual(1.0, 1.0)

    def test_rasterized_line_keeps_endpoints(self):
        line = rasterized_line((20, 20), (20, 35))
        self.assertTupleEqual(tuple(line[0]), (20, 20))
        self.assertTupleEqual(tuple(line[-1]), (20, 35))
        self.assertTrue(np.all(np.diff(line[:, 1]) >= 0))

    def test_build_transport_field_follows_track(self):
        track = rasterized_line((10, 10), (10, 14))
        velocity = build_transport_field(track, (20, 20), support_radius=0)
        self.assertGreater(velocity[1, 10, 10], 0.0)
        self.assertAlmostEqual(velocity[0, 10, 10], 0.0)

    def test_conservative_upwind_advection_is_mass_conservative(self):
        field = np.zeros((5, 5))
        field[2, 1] = 1.0
        velocity = np.zeros((2, 5, 5))
        velocity[1, :, :] = 1.0
        advection = conservative_upwind_advection(field, velocity, dL=1.0)

        self.assertAlmostEqual(float(advection.sum()), 0.0)
        self.assertLess(advection[2, 1], 0.0)
        self.assertGreater(advection[2, 2], 0.0)


if __name__ == "__main__":
    unittest.main()
