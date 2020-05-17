import unittest
import numpy as np
# import the functions


class TestGeneral(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        np.random.RandomState(0)

    def test_init(self):
        # tests should be written here
        self.assertAlmostEqual(1., 1.)
    

if __name__ == '__main__':
    unittest.main()
