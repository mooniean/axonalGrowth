import unittest
import numpy as np
from pyaxon.functions import *


class TestGeneral(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        np.random.RandomState(0)

    def test_init(self):
        # tests should be written here
        self.assertEqual(h(1.),h(1.))
    

if __name__ == '__main__':
    unittest.main()
