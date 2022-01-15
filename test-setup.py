import unittest
import mainfile as m

class TestPotts4Wolff(unittest.TestCase):

    def test_H(self):
        result = m.H(m.grid)
        self.assertLessEqual(result, 0)
        self.assertGreaterEqual(result, -2 * m.L**2)

    
    def test_localH(self):
        result = m.localH(m.grid, m.i, m.j)
        self.assertLessEqual(result, 0)
        self.assertGreaterEqual(result, -4)


    def test_m(self):
        result = m.m(m.grid)
        self.assertLessEqual(result, 1)
        self.assertGreaterEqual(result, 0)        

    
    def test_p(self):
        self.assertGreaterEqual(m.Wolff.p, 0)
        self.assertLessEqual(m.Wolff.p, 1)


    def test_sites(self):
        for i, j in m.sites_in_cluster:
            self.assertGreaterEqual(m.i, 0)
            self.assertGreaterEqual(m.j, 0)
            self.assertLessEqual(m.i, 9)
            self.assertLessEqual(m.j, 9)


if __name__ == '__main__':
    unittest.main()