import unittest

import rasotools
import numpy as np
import xarray as xr

x = np.arange(-179.75, 180, 0.5)
y = np.arange(-89.75, 90, 0.5)
yy, xx = np.meshgrid(y, x)
phi = -25 * 9.81 * np.sin(2 * np.pi * xx / 360.) * np.sin(2 * np.pi * yy / 180.)
# print(phi.shape, x.shape, y.shape)
griddata = xr.DataArray(np.round(phi, 2), coords=[x, y], dims=('lon', 'lat'))


class Grid(unittest.TestCase):

    def setUp(self):
        pass

    def test_numpy(self):
        with self.assertRaises(ValueError):
            rasotools.grid.extract_locations(np.zeros(10), 45, 56)

    def test_lon(self):
        with self.assertRaises(ValueError):
            rasotools.grid.extract_locations(griddata.rename({'lon': 'x'}), 45, 56)

    def test_lat(self):
        with self.assertRaises(ValueError):
            rasotools.grid.extract_locations(griddata.rename({'lat': 'y'}), 45, 56)

    def test_method(self):
        with self.assertRaises(ValueError):
            rasotools.grid.extract_locations(griddata, 45, 56, method=None)


class Extraction(unittest.TestCase):

    def test_gridpoint(self):
        k = 0
        for i, j in zip(np.random.choice(x, 100), np.random.choice(y, 100)):
            with self.subTest(k=k):
                print(i,j)
                res = rasotools.grid.extract_locations(griddata, i, j, method='point', raw=True)
                print(res, griddata.sel(lon=i, lat=j))
                self.assertEqual(float(res), float(griddata.sel(lon=i, lat=j)))
                k += 1

    def test_gridnearest(self):
        pass

    def test_gridbilinear(self):
        pass

    def test_gridborders(self):
        pass

    def test_griddist(self):
        pass


if __name__ == '__main__':
    unittest.main()
