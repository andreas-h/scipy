from __future__ import division, print_function, absolute_import

from numpy.testing import (assert_, assert_equal, assert_almost_equal,
        assert_array_almost_equal, assert_raises, assert_array_equal,
        dec, TestCase, run_module_suite, assert_allclose)
import numpy as np

from scipy.interpolate import InverseDistanceInterpolator


class TestInverseDistanceInterpolator(TestCase):

    def test_1neighbor(self):
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x[::-1]  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = x[1:-1]
        actual = idw(x + .4, nnear=1)
        expected = y
        assert_array_equal(actual, expected)

    def test_1neighbor_1d(self):
        x = np.arange(10)
        y = x[::-1]
        idw = InverseDistanceInterpolator(x, y)
        xi = x[1:-1]
        actual = idw(x + .4, nnear=1)
        expected = y
        assert_array_equal(actual, expected)

    def test_query_at_knots(self):
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x[::-1]  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        actual = idw(x)
        expected = y
        assert_array_equal(actual, expected)

    def test_query_at_knots_1d(self):
        x = np.arange(10)
        y = x[::-1]
        idw = InverseDistanceInterpolator(x, y)
        actual = idw(x)
        expected = y
        assert_array_equal(actual, expected)

    def test_query_at_knots_2d(self):
        x = np.column_stack((np.arange(10), np.arange(10)))  # shape (10, 2)
        y = x[::-1, 0]  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        actual = idw(x)
        expected = y
        assert_array_equal(actual, expected)

    def test_1neighbor_2d(self):
        x = np.column_stack((np.arange(10), np.arange(10)))  # shape (10, 2)
        y = x[::-1, 0]  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        actual = idw(x - 0.4999, nnear=1)
        expected = y
        assert_array_equal(actual, expected)

    def test_query_between_knots_half(self):
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = (x[1:] + x[:-1]) / 2.
        actual = idw(xi, nnear=2)
        expected = (y[1:] + y[:-1]) / 2.
        assert_array_equal(actual, expected)

    def test_query_between_knots_third(self):
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = (x[1:] + 2. * x[:-1]) / 3.
        actual = idw(xi, nnear=2)
        expected = xi
        assert_array_almost_equal(actual, expected)

    def test_query_between_knots_quarter(self):
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = x[:-1] + .25
        actual = idw(xi, nnear=2)
        expected = xi
        assert_array_almost_equal(actual, expected)

    def test_query_between_knots_half_4n(self):
        print("test_query_between_knots_half_4n")
        x = np.atleast_2d(np.arange(10)).T  # shape (10, 1)
        y = x[::-1]  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = x[1:-2] + .25
        actual = idw(xi, nnear=4)
        #print()
        #print("x",x.ravel())
        #print("y",y.ravel())
        #print("xi",xi.ravel())
        #print("actual",actual.ravel())
        #print("expected",expected.ravel())
        expected = (y[1:] + y[:-1]) / 2.
        assert_array_equal(actual, expected)

    def test_1neighbor_area_2d(self):
        x = np.column_stack((np.arange(10), np.arange(10)))  # shape (10, 2)
        y = np.arange(20).reshape((10, 2))  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        actual = idw(x - 0.4999, nnear=1)
        expected = y
        assert_array_equal(actual, expected)

    def test_4neighbors_area_2d(self):
        x = np.column_stack((np.arange(10), np.arange(10)))  # shape (10, 2)
        y = np.arange(20).reshape((10, 2))  # shape (10, 1)
        idw = InverseDistanceInterpolator(x, y)
        xi = x[:-1] + .5
        actual = idw(xi, nnear=4, p=2.)
        expected = np.asarray([y[i:i+2].mean() for i in range(xi.shape[0])])
        #print()
        #print("test_4neighbors_area_2d")
        #print("x",x.T)
        #print("y",y.T)
        #print("xi",xi.T)
        #print("actual",actual.T)
        #print("expected",expected.T)
        assert_array_equal(actual, expected)


if __name__ == "__main__":
    run_module_suite()
