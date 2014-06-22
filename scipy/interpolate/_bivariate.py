from __future__ import division, print_function, absolute_import

import numpy as np

try:
    from ._akima2d import mod_akima_2d
except ValueError:
    from _akima2d import mod_akima_2d


__all__ = ["Akima2DInterpolator"]


class Akima2DInterpolator(object):
    """Bivariate cubic interpolation based on local procedures, after H. Akima


    Parameters
    ----------
    x, y : array_like
        1-D arrays of shape (m, ), (n, ), respectively, or 2-D arrays
        both of shape (m, n), indicating the points of the data to be
        interpolated in x- and y- direction.

    z : array_like
        Array of values to be interpolated, of shape (m, n, ...)

    bounds_error : bool

    fill_value : float

    periodic_x : float

    periodic_y : float

    Methods
    -------
    __call__

    Notes
    -----
    .. versionadded:: 0.15

    Use only for precise data, as the fitted curve passes through the
    given points exactly.

    References
    ----------
    [1] A Method of Bivariate Interpolation and Smooth Surface Fitting Based
        on Local Procedures. Hiroshi Akima, Commun. ACM, 1974, 17(1), 18-20.
    [2] http://sosie.sourceforge.net/

    """

    def __init__(self, x, y, z, bounds_error=True, fill_value=None,
                 periodic_x=False, periodic_y=False):
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        if x.ndim != y.ndim:
            raise ValueError("x and y must have the same number of "
                             "dimensions.")
        if x.ndim > 2:
            raise ValueError("x and y must have either one or two dimensions.")
        if z.ndim < 2:
            raise ValueError("z must have at least two dimensions.")
        if x.ndim == 1 and z.shape[:2] != [x.shape[0], y.shape[0]]:
            raise ValueError("The first two dimensions of z must have the "
                             "same number of elements as x and y, "
                             "respectively")
        if x.ndim == 2 and not x.shape == y.shape == z.shape[:2]:
            raise ValueError("x and y must have the same shape as the first "
                             "two dimensions of z.")


if __name__ == "__main__":
    x, y = np.arange(3), np.arange(4)
    y, x = np.meshgrid(y, x)
    z = np.ones_like(x, dtype=float)

    poly, lon_in, lat_in = mod_akima_2d.akima_2d_init(False, x, y, z)

    xi, yi = [2.5, 5.5, 6.0], [6.7, 8.9]
    yi, xi = np.meshgrid(yi, xi)
    xi, yi = xi/10., yi/10.

    #            !integer intent(hide) :: nsys = 16
    #            !integer intent(hide),depend(poly) :: nx1 = size(poly, 1) - 3
    #            !integer intent(hide),depend(poly) :: ny1 = size(poly, 2) - 3
    #            !integer intent(hide),depend(z2) :: nx2 = size(x2, 1)
    #            !integer intent(hide),depend(z2) :: ny2 = size(y2, 2)

    res = mod_akima_2d.akima_2d_eval(poly, xi, yi, lon_in, lat_in,
                                    16, poly.shape[0] - 3, poly.shape[1] - 3, xi.shape[0], xi.shape[1])



    #res = mod_akima_2d.akima_2d_eval(poly, xi, yi, lon_in, lat_in)
