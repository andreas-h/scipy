# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

__all__ = ['InverseDistanceInterpolator']

import numpy as np
from scipy.spatial import cKDTree as KDTree


class InverseDistanceInterpolator:

    """
    Inverse-distance-weighted interpolation

    .. versionadded:: 0.15

    Parameters
    ----------
    X : array_like, shape (n, m)
        The n points in m dimensions at which the data to be
        interpolated is defined.  This array is not copied unless this
        is necessary to produce a contiguous array of doubles, and so
        modifying this data will result in bogus results.

    z : array_like, shape (n, p)
        The values, defined in p dimensions, at the n points.

    leafsize : integer
        The number of points at which the KDtree algorithm switches
        over to brute-force.

    stat : bool
        accumulate wsum, wn for average weights

    scale : bool
        perform scaling before interpolation?

    Methods
    -------
    __call__

    Notes
    -----

    How many nearest neighbors should one take?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - start with 8, 11, 14, ..., 28 in 2d, 3d, 4d, ..., 10d; see
      Wendel's formula
    - make 3 runs with nnear, e.g., 6, 8, 10, and look at the results
      -- |interpol 6 - interpol 8| etc., or |f - interpol*| if you
      have f(q).

    Which value should one choose for p?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - p=2 weights nearer points more, farther points less.
    - In 2d, the circles around query points have areas ~ distance**2,
      so p=2 is inverse-area weighting. For example,

            (z1/area1 + z2/area2 + z3/area3)
            / (1/area1 + 1/area2 + 1/area3)
            = .74 z1 + .18 z2 + .08 z3  for distances 1, 2, 3.

    - Similarly, in 3d, p=3 is inverse-volume weighting.

    Scaling
    ~~~~~~~
    If different X coordinates measure different things, Euclidean
    distance can be way off.  For example, if X0 is in the range 0 to
    1 but X1 0 to 1000, the X1 distances will swamp X0; rescale the
    data, i.e. make X0.std() ~= X1.std().

    A nice property of IDW is that it's scale-free around query
    points: If I have values z1, z2, z3 from 3 points at distances d1,
    d2, d3, the IDW average

        (z1/d1 + z2/d2 + z3/d3) / (1/d1 + 1/d2 + 1/d3)

    is the same for distances 1, 2, 3, or 10, 20, 30 -- only the
    ratios matter. In contrast, the commonly-used Gaussian kernel
    exp(-(distance/h)**2) is exceedingly sensitive to distance and to
    h.

    See also
    --------
    Rbf : Radial basis functions for interpolation/smoothing scattered Nd data

    Example
    -------
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    Todo
    ----
    - anykernel( dj / av dj ) is also scale-free
    - error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    """

    # originally by @denis-bz, see http://stackoverflow.com/a/3119544/152439

    def __init__(self, X, z, leafsize=10, stat=False, scale=False):
        X, z = np.asarray(X), np.asarray(z)
        if X.ndim == 1:
            X, z = np.atleast_2d(X), np.atleast_2d(z)
        if not X.shape[0] == z.shape[0]:
            raise ValueError('len(X) %d != len(z) %d' % (len(X), len(z)))
        if not leafsize > 0:
            raise ValueError("leafsize has to be a positive integer")
        if scale:
            raise NotImplementedError()
        self.tree = KDTree(X, leafsize=leafsize)
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None

    def __call__(self, xi, nnear=6, eps=0.0, p=1, weights=None):
        """
        Interpolation at coordinates

        Parameters
        ----------
        xi : array_like, of shape (q, m)
            An array of points to query.  q may be one point, or a batch of points.

        nnear : int
            number of nearest neighbors to consider for each sample point

        eps : float
            eps: approximate nearest, dist <= (1 + eps) * true nearest

        p : float
            p: use 1 / distance**p

        weights :
            weights: optional multipliers for 1 / distance**p, of the same shape as q

        """

        if weights is not None and not np.all(weights >= 0.):
            raise ValueError("weights array has to be non-negative")

        xi = np.asarray(xi)
        qdim = xi.ndim
        if qdim == 1:
            xi = np.array([xi])
        #xi = np.atleast_2d(xi)

        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        self.distances, ix = self.tree.query(xi, k=nnear, eps=eps)
        interpol = np.zeros((len(self.distances), ) + np.shape(self.z[0]))
        for j, (dist, i) in enumerate(zip(self.distances, ix)):
            if nnear == 1:  # if nnear == 1, this is nearest-neighbor interp.
                wz = self.z[i]
            elif dist[0] < 1e-10:  # if the closest sample point is AT a node
                wz = self.z[i[0]]
            else:
                # weight z s by 1/dist
                w = 1 / dist ** p
                if weights is not None:
                    w *= weights[i]
                w /= np.sum(w)
                wz = np.dot(w, self.z[i])
                if self.stat:
                    self.wn += 1
                    self.wsum += w
            interpol[j] = wz
        return interpol if qdim > 1 else interpol[0]


if __name__ == '__main__':
    import sys

    N = 10000
    Ndim = 2
    Nask = N  # N Nask 1e5: 24 sec 2d, 27 sec 3d on mac g4 ppc
    Nnear = 8  # 8 2d, 11 3d => 5 % chance one-sided -- Wendel, mathoverflow.com
    leafsize = 10
    eps = .1  # approximate nearest, dist <= (1 + eps) * true nearest
    p = 1  # weights ~ 1 / distance**p
    cycle = .25
    seed = 1

    np.random.seed(seed)
    np.set_printoptions(3, threshold=100, suppress=True)  # .3f

    print('\nInverseDistanceInterpolator:  N %d  Ndim %d  Nask %d  Nnear %d  leafsize %d  eps %.2g  p %.2g' \
        % (N, Ndim, Nask, Nnear, leafsize, eps, p,))


    def terrain(x):
        """ ~ rolling hills """

        return np.sin(2 * np.pi / cycle * np.mean(x, axis=-1))


    known = np.random.uniform(size=(N, Ndim)) ** .5  # 1/(p+1): density x^p
    z = terrain(known)
    ask = np.random.uniform(size=(Nask, Ndim))

    invdisttree = InverseDistanceInterpolator(known, z, leafsize=leafsize, stat=True)
    interpol = invdisttree(ask, nnear=Nnear, eps=eps, p=p)

    print('average distances to nearest points: %s' \
        % np.mean(invdisttree.distances, axis=0))
    print('average weights: %s' % (invdisttree.wsum / invdisttree.wn))

        # see Wikipedia Zipf's law

    err = np.abs(terrain(ask) - interpol)
    print('average |terrain() - interpolated|: %.2g' % np.mean(err))

    # print "interpolate a single point: %.2g" % \
    #     invdisttree( known[0], nnear=Nnear, eps=eps )
