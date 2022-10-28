'''

This code has been borrowed from the open-source Python tool NOISI
under the LGPL-3.0 license.

For the original code, refer to:
    https://github.com/lermert/noisi/blob/master/noisi/util/geo.py

'''

import numpy as np
try:
    import cartopy.io.shapereader as shpreader
    import shapely.geometry as sgeom
    from shapely.ops import unary_union
    from shapely.prepared import prep
except ImportError:
    pass


def is_land(x, y, res="110m"):

    if 'prep' not in globals():
        raise ImportError("cartopy is needed to design ocean-only source.")
    assert(res in ["10m", "50m", "110m"]), "Resolution must be 10m, 50m, 110 m"

    land_shp_fname = shpreader.natural_earth(resolution=res,
                                             category='physical',
                                             name='land')

    land_geom = unary_union(list(shpreader.Reader(land_shp_fname).
                                 geometries()))
    land = prep(land_geom)
    is_land = np.zeros(len(x))
    for i in range(len(x)):
        is_land[i] = land.contains(sgeom.Point(x[i], y[i]))
    return is_land

