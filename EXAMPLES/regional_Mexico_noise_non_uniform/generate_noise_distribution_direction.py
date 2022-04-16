#!/usr/bin/env python3
#
# Generates noise distribution and direction binary files
# according to parfile_noise.yaml.
#
# The output binary files contain float values for all surface
# GLL points. They must be stored in DATABASES_MPI.
#
# IMPORTANT: This script reads the coordinates of all surface GLL
# points from "mask_noise.bin", which is stored in OUTPUT_FILES
# when running a noise simulation. Therefore, to use this script,
# you first need to generate "mask_noise.bin" by running a noise
# simulation using the default noise distribution.
#


import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import os
import yaml
from obspy.geodetics import gps2dist_azimuth
from geo import is_land


# CONSTANTS, DO NOT MODIFY
DTYPE = 'float32'

# from specfem3d_globe constants.h
CUSTOM_REAL = 4
PI = 3.141592653589793
PI_OVER_TWO = PI / 2.0
RADIANS_TO_DEGREES = 180.0 / PI
EARTH_FLATTENING_F = 1.0 / 299.8
EARTH_ONE_MINUS_F_SQUARED = (1.0 - EARTH_FLATTENING_F) ** 2
FACTOR_TAN = 1.0 / EARTH_ONE_MINUS_F_SQUARED


def gaussian_distribution(distribution, mask_noise, lon, lat):
    center_lon = distribution['center_lon']
    center_lat = distribution['center_lat']
    dist = np.zeros((lon.size))

    for i in range(0, dist.size):
        # great circle distance in meters using WGS84 ellipsoid
        dist[i] = gps2dist_azimuth(center_lat, center_lon, lat[i], lon[i])[0]
        dist[i] = np.abs(dist[i])

    gaussian = np.exp(- (dist ** 2) / (2 * distribution['sigma_m'] ** 2))
    mask_noise[:] += gaussian * distribution['weight']
    return mask_noise


def uniform_distribution(distribution, mask_noise):
    mask_noise[:] += distribution['weight']
    return mask_noise


def ocean_distribution(distribution, mask_noise, lon, lat):
    land = is_land(lon, lat, '50m')
    ocean_idx = np.where(land == 0)[0]

    if not ocean_idx.any():
        print('Skipping ocean distribution: There are no gll points on\
               the ocean.')
    else:
        mask_noise[ocean_idx] += distribution['weight']
    return mask_noise


def geoc_colat_to_lat(geoc_colat, perfect_sphere):
    '''
    from specfem3d_globe rthetaphi_xyz.f90/xyz_2_rlatlon_dble()

    converts geocentric colatitude in radians to
    geographic latitude in degrees
    '''

    # convert geocentric colatitude to geographic colatitude
    if perfect_sphere:
        colat = geoc_colat
    else:
        colat = PI_OVER_TWO -\
                np.arctan(FACTOR_TAN * np.cos(geoc_colat) / np.sin(geoc_colat))

    # convert geographic colatitude to latitude
    lat = PI_OVER_TWO - colat

    # convert to degrees
    lat *= RADIANS_TO_DEGREES
    return lat


def radians_to_degrees(angle):
    angle *= RADIANS_TO_DEGREES
    return angle


def read_mask_noise(idir, proc_ngll_surface, nproctot_val):
    '''
    mask_noise.bin includes 6 arrays:
     val_x_all  (r coordinate normalized over Earth radius)
     val_y_all  (theta coordinate radians, geocentric colatitude)
     val_z_all  (phi coordinate radians, longitude)
     val_ux_all (noise mask)
     val_uy_all (noise mask)
     val_uz_all (noise mask)

    dimension is (num_noise_surface_points, NPROCTOT_VAL)
    data type is real(kind=CUSTOM_REAL)

    val_ux_all, val_uy_all, val_uz_all have the same values
    according to noise_tomography.f90/save_mask()
    '''

    ifile = os.path.join(idir, 'mask_noise.bin')
    num_noise_surface_points = proc_ngll_surface
    ngll_surface = num_noise_surface_points * nproctot_val

    offset = CUSTOM_REAL * ngll_surface + 4
    fline = 0

    with open(ifile, 'rb') as _file:
        # read r coordinates (normalized)
        fline += 4
        _file.seek(fline)
        val_x = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_x = np.reshape(val_x, (num_noise_surface_points, nproctot_val),
                           order='F')
        fline += offset

        # read theta coordinates (geocentric colatitude)
        fline += 4
        _file.seek(fline)
        val_y = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_y = np.reshape(val_y, (num_noise_surface_points, nproctot_val),
                           order='F')
        fline += offset

        # read phi coordinates (geographic longitude)
        fline += 4
        _file.seek(fline)
        val_z = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_z = np.reshape(val_z, (num_noise_surface_points, nproctot_val),
                           order='F')
        fline += offset

        # read mask coordinates (x direction)
        fline += 4
        _file.seek(fline)
        val_ux = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_ux = np.reshape(val_ux, (num_noise_surface_points, nproctot_val),
                            order='F')
        fline += offset

        # read mask coordinates (y direction)
        fline += 4
        _file.seek(fline)
        val_uy = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_uy = np.reshape(val_uy, (num_noise_surface_points, nproctot_val),
                            order='F')
        fline += offset

        # read mask coordinates (z direction)
        fline += 4
        _file.seek(fline)
        val_uz = np.fromfile(_file, DTYPE, count=ngll_surface)
        val_uz = np.reshape(val_uz, (num_noise_surface_points, nproctot_val),
                            order='F')
        fline += offset

        return val_x, val_y, val_z, val_ux, val_uy, val_uz


def write_files(proc, mask_noise, normal_x_noise, normal_y_noise,
                normal_z_noise, outdir):
    _write(mask_noise,
           os.path.join(outdir, 'proc{:06}_reg1_mask_noise.bin'.format(proc)))
    _write(normal_x_noise,
           os.path.join(outdir, 'proc{:06}_reg1_normal_x_noise.bin'.format(proc)))
    _write(normal_y_noise,
           os.path.join(outdir, 'proc{:06}_reg1_normal_y_noise.bin'.format(proc)))
    _write(normal_z_noise,
           os.path.join(outdir, 'proc{:06}_reg1_normal_z_noise.bin'.format(proc)))
    return


def _write(x, filename):
    n = np.array([4 * len(x)], dtype='int32')
    x = np.array(x, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        x.tofile(file)
        n.tofile(file)


if __name__ == '__main__':
    try:
        with open('parfile_noise.yaml', 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: parfile_noise.yaml not found.')

    if par['PLOT_MASK']:
        all_gll_lon = np.array([])
        all_gll_lat = np.array([])
        all_mask_noise = np.array([])

    proc_ngll_surface = par['NSPEC2D_TOP_CM'] * par['NGLLX'] * par['NGLLY']

    # read mask_noise.bin
    gll_r, theta, phi, _, _, _ = read_mask_noise(par['INDIR'],
                                                 proc_ngll_surface,
                                                 par['NPROCTOT_VAL'])

    # specfem3d_globe uses geocentric colatitude in radians, convert to
    # geographic latitude in degrees
    gll_lat = geoc_colat_to_lat(theta, perfect_sphere=False)

    # specfem3d_globe uses geographic longitude in radians, convert to degrees
    gll_lon = radians_to_degrees(phi)

    for p in range(0, par['NPROCTOT_VAL']):
        print('Generating files for process {}'.format(p))

        mask_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_x_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_y_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_z_noise = np.zeros(proc_ngll_surface, dtype='float32')

        # set noise direction
        normal_x_noise[:] = par['XDIR']
        normal_y_noise[:] = par['YDIR']
        normal_z_noise[:] = par['ZDIR']

        # set noise distribution
        for d in par['DISTRIBUTIONS']:
            if d['type'] == 'uniform':
                mask_noise = uniform_distribution(d, mask_noise)
            elif d['type'] == 'ocean':
                mask_noise = ocean_distribution(d,
                                                mask_noise,
                                                gll_lon[:, p],
                                                gll_lat[:, p])
            elif d['type'] == 'gaussian':
                mask_noise = gaussian_distribution(d,
                                                   mask_noise,
                                                   gll_lon[:, p],
                                                   gll_lat[:, p])
            else:
                print('Undefined noise distribution.')

        if par['WRITE_FILES']:
            write_files(p, mask_noise, normal_x_noise, normal_y_noise,
                        normal_z_noise, par['OUTDIR'])

        if par['PLOT_MASK']:
            all_gll_lon = np.append(all_gll_lon, gll_lon[:, p])
            all_gll_lat = np.append(all_gll_lat, gll_lat[:, p])
            all_mask_noise = np.append(all_mask_noise, mask_noise)

    if par['PLOT_MASK']:
        mercator = ccrs.PlateCarree()

        ax = plt.axes(projection=mercator)
        ax.coastlines()
        im = ax.scatter(all_gll_lon, all_gll_lat, c=all_mask_noise,
                        transform=mercator)

        ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                     y_inline=False)
        ax.set_title('Noise distribution mask')

        cbar = plt.colorbar(im)
        cbar.set_label('Weight')

        plt.show()
        plt.close()
