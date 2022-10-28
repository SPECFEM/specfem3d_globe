#!/usr/bin/env python
#
# script to create crustmap files for a generic 2/3-layer mars or moon crust
#
# required python modules:
#   - numpy
#   - pyshtools
#
# There are 4 types of crustmaps:
#    topography crustmaps (e.g., marscrustt3.cmap - marscrustt7.cmap)
#    s-velocity crustmaps (      marscrusts3.cmap - marscrusts7.cmap)
#    p-velocity crustmaps (      marscrustp3.cmap - marscrustp7.cmap)
#    density    crustmaps (      marscrustr3.cmap - marscrustr7.cmap)
#
# These are simple ascii data files and you can change them easily.
# For Mars crustmaps, the current crustmap resolution is set to 1/6 degree.
#
# Values are ordered on a grid that runs from 90 degrees lat (0 degrees lng to 180 degrees lng ) to -90 degrees lat.
# Every line in a crustmap files sets all values for one specific latitude:
#
#    (length of one line = 360 * 4 entries)
#    (total number of lines = 180 * 4)
#
# The value for (lat = 89.875, long = 0.125) is stored in entry (line 1,column 1).
# You get the value that is north west of your (lat, long) point by using:
#
#    line = int(0.5+((90.0-lat)*CRUSTMAP_RESOLUTION))
#    column = int0.5+(xlng*CRUSTMAP_RESOLUTION)
#
# Crustal model generation:
#
# 1. produce Moho depth maps:
#
#    examples are provide in folder moho/ for different Mars and Moon crusts
#
#    check out Mark's ctplanet tool for more:
#    see: https://github.com/MarkWieczorek/ctplanet
#
# 2. modify the parameters in this script:
#      - choose number of crustal layers
#      - choose crustmap resolution
#      - choose if a perturbation should be added to the upper crustal layer
#
# 3. run script, for example
#    > ./run_create_a_crustmap.py moho/Moho-Mars-DWAK-39-2550.sh
#
# this will create all crustmap files and put them into a temporary subfolder tmp/
#
# to run a SPECFEM simulation with these new crustmap model, you can copy them into DATA/crustmap/ folder
# and append **_crustmap to the MODEL name in Par_file (in case of Mars models, **_3D_crust will have the same effect)
#
#
from __future__ import print_function

import os
import sys

## numpy
try:
    import numpy as np
except ImportError:
    print('This script requires NumPy.')
    sys.exit()

## shtools
try:
    import pyshtools as pysh
except ImportError:
    print('This script requires pyshtools. Please install, e.g., by: pip install pyshtools')
    sys.exit()

# matplotlib
import matplotlib.pyplot as plt

# print-function always flush to buffer (requires python3)
import functools
print = functools.partial(print, flush=True)

####################################################################
# User parameters

# number of crustal layers ( 2 == upper/lower crust, 3 == upper/middle/lower crust)
num_layers_crust = 3

# output crustmaps resolution
CRUSTMAP_RESOLUTION = 6     # e.g. 6 means 1/6 degrees - default for mars crustmaps
                            # earth 1 degree ~ 111 km; mars 1 degree ~ 59 km; moon 1 degree ~ 30 km

# adds perturbations to velocity values
add_perturbations = False

# perturbation distribution (0 == power-law / 1 == von Karman)
perturbation_dist = 1

# perturbation strength
perturbation_max_amplitude = 0.1   # 0.1 means 10% maximum perturbation

# correlation wavelength (in km)
perturbation_corr_wavelength = 5.0     # 1 == 1km wavelength, 10 == 10km correlation length

# planet size for conversion (km -> degree)
# for Earth radius 6371 km, for Mars 3389.5 km, for Moon 1737.1 km
R_PLANET = 3389.5    # mars

# plotting: show figures
show_maps = True

####################################################################

def psd_powerlaw_1D(k):
    # power-law w/ degree -2
    psd = k**(-2)
    return psd

def psd_vonKarman_1D(a_in,H_in,sigma_in,k):
    # von Karman distribution
    #
    # scattering regimes: k * a >> 1 (lambda << a) - high-frequency scattering
    #                     k * a ~ 1  (lambda ~  a) - diffraction condition
    #                     k * a << 1 (lambad >> a) - low-frequency scattering

    # Frankel and Clayton (1986): 1-D Fourier transform:
    #    P = a / (1 + k**2 * a**2)**0.5
    # such that correlation distance k*a ~ 1 for characteristic signal
    #
    # 1D power spectral density from Frankel & Clayton:
    #a = a_in / (2.0 * np.pi)
    #psd = a / (1.0 + k**2 * a**2)**0.5
    #return psd

    # von Karman 1D
    # 1D power spectral density: Imperatori & Mai (2012)
    #
    # P(k) = sigma**2 ( 2 sqrt(pi * a) )**E gamma(H + E/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + E/2) )
    #      = sigma**2 ( 2 sqrt(pi * a) )**1 gamma(H + 1/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 1/2) )
    #      = sigma**2 * 2 * sqrt(pi * a) * gamma(H + 1/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 1/2) )
    #      =             amp             *        g1  / (   g2    * (1 + k**2 a**2 )**(H + 1/2)
    from scipy.special import gamma

    a = a_in          # correlation length
    H = H_in          # Hurst exponent ~ 0.1-0.3
    sigma = sigma_in  # standard deviation ~ 10%

    g1 = gamma(H + 0.5)                                   # gamma function: gamma(H+1/2)
    g2 = gamma(H)                                         # gamma function: gamma(H)
    amp = sigma**2 * ( 2 * np.sqrt(np.pi * a) )           # coefficient   : sigma**2 * ( 2 * sqrt(pi * a) )

    psd = amp * g1 / ( g2 * (1.0 + k**2 * a**2)**(H + 0.5) )
    return psd



def get_map_pertub_powerlaw(degrees,lmax_calc):
    """
    creates a map based on spherical harmonics with a power-law distribution
    """
    global perturbation_max_amplitude,show_maps

    # power-law distribution
    # power-law w/ degree -2
    psd_power = psd_powerlaw_1D(degrees)

    # psd plot
    if show_maps:
        plt.plot(degrees,psd_power,label="psd - power-law 1D"); plt.xscale('log');
        plt.xlabel("spherical harmonic degree (l)"); plt.legend(); plt.show()

    # generates spherical harmonics following power-law distribution
    clm = pysh.SHCoeffs.from_random(psd_power, seed=12345)

    #clm.plot_spectrum()
    #plt.show()

    # spectrum
    power = clm.spectrum()
    total_power = power.sum()
    print("  power-law: power = ",total_power)

    # discrete map
    pert_grid = clm.pad(lmax_calc).expand(grid='DH2')

    #debug
    #print("debug: ",pert_grid.data.sum(),pert_grid.data.shape,pert_grid.data.size)

    # removes the mean to have centered at zero
    avg = pert_grid.data.sum()/ (pert_grid.data.size)
    print("  power-law: mean = ",avg)
    pert_grid.data -= avg
    print("  power-law: mean corrected = ",pert_grid.data.sum()/(pert_grid.data.size))

    # normalizes range
    pert_max = np.absolute(pert_grid.data).max()
    print("  power-law: absolute max = ",pert_max)

    # sets maximum perturbation value to range [-1,1]
    pert_grid.data /= pert_max
    # scales with maximum strength
    pert_grid.data *= perturbation_max_amplitude
    # checks
    print("  power-law: normalized min/max = {} / {}".format(pert_grid.min(),pert_grid.max()))
    # plotting
    fig, ax = pert_grid.plot(colorbar='bottom',cb_label='Perturbation',title='Perturbation - power-law')
    plt.savefig("tmp/pert_powerlaw.jpg")
    if show_maps: plt.show()

    return pert_grid


def get_map_pertub_vonKarman(degrees,lmax_calc):
    """
    creates a map based on spherical harmonics with a simple von Karman distribution
    """
    global perturbation_max_amplitude,show_maps,perturbation_corr_wavelength,R_PLANET

    # von Karman distribution
    # k wavenumber: k = omega/c = 2 pi/lambda = 2 pi / (T * c)
    #
    #               for a wavelength lambda ~ 5 s * 2.0 km/s = 10 km for an S-wave
    #               k = 2 pi / 10 km -> determine a such that a * k ~ 1 with k in degrees?
    #
    # scattering regimes: k * a >> 1 (lambda << a) - high-frequency scattering
    #                     k * a ~ 1  (lambda ~  a) - diffraction condition
    #                     k * a << 1 (lambad >> a) - low-frequency scattering
    #
    # wavelength lamda: convert from km to degree -> 1 degree in km: 2 * pi * R_planet / 360
    #                   lambda_in_degrees = L / (2 * pi * R_PLANET / 360.0)
    lambda_in_degrees = perturbation_corr_wavelength / (2.0 * np.pi * R_PLANET / 360.0)

    print("  von Karman: planet radius (in km) = ",R_PLANET)
    print("  von Karman: correlation wavelength (in km)  = ",perturbation_corr_wavelength)
    print("  von Karman:                        (in deg) = ",lambda_in_degrees)


    # 1D power spectral density: Imperatori & Mai (2012)
    #
    # P(k) = sigma**2 ( 2 sqrt(pi * a) )**E gamma(H + E/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + E/2) )
    #      = sigma**2 ( 2 sqrt(pi * a) )**1 gamma(H + 1/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 1/2) )
    #      = sigma**2 * 2 * sqrt(pi * a) * gamma(H + 1/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 1/2) )
    #
    a = lambda_in_degrees / (2.0 * np.pi) # correlation length
    H = 0.3                               # Hurst exponent
    sigma = 0.1                           # standard deviation: 10%

    print("  von Karman: Hurst exponent     = ",H)
    print("  von Karman: standard deviation = ",sigma)
    print("")

    psd_vonKarman = psd_vonKarman_1D(a,H,sigma,degrees)

    # psd plot
    if show_maps:
        plt.plot(degrees,psd_vonKarman,label="psd - von Karman 1D - a={}".format(a)); plt.xscale('log');
        plt.xlabel("spherical harmonic degree (l)"); plt.legend(); plt.show()

    # generates spherical harmonics following power-law distribution
    clm = pysh.SHCoeffs.from_random(psd_vonKarman, seed=12345)

    # spectrum
    power = clm.spectrum()
    total_power = power.sum()
    print("  von Karman: spectral power = ",total_power)

    #clm = clm / total_power
    #clm = clm / (np.pi)
    #clm.plot_spectrum()
    #plt.show()

    # discrete map
    pert_grid = clm.pad(lmax_calc).expand(grid='DH2')

    #debug
    #print("debug: ",pert_grid.data.sum(),pert_grid.data.shape,pert_grid.data.size)

    # removes the mean to have centered at zero
    avg = pert_grid.data.sum()/(pert_grid.data.size)
    print("  von Karman: mean = ",avg)
    pert_grid.data -= avg
    print("  von Karman: mean corrected = ",pert_grid.data.sum()/(pert_grid.data.size))

    # normalizes
    pert_max = np.absolute(pert_grid.data).max()
    print("  von Karman: absolute max = ",pert_max)

    # sets maximum perturbation value to range [-1,1]
    pert_grid.data /= pert_max
    # scales with maximum strength
    pert_grid.data *= perturbation_max_amplitude
    # checks
    print("  von Karman: normalized min/max = {} / {}".format(pert_grid.min(),pert_grid.max()))
    # plotting
    fig, ax = pert_grid.plot(colorbar='bottom',cb_label='Perturbation',title='Perturbation - von Karman')
    plt.savefig("tmp/pert_vonKarman.jpg")
    if show_maps: plt.show()

    return pert_grid


def create_crustmaps(moho_sh_file):
    """
    creates crustmaps file for crustal model in SPECFEM3D_GLOBE format
    """
    global CRUSTMAP_RESOLUTION,num_layers_crust,show_maps, \
           add_perturbations,perturbation_max_amplitude,perturbation_corr_wavelength

    print("crustmaps:")
    print("  input moho file     : %s" % moho_sh_file)
    print("")
    print("  number of crustal layers: ",num_layers_crust)
    print("  crustmaps resolution    : ",CRUSTMAP_RESOLUTION)
    print("")
    if add_perturbations:
        print("  adding perturbations    : distribution           = ",perturbation_dist,"(0==power-law/1==vonKarman)")
        print("                            amplitude max          = ",perturbation_max_amplitude)
        if perturbation_dist == 1:
            # von Karman
            print("                            correlation wavelength = ",perturbation_corr_wavelength)
        print("")

    # work directory
    path = os.getcwd()

    # creates temporary directory ./tmp for files
    if not os.path.exists("tmp"): os.mkdir("tmp")

    # sets basename
    crust_model_type = 0
    if "mars" in moho_sh_file.lower():
        # Mars
        rootname = "marscrust"
        crust_model_type = 1
    elif "moon" in moho_sh_file.lower():
        # Moon
        rootname = "mooncrust"
        crust_model_type = 2
    else:
        # something else
        rootname = "tmpcrust"
        crust_model_type = 0
        print("Topography not specified yet for this crustal model, please use Mars or Moon crustal models for now...")
        sys.exit(1)

    # reads in moho depth file
    # mars moho thickness map (spherical harmonics)
    # e.g., moho_sh_file = "mars/Moho-Mars-DWAK-20-2550.sh"
    #
    # assumes shtools file format
    # see: https://shtools.github.io/SHTOOLS/file-formats.html
    #
    # shtools-format will look like:
    # #spherical_harmonic_degree (l) #spherical_harmonic_order (m) #coeff_0 #coeff_1
    # 0, 0, 3.3601982815853325e+06, 0.0000000000000000e+00
    # 1, 0, 5.4091504727215333e+03, 0.0000000000000000e+00
    # 1, 1, 3.1830716273821452e+02, 2.4692837341857594e+03
    # ..
    moho = pysh.SHCoeffs.from_file(moho_sh_file, format='shtools')

    # stats
    lmax = moho.lmax

    # limits expansion
    # if lmax > 90 like for Moon expansion, this becomes too expensive to compute the crustmaps grid locations
    # thus, we will limit to lmax = 90 as used for Mars expansions by default
    if lmax > 90:
        lmax = 90
        moho = moho.pad(lmax)

    # ctplanets will generate absolute moho depths (in m), with up to spherical degree lmax (e.g. = 90)
    # thus, to get the corresponding moho thickness, we need the elevation info as well.
    #
    if crust_model_type == 1:
        # Mars
        # we will read in mars topo using the pyshtools dataset (MarsTopo2600), up to a maximum degree lmax * 4,
        # which was also used to create the moho maps.
        lmax_topo = lmax * 4
        topo = pysh.datasets.Mars.MarsTopo2600(lmax=lmax_topo)
    elif crust_model_type == 2:
        # Moon
        # Moon topo was used with a maximum degree 900 to create the moho maps.
        # we will limit to 4 * lmax = 4 * 90 = 360 maximum.
        lmax_topo = lmax * 4
        topo = pysh.datasets.Moon.MoonTopo2600p(lmax=lmax_topo)
    else:
        # something else
        print("Topography not specified yet for this crustal model")
        sys.exit(1)

    # note: moho values are given in m as elevation (absolute radius) when produced by ctplanets.
    #       we will convert to km and calculate moho thickness rather than actual moho elevation.
    #
    # convert to km
    moho /= 1000.0
    topo /= 1000.0

    # lmax_calc uses a resolution to give the same lat/lon grid as with the crustmaps resolution
    # when using the expand(..) function
    #
    # note: the DH2 grid uses a nlat x nlon grid with N x 2N points (L == N/2 - 1 degree)
    #       for example, to have a grid spacing of 1x1 degree,
    #       a spherical harmonic degree L = 180/2 - 1 = 89 would be used.
    N = 180 * CRUSTMAP_RESOLUTION
    lmax_calc = int(N / 2 - 1)

    # moho: average value (l=0,m=0)
    avg = moho.coeffs[0,0,0]

    # moho: determine min/max values based on a predefined global grid (DH2)
    moho_grid = moho.pad(lmax_calc).expand(grid='DH2')

    # topo: determine min/max values based on a predefined global grid (DH2)
    topo_grid = topo.pad(lmax_calc).expand(grid='DH2')

    # thickness on lower lmax degree
    thick = topo.pad(lmax_calc) - moho.pad(lmax_calc)

    # map on global grid
    thick_grid = topo_grid - moho_grid

    # plotting
    #fig, ax = moho.plot_spectrum()
    #plt.title('Moho power spectrum')
    #plt.show()
    fig, ax = moho_grid.plot(colorbar='bottom',cb_label='Moho radius (km)',title='Moho')
    # saves as JPEG file
    plt.savefig("tmp/moho.jpg")
    if show_maps: plt.show()

    # plotting
    thick_grid.plot(colorbar='bottom',cb_label='Moho thickness (km)',title='Moho thickness')
    #saves as JPEG file
    plt.savefig("tmp/moho_thickness.jpg")
    if show_maps: plt.show()

    print("moho map:")
    print("  maximum spherical harmonic degree lmax = ",lmax)
    print("  spherical harmonic degree lmax_calc    = ",lmax_calc)
    print("")
    print("  grid points nlat / nlon                = {} / {}".format(thick_grid.nlat,thick_grid.nlon))
    print("")
    print("  average value                          = {:.2f} (km)".format(avg))
    print("  min/max value                          = {:.2f} / {:.2f} (km)".format(moho_grid.min(),moho_grid.max()))
    print("")
    print("  topo average value                     = {:.2f} (km)".format(topo.coeffs[0,0,0]))
    print("  topo min/max value                     = {:.2f} / {:.2f} (km)".format(topo_grid.min(),topo_grid.max()))
    print("")
    print("  moho thickness averag value            = {:.2f} (km)".format(thick.coeffs[0,0,0]))
    print("  moho thickness min/max value           = {:.2f} / {:.2f} (km)".format(thick_grid.min(),thick_grid.max()))
    print("")

    # moho depth at lander site
    if crust_model_type == 1:
        # Mars
        # XB.ELYSE: 4.502384        135.623444
        # from ctplanets script:
        lat_InSight = 4.502384
        lon_InSight = 135.623447
        #r_topo = topo.expand(lat=lat_InSight,lon=lon_InSight)  # from higher resolution topo
        #r_moho = moho.expand(lat=lat_InSight,lon=lon_InSight)  # from moho file
        #thickness = r_topo - r_moho
        thickness = thick.expand(lat=lat_InSight,lon=lon_InSight)
        print("  moho thickness at InSight lander       = {:.2f} (km)".format(thickness))
        print("")
    elif crust_model_type == 2:
        # Moon
        # Apollo 12 landing site S12 -3.01084     -23.42456
        lat_Apollo_12 = -3.01084
        lon_Apollo_12 = -23.42456
        #r_topo = topo.expand(lat=lat_Apollo_12,lon=lon_Apollo_12)  # from higher resolution topo
        #r_moho = moho.expand(lat=lat_Apollo_12,lon=lon_Apollo_12)  # from moho file
        #thickness = r_topo - r_moho
        thickness = thick.expand(lat=lat_InSight,lon=lon_InSight)
        print("  moho thickness at Apollo 12 site       = {:.2f} (km)".format(thickness))
        print("")

    # random perturbations
    if add_perturbations:
        # creates a random distribution map based on a defined power spectrum
        print("adding perturbations:")
        # power spectrum
        # see: https://nbviewer.org/github/SHTOOLS/SHTOOLS/blob/master/examples/notebooks/grids-and-coefficients.ipynb
        # up to degree 100
        #degrees = np.arange(101, dtype=float)
        # up to degree lmax * 4 = 90 * 4 = 360
        degrees = np.arange(lmax * 4, dtype=float)
        if perturbation_dist == 0:
            # power-law w/ degree -2
            degrees[0] = np.inf
            pert_grid = get_map_pertub_powerlaw(degrees,lmax_calc)
        elif perturbation_dist == 1:
            # von Karman
            pert_grid = get_map_pertub_vonKarman(degrees,lmax_calc)
        else:
            print("Invalid perturbation distribution, choose 0 == power-law / 1 == von Karman")
            sys.exit(1)
        print("")

    # loops over lat/lon
    #
    # Values are ordered on a grid that runs from 90 degrees lat (0 degrees lng to 360 degrees lng ) to -90 degrees lat.
    # Every line in a crustmap files sets all values for one specific latitude:
    #
    #    (length of one line = 360 * 4 entries)
    #    (total number of lines = 180 * 4)
    #
    # The value for (lat = 89.875, long = 0.125) is stored in entry (line 1,column 1).
    # You get the value that is north west of your (lat, long) point by using:
    #
    #    line = int(0.5+((90.0-lat)*CRUSTMAP_RESOLUTION))
    #    column = int0.5+(xlng*CRUSTMAP_RESOLUTION)
    #
    print("looping over lat/lon grid")

    # crustmaps:
    #
    #   topography crustmaps (marscrustt3.cmap - marscrustt7.cmap)
    #   s-velocity crustmaps (marscrusts3.cmap - marscrusts7.cmap)
    #   p-velocity crustmaps (marscrustp3.cmap - marscrustp7.cmap)
    #   density    crustmaps (marscrustr3.cmap - marscrustr7.cmap)
    #
    # with 5-layers:
    #   marscrust**3 - soft sediment layer thickness
    #   marscrust**4 - hard sediment layer
    #   marscrust**5 - upper crust layer
    #   marscrust**6 - middle crust layer
    #   marscrust**7 - lower crust layer
    #
    # numpy arrays to hold crustal values for all lat/longitudinal values at a certain lat
    # note: numpy arrays use C-order indexing, i.e., the last index represents the most rapidly changing one,
    #       which is opposite for Fortran arrays, where the first index is the inner-loop index.
    #       -> i use 8 for the left, and ilon on the last index to have a faster copy of a single line
    #          when writing out to files: arr = lon_line_t[ilayer,:]
    all_val_t = np.zeros((5,180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION))
    all_val_p = np.zeros((5,180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION))
    all_val_s = np.zeros((5,180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION))
    all_val_r = np.zeros((5,180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION))

    for ilat in range(0,180*CRUSTMAP_RESOLUTION):
        # values in fortran: [1,180*CRUSTMAP_RESOLUTION],
        #        in python : [0,180*CRUSTMAP_RESOLUTION-1]

        # numpy arrays to hold crustal values for all longitudinal values at a certain lat
        # note: numpy arrays use C-order indexing, i.e., the last index represents the most rapidly changing one,
        #       which is opposite for Fortran arrays, where the first index is the inner-loop index.
        #       -> i use 8 for the left, and ilon on the last index to have a faster copy of a single line
        #          when writing out to files: arr = lon_line_t[ilayer,:]
        lon_line_t = np.zeros((5,360*CRUSTMAP_RESOLUTION))
        lon_line_p = np.zeros((5,360*CRUSTMAP_RESOLUTION))
        lon_line_s = np.zeros((5,360*CRUSTMAP_RESOLUTION))
        lon_line_r = np.zeros((5,360*CRUSTMAP_RESOLUTION))

        for ilon in range(0,360*CRUSTMAP_RESOLUTION):
            # values in fortran: [1,360*CRUSTMAP_RESOLUTION]

            # increment in degrees
            incr = 1.0 / CRUSTMAP_RESOLUTION

            # converts to ilat from lat
            # see model_crustmaps.f90, ibilinearmap() routine
            # range: lat in [-90,90]
            # fortran: ilat = 0.5 + ((90.0-lat)*CRUSTMAP_RESOLUTION)
            #
            # -> converts to lat from ilat
            lat = 90.0 - ((ilat+1) - 0.5)/CRUSTMAP_RESOLUTION

            # converts to ilon from lon
            # see model_crustmaps.f90, ibilinearmap() routine
            # range: xlng in [0,360]
            # fortran: ilon = 0.5 + (xlng*CRUSTMAP_RESOLUTION)
            #
            # -> converts to lon from ilon
            lon = ((ilon+1) - 0.5)/CRUSTMAP_RESOLUTION

            # user output
            if ilon == 0 and ilat % 10*CRUSTMAP_RESOLUTION == 0: print("  ...lat ",lat)

            # extracts moho depth (absolute value as radius in m)
            #r_moho = moho.expand(lat=lat,lon=lon)
            #r_moho = moho_grid.data[ilat,ilon]
            # extracts topo elevation
            #r_topo = topo.expand(lat=lat,lon=lon)
            #r_topo = topo_grid.data[ilat,ilon]
            # moho thickness in m
            #moho_thickness = r_topo - r_moho

            # moho thickness
            moho_thickness = thick_grid.data[ilat,ilon]

            #debug
            #print("debug: lat/lon {}/{} - moho thickness {} - r_topo {} r_moho {}".format(lat,lon,moho_thickness,r_topo,r_moho))
            # debug compare lower degree thick map vs. higher degree topo value
            #print("debug: moho thickness: ",moho_thickness,thick.expand(lat=lat,lon=lon))

            # checks - r_topo should always be higher than r_moho
            if moho_thickness <= 0.0:
                print("Error: invalid moho thickness ",moho_thickness)
                print("       lat/lon = {} / {}".format(lat,lon))
                #print("       radius topo    = ",r_topo)
                #print("       radius moho    = ",r_moho)
                sys.exit(1)

            # layer values
            #
            # with 5-layers:
            #   marscrust**3 - soft sediment layer thickness
            #   marscrust**4 - hard sediment layer
            #   marscrust**5 - upper crust layer
            #   marscrust**6 - middle crust layer
            #   marscrust**7 - lower crust layer
            #
            # we won't assign sediments
            # -> will put sediment layer thickness and values to zero
            #
            # initializes layer values
            val_t = np.zeros(5)
            val_p = np.zeros(5)
            val_s = np.zeros(5)
            val_r = np.zeros(5)

            if num_layers_crust == 2:
                # assigns upper/lower crust
                # according to plan A, fig. 4 in Knapmeyer-Endrun (2021, Science)
                #
                # ctplanets - constant InSight crust output
                # summary-20.txt:
                # Model    rho_c    rho_mantle    t_insight    t_ave     t_min    t_max
                # DWAK    2550.000000    3380.160000    20.000184    29.278954    2.053932    60.580403
                #
                # thickness assigned proportional to lander site estimates:
                # plan A, fig. 4 uses an upper thickness of 9km,
                #                        lower           of 20-9 = 11km
                #
                # -> factor upper = 9 km / 20 km
                #    factor lower = 11 km / 20 km
                factor_up = 9.0 / 20.0
                factor_low = 11.0 / 20.0

                t_up = moho_thickness * factor_up
                t_low = moho_thickness * factor_low

                # vp,vs,rho values based on mars_1D.dat model
                # upper crust - layer 5 -> index 5-3 = 2
                val_t[2] = t_up    # thickness km
                val_p[2] = 3.80    # vp km/s
                val_s[2] = 1.85    # vs km/s
                val_r[2] = 2.30    # density kg/m^3
                # skips middle crust
                # lower crust - layer 6 -> index 6-3 = 3
                val_t[4] = t_low   # km      - values from middle crust of mars_1D.dat
                val_p[4] = 4.50    # km/s
                val_s[4] = 2.80    # km/s
                val_r[4] = 2.57    # kg/m^3

            elif num_layers_crust == 3:
                # assigns upper/middle/lower crust
                # according to plan B, fig. 4 in Knapmeyer-Endrun (2021, Science)
                #
                # ctplanets - constant InSight crust output
                # summary-39.txt:
                # Model    rho_c    rho_mantle    t_insight    t_ave     t_min    t_max
                # DWAK    2550.000000    3380.160000    38.997634    48.727003    18.837998    81.342647
                #
                # thickness assigned proportional to lander site estimates:
                # plan B, fig. 4 uses an upper thickness of 9km,
                #                        middle of 24-9 = 15km,
                #                        lower of 39-24 = 15km
                # -> factor upper = 9 km / 39 km
                #    factor middle = 15 km / 39 km
                #    factor lower = 15 km / 39 km
                factor_up = 9.0 / 39.0
                factor_mid = 15.0 / 39.0
                factor_low = 15.0 / 39.0

                t_up = moho_thickness * factor_up
                t_mid = moho_thickness * factor_mid
                t_low = moho_thickness * factor_low

                # vp,vs,rho values based on mars_1D.dat model
                # upper - layer 5 -> index 5-3 = 2
                val_t[2] = t_up    # km
                val_p[2] = 3.80    # km/s
                val_s[2] = 1.85    # km/s
                val_r[2] = 2.30    # kg/m^3
                # middle - layer 6 -> index 5-3 = 3
                val_t[3] = t_mid   # km
                val_p[3] = 4.50    # km/s
                val_s[3] = 2.80    # km/s
                val_r[3] = 2.57    # kg/m^3
                # lower - layer 7 -> index 7-3 = 4
                val_t[4] = t_low   # km
                val_p[4] = 6.22    # km/s
                val_s[4] = 3.75    # km/s
                val_r[4] = 2.86    # kg/m^3

            else:
                print("Invalid number of crustal layers ",num_layers_crust," - not implemented yet")
                sys.exit(1)

            # perturbs velocities
            if add_perturbations:
                # perturbation
                #pert = clm.expand(lat=lat,lon=lon)
                pert = pert_grid.data[ilat,ilon]

                # upper crust layer - perturbation
                # vp
                val_p[2] = (1.0 + pert) * val_p[2]
                # vs
                val_s[2] = (1.0 + pert) * val_s[2]

            # sets array values for line
            for ilayer in range(0,5):
                lon_line_t[ilayer,ilon] = val_t[ilayer]
                lon_line_p[ilayer,ilon] = val_p[ilayer]
                lon_line_s[ilayer,ilon] = val_s[ilayer]
                lon_line_r[ilayer,ilon] = val_r[ilayer]

        # stores all ilat values
        for ilayer in range(0,5):
            all_val_t[ilayer,ilat,:] = lon_line_t[ilayer,:]
            all_val_p[ilayer,ilat,:] = lon_line_p[ilayer,:]
            all_val_s[ilayer,ilat,:] = lon_line_s[ilayer,:]
            all_val_r[ilayer,ilat,:] = lon_line_r[ilayer,:]

    # crustmaps:
    #
    #   topography crustmaps (marscrustt3.cmap - marscrustt7.cmap)
    #   s-velocity crustmaps (marscrusts3.cmap - marscrusts7.cmap)
    #   p-velocity crustmaps (marscrustp3.cmap - marscrustp7.cmap)
    #   density    crustmaps (marscrustr3.cmap - marscrustr7.cmap)
    #
    # with 5-layers:
    #   marscrust**3 - soft sediment layer thickness
    #   marscrust**4 - hard sediment layer
    #   marscrust**5 - upper crust layer
    #   marscrust**6 - middle crust layer
    #   marscrust**7 - lower crust layer
    print("")
    print("saving files...")
    os.chdir("tmp")

    for ilayer in range(0,5):
        # layer start at 3 -> array index start at 0: layer index 3 - 3 = 0 ilayer
        index = ilayer + 3

        name = rootname+"t{}.cmap".format(index)
        print("  file ",name)
        arr = all_val_t[ilayer,:,:]
        np.savetxt(name,arr,delimiter=' ',fmt='%10.3f')

        name = rootname+"p{}.cmap".format(index)
        print("  file ",name)
        arr = all_val_p[ilayer,:,:]
        np.savetxt(name,arr,delimiter=' ',fmt='%10.3f')

        name = rootname+"s{}.cmap".format(index)
        print("  file ",name)
        arr = all_val_s[ilayer,:,:]
        np.savetxt(name,arr,delimiter=' ',fmt='%10.3f')

        name = rootname+"r{}.cmap".format(index)
        print("  file ",name)
        arr = all_val_r[ilayer,:,:]
        np.savetxt(name,arr,delimiter=' ',fmt='%10.3f')

        print("")

    print("see files written in: ./tmp")

    # changes back to work directory
    os.chdir(path)

    print("")
    print("stats:")
    print("  soft sediment layer: vp min/max = {} / {} (km/s)".format(all_val_p[0,:,:].min(),all_val_p[0,:,:].max()))
    print("  hard sediment layer: vp min/max = {} / {} (km/s)".format(all_val_p[1,:,:].min(),all_val_p[1,:,:].max()))
    print("  upper  crust layer : vp min/max = {} / {} (km/s)".format(all_val_p[2,:,:].min(),all_val_p[2,:,:].max()))
    print("  middle crust layer : vp min/max = {} / {} (km/s)".format(all_val_p[3,:,:].min(),all_val_p[3,:,:].max()))
    print("  lower  crust layer : vp min/max = {} / {} (km/s)".format(all_val_p[4,:,:].min(),all_val_p[4,:,:].max()))
    print("")

    # plot vp
    arr = all_val_p[2,:,:]
    total_max = np.absolute(arr).max()
    total_min = np.absolute(arr).min()
    avg = arr.sum()/arr.size
    print("  vp - upper crust layer: min/max = {}/{}".format(total_min,total_max))
    print("                          average = {}".format(avg))

    vp_grid = pysh.SHGrid.from_array(arr)
    vp_grid.plot(colorbar='bottom',cb_label='Vp (km/s)',title='Vp - upper crust',
                 cmap="seismic",cmap_limits=[total_min,total_max])
    plt.savefig("tmp/vp_upper_crust_grid.jpg")
    if show_maps: plt.show()

    print("")
    print("all done")


if __name__ == '__main__':
    # gets arguments
    if '--help' in sys.argv or '-h' in sys.argv or len(sys.argv) != 2:
        print("usage: ./run_create_a_crustmap.py moho_sh_file [Moho depths in spherical harmonics format, e.g. == mars/Moho-Mars-DWAK-20-2550.sh]")
        sys.exit(1)
    else:
        moho_sh_file = sys.argv[1]

    create_crustmaps(moho_sh_file)

