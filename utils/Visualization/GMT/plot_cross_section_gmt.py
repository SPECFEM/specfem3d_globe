#!/usr/bin/env python
#
#
import os
import sys

import pygmt
import numpy as np
import pandas as pd
#from io import StringIO

######################################################################
# USER Parameters

# True for surface interpolation, False for plotting dots
use_surface_interpolation = True

# small perspective globe (on left top corner)
add_perspective = True

# use hillshade topo as background
use_topo = True

# True for vertical, False for horizontal cross-section
plot_vert = False

# show figure
show_plot = False


######################################################################

# globals
color_range = None   # specify a colormap min/max
fix_region = None    # user specified region code

# to convert radius to depth
EARTH_SURFACE_RADIUS = 6371.0  # PREM surface radius in km


def process_data(filename):
    """
    reads vertical/horizontal cross-section data
    """
    global plot_vert,use_surface_interpolation
    global EARTH_SURFACE_RADIUS

    # pandas read_csv doesn't work if file has spacings before '#'...
    #data = pd.read_csv(filename, sep=' ', comment='#', header=None)

    # numpy reads tabled data
    data_array = np.loadtxt(filename)

    # check length
    if len(data_array) == 0:
        print("no data points")
        print("Please check file...")
        sys.exit(1)

    # convert to pandas dataframe
    data = pd.DataFrame(data_array)

    # file format: (lon, lat, radius, value)
    if plot_vert:
        # Project data onto vertical plane
        lon_min, lon_max = data_array[:,0].min(), data_array[:,0].max()
        lat_min, lat_max = data_array[:,1].min(), data_array[:,1].max()

        data_xyz = data_array[:,0:3]    # lon/lat/radius
        #print("data: xyz",data_xyz)

        # project coordinates on (theta,r) plane
        data_xyz = pygmt.project(data=data_xyz, center=[lon_min, lat_min], endpoint=[lon_max, lat_max])
        #print("data: xy projected ",data_xyz)

        # Extract data (lon, lat, radius, value)
        data = data[[0, 1, 2, 3]]
        data.columns = ['lon', 'lat', 'radius', 'value']

        # (p,q) coordinates in plane; (r,s) coordinates in x/y system
        data_xyz = data_xyz[[0, 1, 2, 3, 4, 5, 6]]
        data_xyz.columns = ['lon', 'lat', 'radius', 'p', 'q', 'r', 's']

        # get epicentral distance
        # adds difference in coordinate s as epicentral distance
        data['distance'] = data_xyz['s'] - data_xyz['s'].min()

        # get depth
        data['depth'] = EARTH_SURFACE_RADIUS - data['radius']

        #epi_dist = data['distance']
        #print("data: epicentral dist ",epi_dist)

    else:
        # Extract data (lon, lat, value)
        data = data[[0, 1, 3]]
        data.columns = ['lon', 'lat', 'value']

    return data


def create_colormap(val_average):
    """
    creates colormap
    """
    # modified spectral w/ white at zero and continuous steps
    # http://soliton.vm.bytemark.co.uk/pub/cpt-city/cb/div/tn/Spectral_11.png.index.html
    txt = """# COLOR_MODEL = RGB\n
-5.50 058 001 016  -4.50 158 001 066\n
-4.50 158 001 066  -3.50 213 062 079\n
-3.50 213 062 079  -2.50 244 109 067\n
-2.50 244 109 067  -1.50 253 204 097\n
-1.50 253 204 097  -0.50 254 244 139\n
-0.50 254 244 139   0.0  255 255 255\n
 0.0  255 255 255   0.50 200 245 152\n
 0.50 200 245 152   1.50 171 221 164\n
 1.50 171 221 164   2.50 102 194 165\n
 2.50 102 194 165   3.50 050 136 199\n
 3.50 050 136 199   4.50 094 079 162\n
 4.50 094 079 162   5.50 054 039 122\n
 """

    cptname = filename + ".cpt"
    with open(cptname,'w') as f:
        f.write(txt)

    # center color map around average value
    if not val_average is None:
        col_min = val_average - 0.5  # 2.0   # val_min
        col_max = val_average + 0.5  # 3.0   # val_max
    else:
        col_min = 2.0   # defaults
        col_max = 3.0

    # fixed
    if not color_range is None:
        col_min = color_range[0]   # val_min
        col_max = color_range[1]   # val_max
        print("fixing colormap range: {} / {}".format(col_min,col_max))
        print("")

    # colormap
    pygmt.makecpt(cmap=cptname, series=[col_min, col_max, 0.01 * (col_max - col_min)],background=True)


def add_topo_background(dim_max,R,fig):
    """
    adds hillshaded topography on land masses
    """
    print("adding topography")
    print("")
    # Here, '@earth_relief_01m' is a built-in dataset in GMT/PyGMT
    #fig.grdimage(grid='@earth_relief_01m', cmap="gray", region=R)

    # hillshade
    # Load sample grid in target area
    if dim_max < 10.0:
        topo = pygmt.datasets.load_earth_relief(resolution="30s", region=R)
    else:
        topo = pygmt.datasets.load_earth_relief(resolution="02m", region=R)

    # calculate the reflection of a light source projecting from west to east
    # (azimuth of 300 degrees) and at a latitude of 40 degrees from the horizon
    dgrid = pygmt.grdgradient(grid=topo, radiance=[300, 40], verbose='q')  # normalized between [-1,1]

    #pygmt.makecpt(cmap="gray", series=[-2.0, 2.0, 0.02])
    fig.grdimage(grid=dgrid, cmap="gray", region=R, verbose='q')
    fig.coast(water="#FFFFFF", verbose='q') # mask out ocean


def add_perspective_view(lon_min,lon_max,lat_min,lat_max,fig):
    """
    adds small globe on top left corner w/ a perspective view
    """
    print("adding perspective globe")
    print("")

    # Define the projection. 'Glon0/lat0' is a perspective projection.
    # Adjust 'lon0' and 'lat0' to the longitude and latitude of your first point.
    projection = f"G{lon_min}/{lat_min}/4c"

    # Define the region. 'g' stands for global.
    region_globe = 'g'

    # Define the frame. '5g15' sets the annotation and gridline interval.
    frame_globe = '5g15'

    with fig.inset(position="jTL+w4c+o0.1c"):
        # Plot the coastlines
        #fig.basemap(region=region, projection=projection, frame=frame)
        fig.coast(region=region_globe, projection=projection, frame=frame_globe, land="lightbrown", water="lightblue")
        if plot_vert:
            # Plot the great-circle line
            fig.plot(x=[lon_min, lon_max], y=[lat_min, lat_max], pen='0.1')

            # Plot the points
            fig.plot(x=[lon_min, lon_max], y=[lat_min, lat_max], style='c0.3', fill='100/100/100', pen='1')

            # Add text annotations
            fig.text(x=lon_min, y=lat_min, text='A', angle=15, justify='RT')
            fig.text(x=lon_max, y=lat_max, text='B', angle=15, justify='RT')
        else:
            # plot lines
            data = np.array([[lon_min, lat_min], [lon_min, lat_max], [lon_max, lat_max], [lon_max, lat_min]])
            fig.plot(data=data, pen="1p", close=True)


def plot_cross_section(filename):
    """
    plot cross section w/ GMT
    """
    global plot_vert,use_surface_interpolation
    global add_perspective

    # user output
    print("")
    print("plotting cross section:")
    print("  file: ",filename)
    print("")
    if plot_vert:
        print("  using vertical cross-section")
    else:
        print("  using horizontal cross-section")
    if use_surface_interpolation:
        print("  using surface interpolation")
    print("")

    # check file
    if not os.path.isfile(filename):
        print("Please check if file exists: ",filename)
        sys.exit(1)

    # read header to determine cross-section type horizontal/vertical
    with open(filename,'r') as f:
        print("file header:")
        lines = f.readlines()

    val_average = None
    val_depth = None

    for line in lines:
        #if "#" in line: print("line: ",line.strip())
        #line = line.strip()
        if "#" in line and "horizontal" in line:
            print("   horizontal cross-section")
            plot_vert = False

        if "#" in line and "vertical" in line:
            print("   vertical cross-section")
            plot_vert = True

        if "#" in line and "depth" in line:
            val_depth = float(line.split('=')[1])
            print("   depth: ",val_depth)

        if "#" in line and "cross-section average" in line:
            val_average = float(line.split('=')[1])
            print("   cross-section average: ",val_average)
    print("")

    # parameter name based on file name format
    #cross_section_***.dat
    name = os.path.basename(filename)    # cross_section_***.dat
    name = name.split('.')[0]            # cross_section_***
    parameter_name = name.split('cross_section_')[1]   # vsv
    parameter_name = parameter_name.strip().capitalize()         # Vsv

    print("parameter: ",parameter_name)
    print("")

    # Read data and extract statistics
    data = process_data(filename)
    #print("data: ",data)

    # main area
    lon_min, lon_max, lat_min, lat_max = data['lon'].min(), data['lon'].max(), \
                                         data['lat'].min(), data['lat'].max()
    if plot_vert:
        # radius/depth
        r_min, r_max, depth_min, depth_max = data['radius'].min(), data['radius'].max(), \
                                             data['depth'].min(), data['depth'].max()
        # epicentral distance
        dist_min, dist_max = data['distance'].min(), data['distance'].max()

    # data values (Vsv,..)
    val_min, val_max = data['value'].min(), data['value'].max()

    print("data statistics: ")
    print("  lon    : min/max = {:6.2f} / {:6.2f}".format(lon_min,lon_max))
    print("  lat    : min/max = {:6.2f} / {:6.2f}".format(lat_min,lat_max))
    if plot_vert:
        print("  radius : min/max = {:12.2f} / {:12.2f}".format(r_min,r_max))
        print("  depth  : min/max = {:12.2f} / {:12.2f}".format(depth_min,depth_max))
        print("  epicentral distance : min/max = {:12.2f} / {:12.2f}".format(dist_min,dist_max))
    print("")
    print("  data   : min/max = {} / {}".format(val_min,val_max))
    print("")

    # Define regions, projections, and annotations based on script and data
    # default global
    if plot_vert:
        # vertical
        R = '0/360/3480/6371'                             # Vertical region (radius CMB to surface)
        J = 'P6i'                                         # Projection (polar view for full globe)
        B = '45a90/10a45WeSn:.\"Cross-Section\":+tle'     # Basemap annotation
        title = 'Vertical Cross-Section'
    else:
        # horizontal
        R = 'd'                                           # Horizontal region
        J = 'M6i'                                         # Projection
        B = '15a15/100a100WeSn:.\"Cross-Section\":+tle'   # Basemap annotation
        title = 'Horizontal Cross-Section'

    # dimensions
    dim_lon = abs(lon_max - lon_min)    # in deg
    dim_lat = abs(lat_max - lat_min)
    dim_max = max(dim_lon,dim_lat)

    if plot_vert:
        dim_radius = abs(r_max - r_min)     # in km
        dim_depth = abs(depth_max - depth_min)     # in km
        dim_dist = abs(dist_max - dist_min) # in km

    print("dimensions:")
    print("  range lon      = ",dim_lon,"(deg)")
    print("  range lat      = ",dim_lat)
    if plot_vert:
        print("  range radius   = ",dim_radius,"(km)")
        print("  range depth    = ",dim_depth,"(km)")
        print("  range distance = ",dim_dist,"(deg)")
    print("")

    # determine if regional plot
    if max(dim_lon,dim_lat) < 180.0:
        is_regional = True
    else:
        is_regional = False

    # regional plots
    if is_regional:
        print("  using regional plot")
        print("")
        if plot_vert:
            # vertical
            R = f"{dist_min}/{dist_max}/{depth_min}/{depth_max}"     # uses depth instead of radius
            J = 'X8i/-4i'                                            # flips depth direction (-4i), having zero depth on top
            B = '5a10/5a10WeSn:.\"Cross-Section\":+tle'
        else:
            # horizontal
            R = f"{lon_min}/{lon_max}/{lat_min}/{lat_max}"
            B = '5a10/5a10WeSn:.\"Cross-Section\":'

    # fixes region
    if not fix_region is None:
        R = fix_region

    print("plot region: ",R)
    print("")

    # adapt data for vertical plots
    if plot_vert:
        # vertical
        if is_regional:
            # regional format: radius/epidistance/value
            data_xyz = data[['distance', 'depth', 'value']]
            data_xyz.columns = ['x', 'y', 'z']
        else:
            # global format: radius/epidistance/value
            data_xyz = data[['distance', 'radius', 'value']]
            data_xyz.columns = ['x', 'y', 'z']
    else:
        # horizontal
        data_xyz = data[['lon', 'lat', 'value']]
        data_xyz.columns = ['x', 'y', 'z']

    # region info
    #region_info = pygmt.info(data=data_xyz,per_column=True)
    #print("region info: ",region_info)
    #print("")

    # Create figure and set basic elements
    fig = pygmt.Figure()

    # basemap
    if plot_vert:
        # vertical
        if is_regional:
            # reverse depth direction to have depth 0 on top
            fig.basemap(region=R, projection=J, frame=['x+l"epicentral distance (deg)"','y+l"depth (km)"',f"+t{title}"])
        else:
            fig.basemap(region=R, projection=J, frame=['x+l"epicentral distance (deg)"','y+l"radius (km)"',f"+t{title}"])
    else:
        # horizontal
        fig.basemap(region=R, projection=J, frame=['a',f"+t{title}"])

    # Add topography as background only for land masses
    if use_topo and not plot_vert:
        add_topo_background(dim_max,R,fig)

        # overlay next image w/ 90% transparency
        if use_surface_interpolation:
            transparency = 30
        else:
            transparency = 90  # points need more transparency
    else:
        # no hillshade topo
        transparency = 0

    # Create colormap based on value range
    create_colormap(val_average)

    if use_surface_interpolation:
        # surface interpolation
        print("surface interpolation:")

        # use plot region for interpolation
        region = R

        # Create a grid
        if plot_vert:
            # vertical
            space_x = dim_dist / 500.0
            space_y = dim_depth / 500.0

            print("  grid: spacing = {} / {}".format(space_x,space_y))
            print("")
            
            grid = pygmt.surface(data_xyz, region=R, spacing=[space_x,space_y], verbose='q')
        else:
            # horizontal lon/lat/radius
            # determine grid spacing
            spacing = max(lon_max-lon_min,lat_max-lat_min) / 2000.0
            spacing = np.round(spacing, -int(np.floor(np.log10(np.abs(spacing))))) # limit to significant digit 0.12 -> 0.1
            if spacing < 0.05: spacing = 0.05

            print("  grid: spacing = ",spacing)
            print("")

            grid = pygmt.surface(data_xyz, region=R, spacing=spacing, verbose='q')
        #print("grid: ",grid)

        # determine which areas have no data points to dimm outside areas
        if not plot_vert:
            from scipy.spatial import cKDTree
            # Create coordinate arrays from the grid
            grid_lon, grid_lat = np.meshgrid(grid['x'], grid['y'])

            # Flatten the coordinate arrays
            grid_lon_flat = grid_lon.ravel()
            grid_lat_flat = grid_lat.ravel()

            # Create a cKDTree object for the data points
            tree = cKDTree(data[['lon', 'lat']])

            # Find the distance to the nearest data point for each grid point
            distances, _ = tree.query(np.vstack([grid_lon_flat, grid_lat_flat]).T)

            # Reshape the distances back to the shape of the grid
            distances = distances.reshape(grid.shape)
            #print("distances: ",distances)

            # Define a distance threshold
            threshold = 1.0  # Adjust this value based on your specific criteria

            # Create a mask that is True for grid points that are too far away from a data point
            mask = (distances > threshold)
            #print("mask: ",mask)

            # Create a new grid that is filled with nan
            nan_grid = np.full(grid.shape, np.nan)

            # Copy the data from your original grid to the new grid, but only for cells with data
            nan_grid = grid.where(~mask)

            # replace grid
            grid = nan_grid

        # Plot the data
        fig.grdimage(grid, cmap=True, transparency=transparency)

    else:
        # point plot
        print("point plot:")

        # determine point size based on data increments
        # difference between 2 consecutive entries
        diffs_x = np.abs(np.diff(data_xyz['x']))
        diffs_y = np.abs(np.diff(data_xyz['y']))
        #print("debug: diffs_x ",diffs_x)
        #print("debug: diffs_y ",diffs_y)

        # get increments
        dx = next((diff for diff in diffs_x if diff > 0.00001), None)
        if dx is None: dx = 1.0
        dy = next((diff for diff in diffs_y if diff > 0.00001), None)
        if dy is None: dy = 1.0
        print("  increments dx/dy = {:.2f} / {:.2f}".format(dx,dy))

        # determine point size
        dmin = min(dx,dy)
        if dmin < 0.5:
            point_size = 0.1
        elif dmin < 1.0:
            point_size = 0.2
        else:
            point_size = 0.5

        if plot_vert:
            # square points
            point_style = f"s{point_size}c"
        else:
            # circle points
            point_style = f"c{point_size}c"

        print("  point size  = {} - style = {}".format(point_size,point_style))
        print("")

        # Plot data points directly
        fig.plot(x=data_xyz['x'], y=data_xyz['y'], fill=data_xyz['z'], style=point_style, cmap=True, transparency=transparency)

    # adds a perspective globe view
    if add_perspective:
        add_perspective_view(lon_min,lon_max,lat_min,lat_max,fig)

    # coast lines
    if not plot_vert:
        fig.coast(shorelines="0.5p,black")

    # Add a color bar
    #fig.colorbar(frame=['a0.5f0.25',f"y+l\"{parameter_name} (km/s)\""],position="JBC+w5c/0.2c+h") # smaller colorbar
    fig.colorbar(frame=['a0.5f0.25',f"y+l\"{parameter_name} (km/s)\""])

    # save figure as jpeg image
    name = filename + ".jpg"
    fig.savefig(name, crop=True, dpi=720)
    print("")
    print("  figure plotted to: ",name)
    print("")

    # show figure plot
    if show_plot:
        fig.show()

    print("")
    print("all done")
    print("")


def usage():
    print("usage: ./plot_cross_section_gmt.py --file=filename [--horiz] [--vert] [--color_range=min,max] [--region=R]")
    print("                                                   [--surf] [--points] [--topo] [--no-topo] [--show]")
    print("  with")
    print("     --file=filename    - input cross-section, for example 'OUTPUT_FILES/cross_section_vpv.dat'")
    print("     --horiz            - (optional) horizontal cross-section plot (used by default)")
    print("     --vert             - (optional) vertical cross-section plot")
    print("     --surf             - (optional) use surface interpolation (default is on)")
    print("     --points           - (optional) use points for plotting (default is off, using surface interpolation by default)")
    print("     --color_range=..   - (optional) fixes colormap range to (min,max) values")
    print("     --region=R         - (optional) use a fixed region specifier R (e.g. 'lonmin/lonmax/latmin/latmax')")
    print("     --show             - (optional) show figure plot (default is off)")
    print("     --topo             - (optional) use hillshaded topography as land background (default is on)")
    print("     --no-topo          - (optional) empty land (no topography as land background), uses only coast lines")
    sys.exit(1)


if __name__ == '__main__':
    # defaults
    filename = ""

    # gets arguments
    if len(sys.argv) <= 1:
        usage()

    # reads arguments
    i = 0
    for arg in sys.argv:
        i += 1
        #print("arg: ",i,arg)
        # get arguments
        if "--file=" in arg:
            filename = arg.split('=')[1]
        elif "--help" in arg:
            usage()
        elif "--horiz" in arg:
            plot_vert = False
        elif "--vert" in arg:
            plot_vert = True
        elif "--surf" in arg:
            use_surface_interpolation = True
        elif "--no-surf" in arg:
            use_surface_interpolation = False
        elif "--points" in arg:
            use_surface_interpolation = False
        elif "--show" in arg:
            show_plot = True
        elif "--topo" in arg:
            use_topo = True
        elif "--no-topo" in arg:
            use_topo = False
        elif "--color_range" in arg:
            str_array = arg.split('=')[1]
            color_range = np.array([float(val) for val in str_array.strip('()[]').split(',')])
        elif "--region" in arg:
            fix_region = arg.split('=')[1]
        elif i > 1:
            print("argument not recognized: ",arg)
            sys.exit(1)

    plot_cross_section(filename)

