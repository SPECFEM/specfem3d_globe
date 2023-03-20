# class to read in moviedata files
# e.g. OUTPUT_FILES/moviedata000200
#
# requires xyz file OUTPUT_FILES/moviedata_xyz.bin
#
from __future__ import print_function

import os
import sys
import array
import time

import numpy as np

## class for moviedata

class moviedata(object):
    """
    class for moviedata
    """
    def __init__(self,verbose=False):
        ## initializations

        # sets default custom_real type
        # defines float 'f' or double precision 'd' for binary values
        self.custom_type = 'f'

        # value type: 1==Z, 2==N, 3==E, 4==norm
        self.use_component = 1

        # smoothing
        self.use_smoothing = False

        # smoothing kernel size
        self.use_smoothing_kernel_size = 2

        # verbosity
        self.verbose = verbose

    def __str__(self):
        info = "helper class for moviedata files"

        return info


    def read_marker(self,file,verbose=False):
        """
        reads array marker from fortran binary file
        (each time a fortran write(*)-routine call is executed, a marker at beginning and end is placed)
        """
        binlength = array.array('i')
        binlength.fromfile(file,1)

        if verbose:
            print("marker length = ",binlength[0])

        return binlength[0]


    def read_binary_file_custom_real_array(self,file,verbose=False):
        """
        reads data array from file
        """
        # gets array length in bytes
        binlength = self.read_marker(file)

        if self.custom_type == 'f':
            # float (4 bytes) for single precision
            num_points = int(binlength / 4)
        else:
            # double precision
            num_points = int(binlength / 8)

        # user output
        if verbose:
            print("  array length = ",binlength," Bytes")
            print("  number of points in array = ",num_points)
            print("")

        # reads in array data
        binvalues = array.array(self.custom_type)
        binvalues.fromfile(file,num_points)

        # fortran binary file: file ends with array-size again
        # checks markers
        binlength_end = self.read_marker(file)
        if binlength_end != binlength:
            print("Error array markers in fortran binary file:",binlength,binlength_end)
            print("start array length = ",binlength," Bytes")
            print("final array length = ",binlength_end," Bytes")
            print("array lengths should match, please check your file")
            raise Exception('array markers invalid')

        data = list()
        data.append(binvalues)

        # returns data from list-output (first item)
        # converts list to numpy array
        data_array = np.array(data[0],dtype=self.custom_type)

        return data_array


    def read_binary_moviedata_file(self,filename,verbose=False):
        """
        reads moviedata000200,.. arrays
        """
        if verbose: print("binary file: ",filename)

        with open(filename,'rb') as f:
            # fortran file format: binary
            # for example by:
            # ! wavefield values (src/specfem3D/write_movie_surface.f90)
            # write(IOUT) store_val_ux_all
            # write(IOUT) store_val_uy_all
            # write(IOUT) store_val_uz_all

            data_x = self.read_binary_file_custom_real_array(f,verbose)
            data_y = self.read_binary_file_custom_real_array(f,verbose)
            data_z = self.read_binary_file_custom_real_array(f,verbose)

        return data_x,data_y,data_z


    def read_binary_moviedata_xyz_file(self,filename,verbose=False):
        """
        reads moviedata_xyz.bin array
        """
        if verbose: print("binary file: ",filename)

        with open(filename,'rb') as f:
            # fortran file format: binary
            # for example by:
            # (see src/specfem3D/write_movie_surface.f90)
            # ! point coordinates
            # ! (given as r theta phi for geocentric coordinate system)
            # write(IOUT) store_val_x_all
            # write(IOUT) store_val_y_all
            # write(IOUT) store_val_z_all

            r = self.read_binary_file_custom_real_array(f,verbose)
            theta = self.read_binary_file_custom_real_array(f,verbose)
            phi = self.read_binary_file_custom_real_array(f,verbose)

        return r,theta,phi


    def get_data_component(self,r,theta,phi,data_x,data_y,data_z,verbose=False):
        """
        determines component for plotting
        """
        # convert r/theta/phi to x/y/z
        # x = r * sin(theta) * cos(phi)
        # y = r * sin(theta) * sin(phi)
        # z = r * cos(theta)
        #
        xcoord = r * np.sin(theta) * np.cos(phi)
        ycoord = r * np.sin(theta) * np.sin(phi)
        zcoord = r * np.cos(theta)

        # select plotting value
        if self.use_component == 1:
            # vertical Z-component
            if verbose: print("data: Z-component\n")
            # compute unit normal vector to the surface
            RRval = np.sqrt(xcoord**2 + ycoord**2 + zcoord**2)
            if np.any(RRval < 1.e-10):
                print("Error in unit normal vector")
                sys.exit(1)
            normal_x = xcoord / RRval
            normal_y = ycoord / RRval
            normal_z = zcoord / RRval

            # vertical component
            displn = data_x * normal_x + data_y * normal_y + data_z * normal_z

        elif self.use_component == 2:
            # N-component
            if verbose: print("data: N-component\n")
            # compute unit tangent vector to the surface (N-S)
            RRval = np.sqrt(xcoord**2 + ycoord**2 + zcoord**2)
            if np.any(RRval < 1.e-10):
                print("Error in unit normal vector")
                sys.exit(1)

            rhoval = np.sqrt(xcoord**2 + ycoord**2)
            # takes care for location at pole
            thetahat_x = np.where(rhoval < 1.e-10,0.0,(zcoord*xcoord) / (rhoval*RRval))
            thetahat_y = np.where(rhoval < 1.e-10,0.0,(zcoord*ycoord) / (rhoval*RRval))

            thetahat_z = - rhoval / RRval

            # N-component
            displn = - (data_x * thetahat_x + data_y * thetahat_y + data_z * thetahat_z)

        elif self.use_component == 3:
            # E-component
            if verbose: print("data: E-component\n")
            #compute unit tangent to the surface (E-W)
            rhoval = np.sqrt(xcoord**2 + ycoord**2)
            # takes care for location at pole
            phihat_x = np.where(rhoval < 1.e-10,0.0,-ycoord / rhoval)
            phihat_y = np.where(rhoval < 1.e-10,0.0,xcoord / rhoval)

            # E-component
            displn = data_x * phihat_x + data_y * phihat_y

        elif self.use_component == 4:
            # norm
            if verbose: print("data: norm-component\n")
            displn = np.sqrt(data_x**2 + data_y**2 + data_z**2)

        else:
            print("component type {} not supported yet".format(self.use_component))
            sys.exit(1)

        # plot data
        data = displn

        return data


    def read_moviedata(self,filename,verbose=False):
        """
        reads data from moviedata file
        """
        # input file directory
        name = os.path.basename(filename)
        dir = os.path.dirname(filename)

        if verbose:
            print("input file : ",filename)
            print("name       : ",name)
            print("directory  : ",dir)
            print("")


        # gets first coordinates
        filename_xyz = dir + "/" + "moviedata_xyz.bin"

        if verbose:
            print("reading in position file ... ",filename_xyz)
            print("")

        r,theta,phi = self.read_binary_moviedata_xyz_file(filename_xyz,verbose=False)

        # convert theta/phi to lat/lon
        #
        lat = 90.0 - theta * 180.0 / np.pi
        lon = phi * 180.0 / np.pi

        # longitude range [-180,180]
        lon = np.where(lon < -180.0,lon + 360.0,lon)
        lon = np.where(lon > 180.0,lon - 360.0,lon)

        if verbose:
            print("mesh dimensions: lat min/max = ",lat.min(),"/",lat.max())
            print("                 lon min/max = ",lon.min(),"/",lon.max())
            print("")

        # movie data values
        if verbose:
            print("reading in data file     ... ",filename)
            print("")

        # reads x/y/z data values
        data_x,data_y,data_z = self.read_binary_moviedata_file(filename,verbose=False)

        if verbose:
            print("moviedata: ",filename)
            print("  values x min/max = ",data_x.min(),"/",data_x.max())
            print("         y min/max = ",data_y.min(),"/",data_y.max())
            print("         z min/max = ",data_z.min(),"/",data_z.max())
            print("")

        # get data component to plot
        data = self.get_data_component(r,theta,phi,data_x,data_y,data_z,verbose=verbose)

        return lat,lon,data


    def interpolate_data_NNglobal(self,lat,lon,data,grid_resolution=100):
        """
        interpolates irregular lat/lon/data values on regular spherical grid using nearest neighbor value
        """
        # timing
        tic = time.perf_counter()

        print("interpolate: nearest neighbor")
        print("             grid resolution ",grid_resolution)

        # defines grid
        # global range
        xi = np.linspace(-180.0, 180.0, 2*grid_resolution) # x -> lon, y -> lat
        yi = np.linspace(-90.0, 90.0, grid_resolution)

        # grid data
        x = lon.reshape((lon.size))    # x -> lon, y -> lat
        y = lat.reshape((lat.size))
        z = data.reshape((data.size))

        # lengths
        nx = len(x)
        ny = len(y)
        nz = len(z)

        # checks to have same input lengths
        if nx != ny or nx != nz or ny != nz:
            print("error in interpolate: needs lat,lon,data with same lengths, but have nx/ny/nz = ",nx,ny,nz)
            sys.exit(1)

        # grid length
        nxi = len(xi)
        nyi = len(yi)

        #debug
        #print("interpolate: lengths nx/ny/nz = ",nx,ny,nz)
        #print("interpolate: lengths gridded nxi/nyi = ",nxi,nyi)

        # finds nearest neighbor value
        data_gridded = np.empty((nyi,nxi),dtype=data.dtype)  # for correct layout with plotting using (yi,xi) instead of (xi,yi)

        for j in range(0,nyi):
            # position
            pos_lat = yi[j]   # x -> lon, y -> lat

            # gets distance squared to position (lat0-lat)**2
            dist_y = (y - pos_lat)**2

            for i in range(0,nxi):
                # position
                pos_lon = xi[i]    # x -> lon, y -> lat

                # gets distance squared to position (lon0-lon)**2
                dist_x = (x - pos_lon)**2

                # gets distance squared to position
                #dist_x = np.power(np.subtract(pos_lat, x), 2)  # slower
                #dist_y = np.power(np.subtract(pos_lon, y), 2)
                #dist_x = (x - pos_lon)**2                      # faster
                #dist_y = (y - pos_lat)**2

                # gets distance squared to position (lon0-lon)**2 + (lat0-lat)**2
                dist_sq = dist_x + dist_y

                # minimum distance index
                ind = np.argmin(dist_sq)

                # stores closest value on grid
                data_gridded[j][i] = z[ind]

                #debug
                #print("interpolate: nearest ",dist_sq,"i,j",i,j,"pos",pos_lon,pos_lat,"ind",ind,"data val ",val)

        #debug
        #print("interpolate: data ",np.shape(data),data.dtype)
        #print("interpolate: data gridded",np.shape(data_gridded),data_gridded.dtype)

        # timing
        toc = time.perf_counter()
        print("interpolate: elapsed time {:0.4f} seconds\n".format(toc - tic))

        return xi,yi,data_gridded


    def interpolate_data(self,lat,lon,data,grid_resolution=100):
        """
        interpolates irregular lat/lon/data values on regular spherical grid
        """
        # timing
        tic = time.perf_counter()

        print("interpolate: linear interpolation")
        print("             grid resolution ",grid_resolution)

        # defines grid
        # range from data
        xi = np.linspace(lon.min(), lon.max(), 2*grid_resolution)  # x -> lon, y -> lat
        yi = np.linspace(lat.min(), lat.max(), grid_resolution)

        # grid data
        x = lon.reshape((lon.size))   # x -> lon, y -> lat
        y = lat.reshape((lat.size))
        z = data.reshape((data.size))

        # lengths
        nx = len(x)
        ny = len(y)
        nz = len(z)

        # checks to have same input lengths
        if nx != ny or nx != nz or ny != nz:
            print("error in interpolate: needs lat,lon,data with same lengths, but have nx/ny/nz = ",nx,ny,nz)
            sys.exit(1)

        # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
        use_matplotlib = True

        if use_matplotlib:
            # matplotlib
            import matplotlib.tri as tri
            print("interpolate: matplotlib LinearTriInterpolator")

            # extend x/y
            # border points at -180 will not be covered fully
            # extends arrays to wrap around lon-boundary by adding values
            mode = 0

            if mode == 0:
                # original range (no array extension)
                triang = tri.Triangulation(x, y)
                interpolator = tri.LinearTriInterpolator(triang, z)

            elif mode == 1:
                # wrap array at edges
                print("interpolate: wrap array at edges")
                x_extended = np.pad(x,(1,1),mode='wrap')  # lon
                y_extended = np.pad(y,(1,1),mode='wrap')  # lat
                z_extended = np.pad(z,(1,1),mode='wrap')  # vals

                # extended range
                triang = tri.Triangulation(x_extended, y_extended)
                interpolator = tri.LinearTriInterpolator(triang, z_extended)

            elif mode == 2:
                # manual array extension:
                #   lat [-90,90]   -> [-90,90]   same range
                #   lon [-180,180] -> [-540,540] by shifting and adding lon - 360, lon + 360
                print("interpolate: extend array at lon")
                x_extended = np.empty((nx + 2*nx),dtype=x.dtype)
                y_extended = np.empty((ny + 2*ny),dtype=y.dtype)
                z_extended = np.empty((nz + 2*nz),dtype=z.dtype)
                # shifts lon
                x_extended[0:nx] = x[:] - 360.0       # range [-540,-180]
                x_extended[nx:2*nx] = x[:]            # range [-180,180]
                x_extended[2*nx:3*nx] = x[:] + 360.0  # range [180,540]
                # keeps lat
                y_extended[0:ny] = y[:]
                y_extended[ny:2*ny] = y[:]
                y_extended[2*ny:3*ny] = y[:]
                # corresponding vals
                z_extended[0:nz] = z[:]
                z_extended[nz:2*nz] = z[:]
                z_extended[2*nz:3*nz] = z[:]

                # extended range
                triang = tri.Triangulation(x_extended, y_extended)
                interpolator = tri.LinearTriInterpolator(triang, z_extended)

            Xi, Yi = np.meshgrid(xi, yi)
            data_gridded = interpolator(Xi, Yi)

            # fill nan values with zero
            data_gridded = np.where(np.isnan(data_gridded),0.0,data_gridded)

        else:
            # scipy
            from scipy.interpolate import griddata
            print("interpolate: scipy griddata")
            data_gridded = griddata((x,y), z, (xi[None, :], yi[:, None]), method='linear',fill_value=0.0)

        # stats
        print("interpolate: min/max = {} / {}".format(data_gridded.min(),data_gridded.max()))

        # smoothing
        if self.use_smoothing:
            from scipy.ndimage import uniform_filter
            print("smoothing  : scipy uniform filter size = ",self.use_smoothing_kernel_size)
            data_gridded = uniform_filter(data_gridded, size=self.use_smoothing_kernel_size, mode='grid-wrap')
            print("smoothing  : min/max = {} / {}".format(data_gridded.min(),data_gridded.max()))

        # timing
        toc = time.perf_counter()
        print("interpolate: elapsed time {:0.4f} seconds\n".format(toc - tic))

        return xi,yi,data_gridded

