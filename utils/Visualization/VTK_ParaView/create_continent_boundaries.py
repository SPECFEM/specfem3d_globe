#!/usr/bin/env python
#
#
import os
import sys
import subprocess
import datetime

#import pygmt
import numpy as np

######################################################################
# USER Parameters

# flag to use either high-resolution or low-resolution coast lines
use_high_resolution = False

# ellipticity correction
use_ellipticity = True

# spherical Earth radius (SPECFEM3D_GLOBE uses normalized radius)
r_earth = 1.0

######################################################################

fix_region = None    # user specified region code

def geo2xyz(lon,lat):
    global r_earth
    global use_ellipticity

    # perl example
    # convert geographic latitude to geocentric colatitude and convert to radians
    # $pi = 3.14159265;
    # $theta = $pi/2. - atan2(0.99329534 * tan($latitude * $pi / 180.),1) ;
    # $phi = $longitude * $pi / 180. ;
    # compute the Cartesian position of the receiver (ignore ellipticity for AVS)
    # assume a sphere of radius one
    # $r_target = 1. ;
    ## DK DK make the radius a little bit bigger to make sure it is
    ## DK DK correctly superimposed to the mesh in final AVS figure
    # $r_target = 1.015 ;
    # $x_target = $r_target*sin($theta)*cos($phi) ;
    # $y_target = $r_target*sin($theta)*sin($phi) ;
    # $z_target = $r_target*cos($theta) ;

    # considers ellipticity to convert latitude to colatitude
    if use_ellipticity:
        # corrects latitude by ellipticity factor
        # using same constants as in constants.h
        EARTH_FLATTENING_F = 1.0 / 299.80
        ONE_MINUS_F_SQUARED = (1.0 - EARTH_FLATTENING_F)**2
        # colatitude/longitude in rad
        theta = np.pi/2.0 - np.arctan2(ONE_MINUS_F_SQUARED * np.tan(lat * np.pi / 180.0),1)
        phi = lon * np.pi / 180.0
    else:
        # geocentric and geodetic latitude is the same for spherical Earth
        # colatitude/longitude in rad
        theta = (90.0 - lat) * np.pi / 180.0
        phi = lon * np.pi / 180.0

    # compute the Cartesian position of the receiver (ignore ellipticity for AVS)
    # assume a sphere of radius one
    # (or make the radius a little bit bigger to make sure it is correctly superimposed to the mesh in final AVS figure)
    r_target = r_earth
    x_target = r_target * np.sin(theta) * np.cos(phi)
    y_target = r_target * np.sin(theta) * np.sin(phi)
    z_target = r_target * np.cos(theta)

    return x_target,y_target,z_target

def check_status(status):
    if status != 0:
        print("error: status returned ",status)
        sys.exit(status)
    return

def create_AVS_file():
    global use_ellipticity
    global use_high_resolution
    global gmt_region

    print("*******************************")
    print("creating AVS border file ...")
    print("*******************************")

    # current directory
    dir = os.getcwd()
    #print("current directory:",dir)
    #print("")

    # resolution
    if use_high_resolution:
        D = '-Dh'
    else:
        D = '-Dl'

    # region
    if not fix_region is None:
        gmt_region = f"-R{fix_region}"
    else:
        # global region
        gmt_region = '-Rg'

    # GMT segment file
    name = "map_segment.dat"
    cmd = f"gmt pscoast {gmt_region}  {D} -W -M > {name} ;"
    print("  cmd: ",cmd)
    status = subprocess.call(cmd, shell=True)
    check_status(status)
    print("  GMT segment file plotted in file: ",name)

    # note: GMT segment file has format
    # > Shore Bin # .., Level ..
    #   lon1 lat1
    #   lon2 lat2
    # > Shore Bin # .., Level ..
    #   ..

    print("  Getting boundaries from file %s ..." % name)

    # reads gmt boundary file
    with open(name,'r') as f:
        content = f.readlines()

    if len(content) == 0:
        print("")
        print("  INFO: no boundaries in file")
        print("")
        return

    # counts segments and points
    numsegments = 0
    numpoints = 0
    for line in content:
        if ">" in line.strip():
            # new segment
            numsegments += 1
        else:
            # point
            numpoints += 1

    print("  There are %i contours" % numsegments)
    print("  There are %i data points" % numpoints)

    # read the GMT file to get the number of individual line segments
    currentelem = 0
    previous_was_comment = 1
    for line in content:
        line = line.strip()
        #debug
        #print("currentelem %i %i %s" % (currentelem,previous_was_comment,line))
        # skip comment lines
        if line[0:1] == "#": continue
        # get line marker (comment in file)
        if ">" in line:
            previous_was_comment = 1
        else:
            if previous_was_comment == 0: currentelem += 1
            previous_was_comment = 0

    num_individual_lines = currentelem
    print("  There are %i individual line segments" % num_individual_lines)

    print("")
    print("converting to .inp format:")
    print("  normalized radius: ",r_earth)
    if use_ellipticity:
        print("  using elliptical Earth (latitude correction)")
    else:
        print("  using spherical Earth")
    print("")

    # output filename
    if use_ellipticity:
        avsfile = "AVS_boundaries_elliptical.inp"
    else:
        avsfile = "AVS_boundaries_spherical.inp"

    with open(avsfile,'w') as f:
        # write header for AVS (with point data)
        f.write("%i %i 1 0 0\n" % (numpoints,num_individual_lines))

        # read the GMT file to get the points
        currentpoint = 0
        for line in content:
            line = line.strip()
            # skip comment lines
            if line[0:1] == "#": continue

            #   get point only if line is not a comment
            if ">" not in line:
                currentpoint += 1
                elem = line.split()

                ## global lon/lat coordinates
                # longitude is the number before the white space
                lon = float(elem[0])
                # latitude is the number after the white space
                lat = float(elem[1])

                # converts to geocentric x/y/z position
                x,y,z = geo2xyz(lon,lat)

                # location
                x_target = x
                y_target = y
                z_target = z

                f.write("%i %f %f %f\n" % (currentpoint,x_target,y_target,z_target))

        # read the GMT file to get the lines
        currentline = 0
        currentelem = 0
        currentpoint = 0
        previous_was_comment = 1
        for line in content:
            line = line.strip()
            # skip comment lines
            if line[0:1] == "#": continue
            #   get line marker (comment in file)
            if ">" in line:
                # check if previous was line was also a segment
                # for example: there can be empty segment lines
                #  > Shore Bin # 4748, Level 1
                #  > Shore Bin # 4748, Level 1
                #  > Shore Bin # 4748, Level 1
                #  136.117036698   36.2541237507
                #  136.121248188   36.2533302815
                #  ..
                if currentline > 0 and previous_was_comment :
                    continue
                else:
                    currentline += 1
                    currentpoint += 1
                    previous_was_comment = 1
                #print("processing contour %i named %s" % (currentline,line))
            else:
                if previous_was_comment == 0:
                    previouspoint = currentpoint
                    currentelem  +=1
                    currentpoint  +=1
                    # new line
                    f.write("%i %i line %i %i\n" % (currentelem,currentline,previouspoint,currentpoint))
                previous_was_comment = 0

        # dummy variable names
        f.write(" 1 1\n")
        f.write(" Zcoord, meters\n")
        # create data values for the points
        for currentpoint in range(1,numpoints+1):
            f.write("%i 255.\n" % (currentpoint))

    # check
    if numpoints != currentpoint:
        print("  WARNING:")
        print("    possible format corruption: total number of points ",numpoints," should match last line point id ",currentpoint)
        print("")

    print("  see file: %s" % avsfile)
    print("")
    return


def usage():
    print("usage: ./create_continent_boundaries.py [--hires] [--lowres] [--ellip] [--sphere] [--region=R]")
    print("  with")
    print("     --hires         - high resolution global coastlines")
    print("     --lowres        - low resolution global coastlines")
    print("     --ellip         - (optional) for an elliptical Earth (default)")
    print("     --sphere        - (optional) for a spherical Earth")
    print("     --region=R      - (optional) use a fixed region specifier R (e.g. 'lonmin/lonmax/latmin/latmax')")
    
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
        if "--help" in arg:
            usage()
        elif "--hires" in arg:
            use_high_resolution = True
        elif "--lowres" in arg:
            use_high_resolution = True
        elif "--ellip" in arg:
            use_ellipticity = True
        elif "--sphere" in arg:
            use_ellipticity = False
        elif "--region" in arg:
            fix_region = arg.split('=')[1]
        elif i > 1:
            print("argument not recognized: ",arg)
            sys.exit(1)

    # logging
    cmd = " ".join(sys.argv)
    filename = './create_continent_boundaries.log'
    with open(filename, 'a') as f:
      print("command call --- " + str(datetime.datetime.now()),file=f)
      print(cmd,file=f)
      print("command logged to file: " + filename)

    # main routine
    create_AVS_file()

