#!/usr/bin/env python
#
# creates fictitious multiple source file
from __future__ import print_function
import sys

# default values
header0 = "PDEW 2008  1  6  5 14 23.70  37.2200   22.6900  75.0 6.1 6.2 SOUTHERN GREECE"
name0 = "200801060514A"

tshift0 = 0.000
hdur0 = 1.0 # original: 3.0000
lat0 = 36.9800
lon0 = 22.8700
dep0 = 92.3900

Mrr0 =  7.740000e+24
Mtt0 = -1.830000e+25
Mpp0 =  1.060000e+25
Mrt0 =  1.290000e+25
Mrp0 = -8.450000e+24
Mtp0 = -4.430000e+24

def print_cmtlist(cmtlist):
    """
    print CMT solution
    """

    length = len(cmtlist)

    for i in range(0,length):
        header = cmtlist[i][0]
        name = "cmt%d" % (i+1)  #cmtlist[i][1]
        tshift = cmtlist[i][2]
        hdur = cmtlist[i][3]
        lat = cmtlist[i][4]
        lon = cmtlist[i][5]
        dep = cmtlist[i][6]
        Mrr = cmtlist[i][7]
        Mtt = cmtlist[i][8]
        Mpp = cmtlist[i][9]
        Mrt = cmtlist[i][10]
        Mrp = cmtlist[i][11]
        Mtp = cmtlist[i][12]

        print(header)
        print("event name:     ",name)
        print("time shift:     ",str(tshift))
        print("half duration:  ",str(hdur))
        print("latitude:       ",str(lat))
        print("longitude:      ",str(lon))
        print("depth:          ",str(dep))
        print("Mrr:       ",str(Mrr))
        print("Mtt:       ",str(Mtt))
        print("Mpp:       ",str(Mpp))
        print("Mrt:       ",str(Mrt))
        print("Mrp:       ",str(Mrp))
        print("Mtp:       ",str(Mtp))


def create_cmts(lat1,lon1,cmtlist):
    """
    creates a line of CMTs at a given lat/lon location down to CMB
    """
    global header0,name0,tshift0,hdur0,Mrr0,Mtt0,Mpp0,Mrt0,Mrp0,Mtp0

    # margin at surface and CMB
    margin = 10.0 # km

    # depth range of CMTs (between surface and CMB, with 10km margin)
    dr = ((6371.0 - margin) - (3480.0 + margin))

    # creates vertical line of CMTs
    num_cmts_per_line = 2
    for j in range(0,num_cmts_per_line):
        # CMTSOLUTION entry
        header = header0
        name = name0
        tshift = j * 0.2  # slightly time shifted from top to bottom
        hdur = hdur0
        lat = lat1
        lon = lon1
        dep = 10.0 + j * dr/(num_cmts_per_line-1)
        Mrr = Mrr0
        Mtt = Mtt0
        Mpp = Mpp0
        Mrt = Mrt0
        Mrp = Mrp0
        Mtp = Mtp0
        cmt = [ header, name, tshift, hdur, lat, lon, dep, Mrr, Mtt, Mpp, Mrt, Mrp, Mtp ]

        # adds to list
        cmtlist.append(cmt)


def read_STATIONS(filename,cmtlist):
    """
    reads in STATIONS file
    """

    # note: sta_info = np.loadtxt(filename) will not work with string and numbers on same line
    f = open(filename, "r")
    lines = f.readlines()
    for i in range(0,len(lines)):
        a = lines[i].split()
        lat = float(a[2])
        lon = float(a[3])
        #print("station: ",a[0],lat1,lon1)

        # create CMT entries
        create_cmts(lat,lon,cmtlist)

def create_CMTSOLUTION(sta_file):
    """
    creates multiple CMT solutions
    """
    cmtlist = []

    # reads in stations and creates
    read_STATIONS(sta_file,cmtlist)

    # output
    print_cmtlist(cmtlist)


def usage():
    print("usage: ./create_CMTSOLUTION_file.py STATIONS")
    print("   where")
    print("       STATIONS - station location file, e.g. DATA/STATIONS")


if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 2:
        usage()
        sys.exit(1)
    else:
        sta_file = sys.argv[1]

    create_CMTSOLUTION(sta_file)

