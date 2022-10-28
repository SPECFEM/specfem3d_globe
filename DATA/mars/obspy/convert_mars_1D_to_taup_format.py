#!/usr/bin/env python
#
# script to convert mars_1D.dat to taup-format mars_1D.nd
#
import sys
import numpy as np


# reads in mars_1D.dat
filename = "../mars_1D.dat"

print("reading file ",filename)
print("")

# model values start at surface, increase down to center of inner core
# format: #r(km) #vp(km/s) #vs(km/s) #density(kg/cm^3) #Qp  #Qs
dat = np.loadtxt(filename)

# taup needs depths for surface, moho/mantle, outer-core, inner-core
index_moho = 0
index_cmb = 0
index_icb = 0

# surface radius
r_surface = dat[0][0]

# loops over depths to determine boundaries
iinterface = 0
r_prev = r_surface
for i in range(1,len(dat)):
    r = dat[i][0]

    # discontinuities will have two entries for above/below interface, both at same depth
    # tolerance to detect identical depth values
    TOL = 1.e-9

    # interfaces in mars_1D.dat are: upper crust, middle crust, Moho, CMB, ICB
    # checks if depths increase
    if abs(r_prev - r) < TOL:
        # interface
        iinterface += 1
        if iinterface == 1:
            # upper crust
            print("  skipping upper  crust interface at ",i)
            continue
        elif iinterface == 2:
            # middle crust
            print("  skipping middle crust interface at ",i)
            continue
        elif iinterface == 3:
            index_moho = i
            print("  found moho at {:3d} - radius {}   depth {}".format(i,dat[i][0],r_surface - dat[i][0]))
        elif iinterface == 4:
            index_cmb = i
            print("  found cmb  at {:3d} - radius {}   depth {}".format(i,dat[i][0],r_surface - dat[i][0]))
        elif iinterface == 5:
            index_icb = i
            print("  found icb  at {:3d} - radius {}   depth {}".format(i,dat[i][0],r_surface - dat[i][0]))
        else:
            print("warning: more interfaces detected than assumed (Moho,CMB,ICB) at line ",i)
            print("  previous line: ",dat[i][0],dat[i][1],dat[i][2],dat[i][3],dat[i][4],dat[i][5])
            print("  current line : ",dat[i][0],dat[i][1],dat[i][2],dat[i][3],dat[i][4],dat[i][5])
            print("  interfaces detected: ",iinterface)
            print("  index moho = {} / index cmb = {} / index icb = {}".format(index_moho,index_cmb,index_icb))
            sys.exit(1)

    # checks if depths increase
    if abs(r_prev - r) > TOL and r_prev - r < 0.0:
        print("warning: depth at line ",i," not increasing")
        print("  previous line: ",dat[i-1][0],dat[i-1][1],dat[i-1][2],dat[i-1][3],dat[i-1][4],dat[i-1][5])
        print("  current line : ",dat[i][0],dat[i][1],dat[i][2],dat[i][3],dat[i][4],dat[i][5])

    # update previous
    r_prev = r

print("")
print("model data: ")
print("  radius surface: ",r_surface)
print("  interface index: moho = ",index_moho)
print("                   cmb  = ",index_cmb)
print("                   icb  = ",index_icb)
print("")

# creates taup-format file
filename_taup = "./mars_1D.nd"

with open(filename_taup,'w') as f:
    # no header

    for i in range(len(dat)):
        # taup-format: #depth #rho #vp #vs
        depth = r_surface - dat[i][0]
        vp = dat[i][1]
        vs = dat[i][2]
        rho = dat[i][3]

        # interface info
        if i == index_moho:
            f.write("mantle\n")
        elif i == index_cmb:
            f.write("outer-core\n")
        elif i == index_icb:
            f.write("inner-core\n")

        # data line
        f.write("{} {} {} {}\n".format(depth,vp,vs,rho))


print("converted to: ",filename_taup)
print("")

