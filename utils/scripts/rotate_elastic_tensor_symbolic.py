#!/usr/bin/env python
#
# script to output symbolic formulas for rotating the elastic tensor
#
# daniel: the elastic tensor rotations in SPECFEM have been quite long expressions.
#         based on code comments in those routines, they have been derived using Mathematica. unfortunately, those original
#         scripts were not added to this packages and it became difficult to bug-fixing those.
#
#         this script here tries to replicate those Mathematica expressions, using python's sympy module.
#         very likely, the expressions were derived using tensor rotations based on Bond's matrix multiplications.
#         thus, starting with a rotation matrix, Bond's matrix will be constructed and used for
#         the tensor rotations, which assume Voigt's notation for a symmetric 4th-order elastic tensor.
#
#         the expressions should be consistent in both directions, from a local (radial) to global (Cartesian) reference
#         and vice versa.
#
#         for the implementation of these rotations in SPECFEM, we will use now Bond's matrix multiplication
#         since they are easier to double-check and fix in case needed. therefore, the symbolic expressions outputted
#         by this script are not needed for the code anymore. still, they have been used to double-check the routines.
#
import sys

from sympy import *

# gets arguments
if len(sys.argv) != 2:
    print("usage: ./rotate_elastic_tensor_symbolic.py  <type>")
    print(" with")
    print("    <type> == 1 to rotate from (lat/lon) to Cartesian (x/y/z) reference (local to global transform)")
    print(" or <type> == 2 to rotate from Cartesian (x/y/z) to (lat/lon)  reference (global to local transform)")
    sys.exit(1)
else:
    type = int(sys.argv[1])

## sympy
#init_printing(use_latex='mathjax')
init_printing(use_unicode=False, wrap_line=False, no_global=True)
#init_printing()

## rotation
# colat,lon
theta,phi = symbols(r'theta phi')

# single axis rotations
print("single axis rotation matrix for theta/phi:")
rotZ = rot_axis3(-phi)        # rotation back to pole
rotY = rot_axis2(-theta)
pprint(rotZ)
pprint(rotY)

# rotation matrix
if type == 1:
    # rotation matrix from spherical (colat/lon) to cartesian (x,y,z)
    print("rotation type:  spherical to Cartesian")
    rotmat = rotZ * rotY
else:
    # rotation from cartesian to spherical
    print("rotation type:  Cartesian to spherical")
    rotmat = rotZ * rotY
    rotmat = rotmat.transpose()

print("rotation matrix:")
pprint(rotmat)
print("")


# bond matrix for rotation of Voigt tensor
#
# see Bond matrix used in Auld, 1973
# or:
#  "OPTIMAL RECOVERY OF ELASTIC PROPERTIES FOR ANISOTROPIC MATERIALS THROUGH ULTRASONIC MEASUREMENTS"
#  Miao Sun, PhD thesis, B.S. University of Shanghai for Science & Technology, 1989
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.544.2565&rep=rep1&type=pdf
#
bond = zeros(6,6)
# First column
bond[0,0] = rotmat[0,0]*rotmat[0,0]         # a11**2
bond[1,0] = rotmat[1,0]*rotmat[1,0]         # a21**2
bond[2,0] = rotmat[2,0]*rotmat[2,0]         # a31**2
bond[3,0] = rotmat[1,0]*rotmat[2,0]         # a21*a31
bond[4,0] = rotmat[2,0]*rotmat[0,0]         # a31*a11
bond[5,0] = rotmat[0,0]*rotmat[1,0]         # a11*a21

# Second column
bond[0,1] = rotmat[0,1]*rotmat[0,1]         # a12**2
bond[1,1] = rotmat[1,1]*rotmat[1,1]         # a22**2
bond[2,1] = rotmat[2,1]*rotmat[2,1]         # a32**2
bond[3,1] = rotmat[1,1]*rotmat[2,1]         # a22*a32
bond[4,1] = rotmat[2,1]*rotmat[0,1]         # a32*a12
bond[5,1] = rotmat[0,1]*rotmat[1,1]         # a12*a22

# Third column
bond[0,2] = rotmat[0,2]*rotmat[0,2]         # a13**2
bond[1,2] = rotmat[1,2]*rotmat[1,2]         # a23**2
bond[2,2] = rotmat[2,2]*rotmat[2,2]         # a33**2
bond[3,2] = rotmat[1,2]*rotmat[2,2]         # a23*a33
bond[4,2] = rotmat[2,2]*rotmat[0,2]         # a33*a13
bond[5,2] = rotmat[0,2]*rotmat[1,2]         # a13*a23

# Fourth column
bond[0,3] = 2.0 * rotmat[0,1]*rotmat[0,2]   # 2 * a12*a13
bond[1,3] = 2.0 * rotmat[1,1]*rotmat[1,2]   # 2 * a22*a23
bond[2,3] = 2.0 * rotmat[2,1]*rotmat[2,2]   # 2 * a32*a33
bond[3,3] = rotmat[1,1]*rotmat[2,2] + rotmat[1,2]*rotmat[2,1]   # a22*a33 + a23*a32
bond[4,3] = rotmat[0,1]*rotmat[2,2] + rotmat[0,2]*rotmat[2,1]   # a12*a33 + a13*a32
bond[5,3] = rotmat[0,1]*rotmat[1,2] + rotmat[0,2]*rotmat[1,1]   # a12*a23 + a13*a22

# Fifth column
bond[0,4] = 2.0 * rotmat[0,2]*rotmat[0,0]
bond[1,4] = 2.0 * rotmat[1,2]*rotmat[1,0]
bond[2,4] = 2.0 * rotmat[2,2]*rotmat[2,0]
bond[3,4] = rotmat[1,0]*rotmat[2,2] + rotmat[1,2]*rotmat[2,0]
bond[4,4] = rotmat[0,2]*rotmat[2,0] + rotmat[0,0]*rotmat[2,2]
bond[5,4] = rotmat[0,2]*rotmat[1,0] + rotmat[0,0]*rotmat[1,2]

# Sixth column
bond[0,5] = 2.0 * rotmat[0,0]*rotmat[0,1]
bond[1,5] = 2.0 * rotmat[1,0]*rotmat[1,1]
bond[2,5] = 2.0 * rotmat[2,0]*rotmat[2,1]
bond[3,5] = rotmat[1,1]*rotmat[2,0] + rotmat[1,0]*rotmat[2,1]
bond[4,5] = rotmat[0,0]*rotmat[2,1] + rotmat[0,1]*rotmat[2,0]
bond[5,5] = rotmat[0,0]*rotmat[1,1] + rotmat[0,1]*rotmat[1,0]

print("bond matrix:")
pprint(bond)
print("")

# Get transpose of bond
# this will lead to an "expand.." error later on, doing it manually as a work-around...
#bond_t = Transpose(bond)
bond_t = zeros(6,6)
for i in range(0,6):
    for j in range(0,6):
        bond_t[i,j] = bond[j,i]

## elastic tensor Cij
if type == 1:
    # spherical (local) to Cartesian (global)
    # uses coefficients symbols dij for local
    c11,c12,c13,c14,c15,c16 = symbols(r'd11 d12 d13 d14 d15 d16')
    c22,c23,c24,c25,c26 = symbols(r'd22 d23 d24 d25 d26')
    c33,c34,c35,c36 = symbols(r'd33 d34 d35 d36')
    c44,c45,c46 = symbols(r'd44 d45 d46')
    c55,c56 = symbols(r'd55 d56')
    c66 = symbols(r'd66')
else:
    # Cartesian (global) to spherical (local)
    # uses coefficients symbols cij for global
    c11,c12,c13,c14,c15,c16 = symbols(r'c11 c12 c13 c14 c15 c16')
    c22,c23,c24,c25,c26 = symbols(r'c22 c23 c24 c25 c26')
    c33,c34,c35,c36 = symbols(r'c33 c34 c35 c36')
    c44,c45,c46 = symbols(r'c44 c45 c46')
    c55,c56 = symbols(r'c55 c56')
    c66 = symbols(r'c66')

# elastic tensor (Voigt notation)
cij = Matrix(([c11,c12,c13,c14,c15,c16],[c12,c22,c23,c24,c25,c26],[c13,c23,c33,c34,c35,c36],
              [c14,c24,c34,c44,c45,c46],[c15,c25,c35,c45,c55,c56],[c16,c26,c36,c46,c56,c66] ))

print("elastic tensor:")
pprint(cij)
print("")
print("symmetric:",cij.is_symmetric())
print("")

## rotates Cij
#(from spherical to cartesian reference system)
# First compute C M^t
tensor_tmp = zeros(6,6)
for j in range(0,6):
     for k in range(0,6):
          for i in range(0,6):
              tensor_tmp[i,j] = tensor_tmp[i,j] + cij[i,k] * bond_t[k,j]
# Second compute M * (C M^t)
tensor_r = zeros(6,6)
for j in range(0,6):
     for k in range(0,6):
          for i in range(0,6):
              tensor_r[i,j] = tensor_r[i,j] + bond[i,k] * tensor_tmp[k,j]


# print rotated tensor
print("rotated tensor:")
for i in range(0,6):
    for j in range(i,6):
        if type == 1:
            # local (spherical) to global (Cartesian)
            print("element:",i+1,j+1,'c{}{}'.format(i+1,j+1))
        else:
            # global (Cartesian) to local (spherical)
            print("element:",i+1,j+1,'d{}{}'.format(i+1,j+1))
        print("")
        ## pure element
        #print(tensor_r[i,j])
        #print("")
        ## sympy simplication (e.g. combines cosine/sine expressions ... not really needed)
        #print(simplify(tensor_r[i,j]))
        ## fortran code ... just in case, but formatting looks ugly
        #print(fcode(tensor_r[i,j], standard=95))
        ## collect
        # (similar to what Dimitri mentioned about the original routines based on Mathematica Collect[Expand(Cij),{dij}]
        print(collect(expand(tensor_r[i,j]),{c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66}))
        # alternative collect
        #print(collect(tensor_r[i,j],{c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66}))
        print("")

print("symmetric:",tensor_r.is_symmetric())
print("")

