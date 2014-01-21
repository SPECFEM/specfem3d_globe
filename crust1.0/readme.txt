The files distributed with this README are the model CRUST1.0 and three
Fortran 77 routines to get a 1D profile at a point (give latitude and 
longitude) (getCN1point.f)  or a set of maps (getCN1maps.f)  or a set 
of xyz datai (getCN1xyz.f)..

CRUST1.0 is a 8 layer model defined as 1x1 degree cells. The cells
average crustal structure over a 1x1 degree cells, i.e. the grid
associated with the model is defined as the center of these cells.
If a model is therefore requested for a certain cell, please use
the midpoint of this cell. For example, the model in the cell
5 to 6 deg latitude and 150 to 151 deg longitude should be inquired at
5.5 deg latitude and 150.5 deg longitude.

Prototype:
=========
CRUST1.0 is a proto type model. Please test the model against your data
and provide feedback to glaske-at-ucsd.edu. We plan to release an update,
CRUST1.1 in 2014.

Initial release:
===============
15 July 2013: this is the initial release of the essential files. As
described on the website, http://igppweb.ucsd.edu/~gabi/crust1.html,
the structure in the crystalline crust is defined using a crustal type
assignment. The crustal type file can be obtained upon request.

The 8 crustal layers:
====================
1) water             
2) ice               
3) upper sediments   (VP, VS, rho not defined in all cells) 
4) middle sediments  "
5) lower sediments   "
6) upper crystalline crust
7) middle crystalline crust
8) lower crystalline crust

a ninth layer gives V_Pn, V_Sn and rho below the Moho. The values
are associated with LLNL model G3Cv3 on continents and a thermal
model in the oceans.

Files:
=====
model:
crust1.bnds: boundary topography.
crust1.vp:   vp
crust1.vs:   vs
crust1.rho:  rho 

Format:
======
The model is defined from 89.5 to -89.5 deg latitude and -179.5 to 179.5 deg
longitude. Longitudes are the inner loop, i.e. all longitudes are stored
for each latitude, then the next latitude is given. The model starts at 
89.5 N and 179.5 W.

Each line has the nine values for one parameter, in one cell, i.e.
there are 360x180=64800 lines in each file. 

Current output files using getCN1maps.f:
=======================================
map-bd[x], where x goes from 1 to 9: 
 1) top of water
 2) bottom of water
 3) bottom of ice
 4) bottom of sediments 1
 5) bottom of sediments 2
 6) bottom of sediments 3
 7) bottom of cryst. crust 1
 8) bottom of cryst. crust 2
 9) bottom of cryst. crust 3 = Moho (depth to Moho, not crustal thickness!)
map-vp[x], x=1 to 8 is vp incrustal layers; x=9 is VPn
map-vs[x], same for vs
map=ro[x], same for density

file sedthk: sediment thickness
file crsthk: crustal thickness (without water)

same for getCN1xyz.f, except that each cell has its own line, 
with longitude, latitude, value.

the maps go from 89.5 to -89.5 latitude and from -179.5 to 179.5 longitude
