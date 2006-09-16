!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 6
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!       (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! create mesh quality data for the slice, to be recombined in postprocessing

  subroutine write_AVS_DX_mesh_quality_data(prname,nspec,xstore,ystore,zstore)

  implicit none

  include "constants.h"

  integer nspec
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  integer ispec

  double precision xelm(8)
  double precision yelm(8)
  double precision zelm(8)

  integer iface,icorner,jcorner

  double precision vectorA_x,vectorA_y,vectorA_z
  double precision vectorB_x,vectorB_y,vectorB_z
  double precision norm_A,norm_B,angle_vectors
  double precision distmin,distmax,dist,dist1,dist2,dist3,dist4
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio

! topology of faces of cube for skewness
  integer faces_topo(6,6)

! list of corners defining the edges
  integer iedge
  integer list_corner(2,12)

! processor identification
  character(len=150) prname

! data and element files are identical to AVS_DXpoints.txt and AVS_DXelements.txt
! created in regular AVS or DX routine, therefore not created again here

! define topology of faces of cube for skewness

! face 1
  faces_topo(1,1) = 1
  faces_topo(1,2) = 2
  faces_topo(1,3) = 6
  faces_topo(1,4) = 5

! face 2
  faces_topo(2,1) = 2
  faces_topo(2,2) = 3
  faces_topo(2,3) = 7
  faces_topo(2,4) = 6

! face 3
  faces_topo(3,1) = 4
  faces_topo(3,2) = 3
  faces_topo(3,3) = 7
  faces_topo(3,4) = 8

! face 4
  faces_topo(4,1) = 1
  faces_topo(4,2) = 5
  faces_topo(4,3) = 8
  faces_topo(4,4) = 4

! face 5
  faces_topo(5,1) = 1
  faces_topo(5,2) = 2
  faces_topo(5,3) = 3
  faces_topo(5,4) = 4

! face 6
  faces_topo(6,1) = 5
  faces_topo(6,2) = 6
  faces_topo(6,3) = 7
  faces_topo(6,4) = 8

! define wraparound for angles for skewness calculation
  faces_topo(:,5) = faces_topo(:,1)
  faces_topo(:,6) = faces_topo(:,2)

! list of corners defining the edges

  list_corner(1,1) = 1
  list_corner(2,1) = 2

  list_corner(1,2) = 2
  list_corner(2,2) = 3

  list_corner(1,3) = 3
  list_corner(2,3) = 4

  list_corner(1,4) = 4
  list_corner(2,4) = 1

  list_corner(1,5) = 5
  list_corner(2,5) = 6

  list_corner(1,6) = 6
  list_corner(2,6) = 7

  list_corner(1,7) = 7
  list_corner(2,7) = 8

  list_corner(1,8) = 8
  list_corner(2,8) = 5

  list_corner(1,9) = 1
  list_corner(2,9) = 5

  list_corner(1,10) = 2
  list_corner(2,10) = 6

  list_corner(1,11) = 3
  list_corner(2,11) = 7

  list_corner(1,12) = 4
  list_corner(2,12) = 8

! writing mesh quality data for each element
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='unknown')

! number of elements in AVS or DX file
  write(10,*) nspec

! output global AVS or DX elements
  do ispec=1,nspec

! define the coordinates of the 8 corners of the element
     xelm(1) = xstore(1,1,1,ispec)
     yelm(1) = ystore(1,1,1,ispec)
     zelm(1) = zstore(1,1,1,ispec)

     xelm(2) = xstore(NGLLX,1,1,ispec)
     yelm(2) = ystore(NGLLX,1,1,ispec)
     zelm(2) = zstore(NGLLX,1,1,ispec)

     xelm(3) = xstore(NGLLX,NGLLY,1,ispec)
     yelm(3) = ystore(NGLLX,NGLLY,1,ispec)
     zelm(3) = zstore(NGLLX,NGLLY,1,ispec)

     xelm(4) = xstore(1,NGLLY,1,ispec)
     yelm(4) = ystore(1,NGLLY,1,ispec)
     zelm(4) = zstore(1,NGLLY,1,ispec)

     xelm(5) = xstore(1,1,NGLLZ,ispec)
     yelm(5) = ystore(1,1,NGLLZ,ispec)
     zelm(5) = zstore(1,1,NGLLZ,ispec)

     xelm(6) = xstore(NGLLX,1,NGLLZ,ispec)
     yelm(6) = ystore(NGLLX,1,NGLLZ,ispec)
     zelm(6) = zstore(NGLLX,1,NGLLZ,ispec)

     xelm(7) = xstore(NGLLX,NGLLY,NGLLZ,ispec)
     yelm(7) = ystore(NGLLX,NGLLY,NGLLZ,ispec)
     zelm(7) = zstore(NGLLX,NGLLY,NGLLZ,ispec)

     xelm(8) = xstore(1,NGLLY,NGLLZ,ispec)
     yelm(8) = ystore(1,NGLLY,NGLLZ,ispec)
     zelm(8) = zstore(1,NGLLY,NGLLZ,ispec)

! compute equiangle skewness (as defined in Fluent/Gambit manual)
     equiangle_skewness = - HUGEVAL
     do iface = 1,6
       do icorner = 1,4

! first vector of angle
         vectorA_x = xelm(faces_topo(iface,icorner)) - xelm(faces_topo(iface,icorner+1))
         vectorA_y = yelm(faces_topo(iface,icorner)) - yelm(faces_topo(iface,icorner+1))
         vectorA_z = zelm(faces_topo(iface,icorner)) - zelm(faces_topo(iface,icorner+1))

! second vector of angle
         vectorB_x = xelm(faces_topo(iface,icorner+2)) - xelm(faces_topo(iface,icorner+1))
         vectorB_y = yelm(faces_topo(iface,icorner+2)) - yelm(faces_topo(iface,icorner+1))
         vectorB_z = zelm(faces_topo(iface,icorner+2)) - zelm(faces_topo(iface,icorner+1))

! norm of vectors A and B
         norm_A = dsqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
         norm_B = dsqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

! angle formed by the two vectors
         angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

! compute equiangle skewness
         equiangle_skewness = dmax1(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

       enddo
     enddo

! compute edge aspect ratio using the 8 corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     do iedge = 1,12
       icorner = list_corner(1,iedge)
       jcorner = list_corner(2,iedge)
       dist = dsqrt((xelm(jcorner) - xelm(icorner))**2 + (yelm(jcorner) - yelm(icorner))**2 + (zelm(jcorner) - zelm(icorner))**2)
       distmin = dmin1(distmin,dist)
       distmax = dmax1(distmax,dist)
     enddo
     edge_aspect_ratio = distmax / distmin

! compute diagonal aspect ratio
     dist1 = dsqrt((xelm(1) - xelm(7))**2 + (yelm(1) - yelm(7))**2 + (zelm(1) - zelm(7))**2)
     dist2 = dsqrt((xelm(2) - xelm(8))**2 + (yelm(2) - yelm(8))**2 + (zelm(2) - zelm(8))**2)
     dist3 = dsqrt((xelm(3) - xelm(5))**2 + (yelm(3) - yelm(5))**2 + (zelm(3) - zelm(5))**2)
     dist4 = dsqrt((xelm(4) - xelm(6))**2 + (yelm(4) - yelm(6))**2 + (zelm(4) - zelm(6))**2)
     distmin = dmin1(distmin,dist1,dist2,dist3,dist4)
     distmax = dmax1(distmax,dist1,dist2,dist3,dist4)
     diagonal_aspect_ratio = distmax / distmin

! write mesh quality information for each element
     write(10,*) ispec,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio


  enddo

  close(10)

  end subroutine write_AVS_DX_mesh_quality_data

