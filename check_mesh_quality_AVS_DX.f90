!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! combine mesh quality data files to check the mesh
! displays statistics on mesh quality
! and creates an AVS or OpenDX file showing a given range of elements

  program mesh_quality_AVS_DX

  implicit none

  include "constants.h"

  integer iproc,nspec,npoin
  integer ispec
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer ipoin,numpoin,iglobpointoffset,ntotpoin,ntotspec
  integer numelem,iglobelemoffset
  integer iformat,idoubling,iregion_code
  integer ntotpoinAVS_DX,ntotspecAVS_DX

  double precision xval,yval,zval

! processor identification
  character(len=150) prname

! for chunk numbering
  integer iproc_read,ichunk,idummy1,idummy2
  integer, dimension(:), allocatable :: ichunk_slice

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_AVS_DX_MESH_FILES,ATTENUATION,IASPEI, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

! for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
               nglob_AB,nglob_AC,nglob_BC

! for quality of mesh
  logical, dimension(:), allocatable :: mask_ibool
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision skewness_AVS_DX_min,skewness_AVS_DX_max

! for histogram
  integer, parameter :: NCLASS = 20
  integer classes_skewness(0:NCLASS-1)
  integer iclass
  double precision current_percent,total_percent

  integer proc_p1,proc_p2

  logical USE_OPENDX

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all mesh quality files for slices'
  print *

  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if(iformat<1 .or. iformat>2) stop 'exiting...'
  if(iformat == 1) then
    USE_OPENDX = .true.
  else
    USE_OPENDX = .false.
  endif

! read range of skewness used for elements
  print *,'enter minimum skewness for AVS or DX (between 0. and 1.):'
  read(5,*) skewness_AVS_DX_min
  if(skewness_AVS_DX_min < 0.d0) skewness_AVS_DX_min = 0.d0
  if(skewness_AVS_DX_min > 0.99999d0) skewness_AVS_DX_min = 0.99999d0

  print *,'enter maximum skewness for AVS or DX (between 0. and 1.):'
  read(5,*) skewness_AVS_DX_max
  if(skewness_AVS_DX_max < 0.d0) skewness_AVS_DX_max = 0.d0
  if(skewness_AVS_DX_max > 0.99999d0) skewness_AVS_DX_max = 0.99999d0

  if(skewness_AVS_DX_min > skewness_AVS_DX_max) stop 'incorrect skewness range'

  print *
  print *,'reading parameter file'
  print *

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_AVS_DX_MESH_FILES, &
          ATTENUATION,IASPEI,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL)

  if(.not. SAVE_AVS_DX_MESH_FILES) stop 'AVS or DX files were not saved by the mesher'

! compute other parameters based upon values read
  call compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
      NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
      NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NCHUNKS,INCLUDE_CENTRAL_CUBE)

! user can specify a range of processors here, enter 0 and 0 for all procs
  print *
  print *,'enter first proc (proc numbers start at 0) = '
  read(5,*) proc_p1
  if(proc_p1 < 0) proc_p1 = 0
  if(proc_p1 > NPROCTOT-1) proc_p1 = NPROCTOT-1

  print *,'enter last proc (enter -1 for all procs) = '
  read(5,*) proc_p2
  if(proc_p2 == -1) proc_p2 = NPROCTOT-1
  if(proc_p2 < 0) proc_p2 = 0
  if(proc_p2 > NPROCTOT-1) proc_p2 = NPROCTOT-1

! set interval to maximum if user input is not correct
  if(proc_p1 <= 0) proc_p1 = 0
  if(proc_p2 < 0) proc_p2 = NPROCTOT - 1

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *

! open file with global slice number addressing
  write(*,*) 'reading slice addressing'
  write(*,*)
  allocate(ichunk_slice(0:NPROCTOT-1))
  open(unit=IIN,file='OUTPUT_FILES/addressing.txt',status='old')
  do iproc = 0,NPROCTOT-1
    read(IIN,*) iproc_read,ichunk,idummy1,idummy2
    if(iproc_read /= iproc) stop 'incorrect slice number read'
    ichunk_slice(iproc) = ichunk
  enddo
  close(IIN)

! set total number of points and elements to zero
  ntotpoin = 0
  ntotspec = 0

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc = proc_p1,proc_p2

  print *,'Reading slice ',iproc,' in region ',iregion_code

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old')
  read(10,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'
  ntotpoin = ntotpoin + npoin
  close(10)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old')
  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'
  ntotspec = ntotspec + nspec
  close(10)

  enddo
  enddo

  print *
  print *,'There is a total of ',ntotspec,' elements in all the slices'
  print *,'There is a total of ',ntotpoin,' points in all the slices'
  print *

! ************* compute min and max of skewness and ratios ******************

! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL

! set global element and point offsets to zero
  iglobelemoffset = 0

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc,' in region ',iregion_code

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old')

  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
      if(numelem /= ispec) stop 'incorrect element number'

! compute minimum and maximum of quality numbers
      equiangle_skewness_min = dmin1(equiangle_skewness_min,equiangle_skewness)
      edge_aspect_ratio_min = dmin1(edge_aspect_ratio_min,edge_aspect_ratio)
      diagonal_aspect_ratio_min = dmin1(diagonal_aspect_ratio_min,diagonal_aspect_ratio)

      equiangle_skewness_max = dmax1(equiangle_skewness_max,equiangle_skewness)
      edge_aspect_ratio_max = dmax1(edge_aspect_ratio_max,edge_aspect_ratio)
      diagonal_aspect_ratio_max = dmax1(diagonal_aspect_ratio_max,diagonal_aspect_ratio)

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo
  enddo

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions'
  print *
  print *,'equiangle skewness: 0. perfect  1. bad'
  print *,'skewness max deviation angle: 0. perfect  90. bad'
  print *,'skewness min mesh angle: 90. perfect  0. bad'
  print *,'edge aspect ratio: 1. perfect  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect  above 1. gives stretching factor'
  print *,'------------'

  print *
  print *,'equiangle skewness max = ',equiangle_skewness_max
  print *,'equiangle skewness min = ',equiangle_skewness_min
  print *
  print *,'skewness max deviation angle = ',90.*equiangle_skewness_max
  print *,'skewness min mesh angle = ',90.*(1. - equiangle_skewness_max)
  print *
  print *,'edge aspect ratio max = ',edge_aspect_ratio_max
  print *,'edge aspect ratio min = ',edge_aspect_ratio_min
  print *
  print *,'diagonal aspect ratio max = ',diagonal_aspect_ratio_max
  print *,'diagonal aspect ratio min = ',diagonal_aspect_ratio_min
  print *

! create statistics about mesh quality

  print *
  print *,'creating histogram and statistics of mesh quality - reading mesh data files'
  print *

! erase histogram of skewness
  classes_skewness(:) = 0

! erase number of elements belonging to skewness range for AVS_DX
  ntotspecAVS_DX = 0

! set global element and point offsets to zero
  iglobelemoffset = 0

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old')

  read(10,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
      if(numelem /= ispec) stop 'incorrect element number'

! store skewness in histogram
    iclass = int(equiangle_skewness * dble(NCLASS))
    if(iclass < 0) iclass = 0
    if(iclass > NCLASS-1) iclass = NCLASS-1
    classes_skewness(iclass) = classes_skewness(iclass) + 1

! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. &
       equiangle_skewness <= skewness_AVS_DX_max) ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo
  enddo

! create histogram of skewness and save in Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file='OUTPUT_FILES/mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_skewness(iclass))/dble(ntotspec)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewness(iclass),' ',sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

! create script for Gnuplot histogram file
  open(unit=14,file='OUTPUT_FILES/plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term x11'
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASS)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  print *
  print *,'total number of elements = ',ntotspec
  print *,'total percentage = ',total_percent,' %'
  print *

! display warning if maximum skewness is too high
  if(equiangle_skewness_max >= 0.75d0) then
    print *
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, mesh is bad (max skewness >= 0.75)'
    print *,'*********************************************'
    print *,'*********************************************'
    print *
  endif


! ************* create AVS or DX file with elements in a certain range of skewness

  print *
  print *,'creating AVS or DX file with subset of elements in skewness range'
  print *

! ************* count number of points without multiples  ******************

! set global element and point offsets to zero
  iglobpointoffset = 0
  iglobelemoffset = 0

! allocate flag to remove multiples
  allocate(mask_ibool(ntotpoin))
  mask_ibool(:) = .false.

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old')
  open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old')
  open(unit=14,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old')

  read(10,*) nspec
  read(12,*) npoin
  read(14,*) nspec

! read local elements in this slice
  do ispec=1,nspec
    read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
    if(numelem /= ispec) stop 'incorrect element number'

    read(14,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
    if(numelem /= ispec) stop 'incorrect element number'

! check if element belongs to requested skewness range
! and flag all the points to remove multiples
      iglob1 = iglob1 + iglobpointoffset
      iglob2 = iglob2 + iglobpointoffset
      iglob3 = iglob3 + iglobpointoffset
      iglob4 = iglob4 + iglobpointoffset
      iglob5 = iglob5 + iglobpointoffset
      iglob6 = iglob6 + iglobpointoffset
      iglob7 = iglob7 + iglobpointoffset
      iglob8 = iglob8 + iglobpointoffset
      if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) then
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob5) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      endif

  enddo

  iglobelemoffset = iglobelemoffset + nspec
  iglobpointoffset = iglobpointoffset + npoin

  close(10)
  close(12)
  close(14)

  enddo
  enddo

  if(ntotspecAVS_DX == 0) stop 'no elements in skewness range, no file created'

! count number of independent points
  ntotpoinAVS_DX = count(mask_ibool(:))

! write AVS or DX header with element data
  if(USE_OPENDX) then
    open(unit=11,file='OUTPUT_FILES/DX_meshquality.dx',status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ntotpoinAVS_DX,' data follows'
  else
    open(unit=11,file='OUTPUT_FILES/AVS_meshquality.inp',status='unknown')
    write(11,*) ntotpoinAVS_DX,' ',ntotspecAVS_DX,' 0 1 0'
  endif

! ************* generate points ******************

! set global point offset to zero
  iglobpointoffset = 0

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old')
  read(10,*) npoin

! read local points in this slice and output global AVS or DX points
  do ipoin=1,npoin
      read(10,*) numpoin,xval,yval,zval
      if(numpoin /= ipoin) stop 'incorrect point number'
! write to AVS or DX global file with correct offset if point has been selected
      if(mask_ibool(numpoin + iglobpointoffset)) then
        if(USE_OPENDX) then
          write(11,"(f8.5,1x,f8.5,1x,f8.5)") xval,yval,zval
        else
          write(11,"(i6,1x,f8.5,1x,f8.5,1x,f8.5)") numpoin + iglobpointoffset,xval,yval,zval
        endif
      endif

  enddo

  iglobpointoffset = iglobpointoffset + npoin

  close(10)

  enddo
  enddo

! ************* generate elements ******************

! set global element and point offsets to zero
  iglobpointoffset = 0
  iglobelemoffset = 0

  if(USE_OPENDX) write(11,*) 'object 2 class array type int rank 1 shape 8 items ',ntotspecAVS_DX,' data follows'

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old')
  open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old')
  open(unit=14,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old')

  read(10,*) nspec
  read(12,*) npoin
  read(14,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
    read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
    if(numelem /= ispec) stop 'incorrect element number'

    read(14,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
    if(numelem /= ispec) stop 'incorrect element number'

! write to AVS or DX global file with correct offset for hexahedra (3-D)
! check if element belongs to requested skewness range
      iglob1 = iglob1 + iglobpointoffset
      iglob2 = iglob2 + iglobpointoffset
      iglob3 = iglob3 + iglobpointoffset
      iglob4 = iglob4 + iglobpointoffset
      iglob5 = iglob5 + iglobpointoffset
      iglob6 = iglob6 + iglobpointoffset
      iglob7 = iglob7 + iglobpointoffset
      iglob8 = iglob8 + iglobpointoffset
      if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) then
! in the case of OpenDX, node numbers start at zero
! in the case of AVS, node numbers start at one
        if(USE_OPENDX) then
! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
          write(11,200) iglob4-1,iglob1-1,iglob8-1,iglob5-1,iglob3-1,iglob2-1,iglob7-1,iglob6-1
        else
          write(11,201) numelem + iglobelemoffset, &
              iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
        endif
      endif
  enddo

  iglobelemoffset = iglobelemoffset + nspec
  iglobpointoffset = iglobpointoffset + npoin

  close(10)
  close(12)
  close(14)

  enddo
  enddo

 200 format(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)
 201 format(i6,' 1 hex ',i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)

! ************* generate element data values ******************

! output AVS or DX header for data
  if(USE_OPENDX) then
! label for hexahedra in OpenDX is "cubes"
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'
  else
    write(11,*) '1 1'
    write(11,*) 'Zcoord, meters'
  endif

! set global element and point offsets to zero
  iglobelemoffset = 0

  do iregion_code = 1,MAX_NUM_REGIONS

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old')

  read(10,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
      if(numelem /= ispec) stop 'incorrect element number'

! write skewness data to AVS or DX global file with correct offset
! scale skewness to [0:255] for AVS or DX color palette
! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) then
      if(USE_OPENDX) then
        write(11,*) 255.*sngl(equiangle_skewness)
      else
        write(11,*) numelem + iglobelemoffset,' ',255.*sngl(equiangle_skewness)
      endif
    endif

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo
  enddo

! define OpenDX field
  if(USE_OPENDX) then
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'
  endif

! close AVS or DX file
  close(11)

  print *
  print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX skewness range ',skewness_AVS_DX_min,skewness_AVS_DX_max
  print *

  end program mesh_quality_AVS_DX

