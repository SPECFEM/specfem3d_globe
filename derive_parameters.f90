!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine derive_parameters(myrank,TMIN,RECORD_LENGTH,NPROC_INPUT, &
     NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSTEP,DT, ONE_CRUST, &
     MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
     NER_CRUST,NER_220_MOHO,NER_400_220, &
     NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
     NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
     NER_DOUBLING_OUTER_CORE, &
     NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB, &
     RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC)

  implicit none

  include 'constants.h'

  integer :: myrank
  double precision :: TMIN,RECORD_LENGTH
  integer ::  NPROC_INPUT,NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSTEP
  logical :: ONE_CRUST
  integer :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
       NER_CRUST,NER_220_MOHO,NER_400_220, &
       NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
       NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
       NER_DOUBLING_OUTER_CORE, &
       NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB
  double precision :: DT
  double precision :: RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC

  integer :: nproc_chunk, nproc1
  character(len=20) :: str

  ! derive xi,eta,radial,dt resolution from input parameters
  if (TMIN < 4.0) then
     call exit_MPI(myrank,'Shortest period < 4.0 s is not supported yet')

  else if (TMIN < 5.0) then
     ! 4.0 -- 5.0 secs.
     NEX_XI                   = 1152
     NEX_ETA                  = 1152
     DT                       = 0.0555555555
     ONE_CRUST                = .false.
     MIN_ATTENUATION_PERIOD   = 4
     MAX_ATTENUATION_PERIOD   = 200
     NER_CRUST                = 4
     NER_220_MOHO             = 8
     NER_400_220              = 7
     NER_600_400              = 7
     NER_670_600              = 3
     NER_771_670              = 4
     NER_TOPDDOUBLEPRIME_771  = 62
     NER_CMB_TOPDDOUBLEPRIME  = 5
     RATIO_TOP_DBL_OC         = 0.48
     RATIO_BOTTOM_DBL_OC      = 0.43
     NER_TOPDBL_CMB           = 38
     NER_ICB_BOTTOMDBL        = 30
     NER_TOP_CENTRAL_CUBE_ICB = 6

  else if (TMIN < 6.75) then
     !  5.0 - 6.75 secs.
     NEX_XI                   = 864
     NEX_ETA                  = 864
     DT                       = 0.072
     MIN_ATTENUATION_PERIOD   = 5
     MAX_ATTENUATION_PERIOD   = 250
     ONE_CRUST                = .false.
     NER_CRUST                = 3
     NER_220_MOHO             = 7
     NER_400_220              = 7
     NER_600_400              = 7
     NER_670_600              = 3
     NER_771_670              = 3
     NER_TOPDDOUBLEPRIME_771  = 46
     NER_CMB_TOPDDOUBLEPRIME  = 4
     RATIO_TOP_DBL_OC         = 0.48
     RATIO_BOTTOM_DBL_OC      = 0.43
     NER_TOPDBL_CMB           = 38
     NER_ICB_BOTTOMDBL        = 30
     NER_TOP_CENTRAL_CUBE_ICB = 6

  else if (TMIN < 8.5) then
     ! 6.75 - 8.5 s
     NEX_XI                   = 640
     NEX_ETA                  = 640
     DT                       = 0.125
     ONE_CRUST                = .false.
     MIN_ATTENUATION_PERIOD   = 5
     MAX_ATTENUATION_PERIOD   = 400
     NER_CRUST                = 2
     NER_220_MOHO             = 5
     NER_400_220              = 4
     NER_600_400              = 4
     NER_670_600              = 2
     NER_771_670              = 2
     NER_TOPDDOUBLEPRIME_771  = 32
     NER_CMB_TOPDDOUBLEPRIME  = 3
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.30
     NER_TOPDBL_CMB           = 22
     NER_ICB_BOTTOMDBL        = 12
     NER_TOP_CENTRAL_CUBE_ICB = 4

  else if (TMIN < 18.0) then
     ! 8.5 - 18 s
     NEX_XI                   = 512
     NEX_ETA                  = 512
     DT                       = 0.125
     ONE_CRUST                = .false.
     MIN_ATTENUATION_PERIOD   = 8
     MAX_ATTENUATION_PERIOD   = 1000
     NER_CRUST                = 2
     NER_220_MOHO             = 5
     NER_400_220              = 4
     NER_600_400              = 4
     NER_670_600              = 2
     NER_771_670              = 2
     NER_TOPDDOUBLEPRIME_771  = 32
     NER_CMB_TOPDDOUBLEPRIME  = 3
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.30
     NER_TOPDBL_CMB           = 22
     NER_ICB_BOTTOMDBL        = 12
     NER_TOP_CENTRAL_CUBE_ICB = 4

  else if (TMIN < 36.0) then
     ! 18 - 36 s
     NEX_XI                   = 240
     NEX_ETA                  = 240
     DT                       = 0.26
     ONE_CRUST                = .true.
     NER_CRUST                = 1
     MIN_ATTENUATION_PERIOD   = 20
     MAX_ATTENUATION_PERIOD   = 1000
     NER_220_MOHO             = 3
     NER_400_220              = 2
     NER_600_400              = 2
     NER_670_600              = 1
     NER_771_670              = 1
     NER_TOPDDOUBLEPRIME_771  = 16
     NER_CMB_TOPDDOUBLEPRIME  = 1
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.27
     NER_TOPDBL_CMB           = 12
     NER_ICB_BOTTOMDBL        = 6
     NER_TOP_CENTRAL_CUBE_ICB = 3

  else if (TMIN < 40.0) then
     ! 36.0 - 40 s
     NEX_XI                   = 160
     NEX_ETA                  = 160
     DT                       = 0.28
     MIN_ATTENUATION_PERIOD   = 20
     MAX_ATTENUATION_PERIOD   = 1000
     ONE_CRUST                = .true.
     NER_CRUST                = 1
     NER_220_MOHO             = 3
     NER_400_220              = 2
     NER_600_400              = 2
     NER_670_600              = 1
     NER_771_670              = 1
     NER_TOPDDOUBLEPRIME_771  = 16
     NER_CMB_TOPDDOUBLEPRIME  = 1
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.27
     NER_TOPDBL_CMB           = 12
     NER_ICB_BOTTOMDBL        = 6
     NER_TOP_CENTRAL_CUBE_ICB = 3

  else if (TMIN < 72.0 ) then
     ! 40.0 - 72.0 s
     NEX_XI                   = 160
     NEX_ETA                  = 160
     DT                       = 0.28
     MIN_ATTENUATION_PERIOD   = 20
     MAX_ATTENUATION_PERIOD   = 1000
     ONE_CRUST                = .true.
     NER_CRUST                = 1
     NER_220_MOHO             = 3
     NER_400_220              = 2
     NER_600_400              = 2
     NER_670_600              = 1
     NER_771_670              = 1
     NER_TOPDDOUBLEPRIME_771  = 16
     NER_CMB_TOPDDOUBLEPRIME  = 1
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.27
     NER_TOPDBL_CMB           = 12
     NER_ICB_BOTTOMDBL        = 6
     NER_TOP_CENTRAL_CUBE_ICB = 3

  else
     ! >= 72.0 s
     NEX_XI                   = 160
     NEX_ETA                  = 160
     DT                       = 0.28
     MIN_ATTENUATION_PERIOD   = 20
     MAX_ATTENUATION_PERIOD   = 1000
     ONE_CRUST                = .true.
     NER_CRUST                = 1
     NER_220_MOHO             = 3
     NER_400_220              = 2
     NER_600_400              = 2
     NER_670_600              = 1
     NER_771_670              = 1
     NER_TOPDDOUBLEPRIME_771  = 16
     NER_CMB_TOPDDOUBLEPRIME  = 1
     RATIO_TOP_DBL_OC         = 0.40
     RATIO_BOTTOM_DBL_OC      = 0.27
     NER_TOPDBL_CMB           = 12
     NER_ICB_BOTTOMDBL        = 6
     NER_TOP_CENTRAL_CUBE_ICB = 3

  endif

  ! check if the number of procs is a multiply of 6.
  if (mod(NPROC_INPUT,6) /=0) &
       call exit_MPI(myrank,'NPROC_INPUT need to be multiple of 6')
  nproc_chunk = NPROC_INPUT/6

  ! calculate number of processors in xi and eta direction.
  nproc1 = nint(sqrt(dble(nproc_chunk)))
  if (nproc_chunk /= nproc1 * nproc1) call exit_MPI(myrank,'NPROC_INPUT/6 should be a square of integer')
  if (mod(NEX_XI/16,nproc1) /= 0) then
     write(str,*) NEX_XI/16
     call exit_MPI(myrank,'Please choose '// &
          'number of input processors nproc_input = 6 * n^2, where n is '// &
          'a factor of '// trim(str))
  endif
  NPROC_XI = nproc1
  NPROC_ETA = NPROC_XI

  ! calculate time_steps
  NSTEP = nint(RECORD_LENGTH/DT)
  if (NSTEP <= 0) &
       call exit_MPI(myrank, 'RECORD_LENGTH should > 0')

  ! check attenuation periods

  if(MIN_ATTENUATION_PERIOD /= 20 .and. MIN_ATTENUATION_PERIOD /= 8 .and. &
       MIN_ATTENUATION_PERIOD /= 5 .and. MIN_ATTENUATION_PERIOD /= 4)  &
       stop 'valid minimum attenuation periods are 20, 8, 5 and 4'

  if(MAX_ATTENUATION_PERIOD /= 1000 .and. MAX_ATTENUATION_PERIOD /= 400) &
       stop 'valid maximum attenuation periods are 1000 and 400'

  ! scale radial mesh parameters according to definitions used in mesher
  ! in order to implement mesh doubling
  NER_220_MOHO = NER_220_MOHO * 2
  NER_400_220 = NER_400_220 * 2
  NER_600_400 = NER_600_400 * 2
  NER_670_600 = NER_670_600 * 2

  NER_771_670 = NER_771_670 * 4
  NER_TOPDDOUBLEPRIME_771 = NER_TOPDDOUBLEPRIME_771 * 4
  NER_CMB_TOPDDOUBLEPRIME = NER_CMB_TOPDDOUBLEPRIME * 4
  NER_TOPDBL_CMB = NER_TOPDBL_CMB * 4
  NER_ICB_BOTTOMDBL = NER_ICB_BOTTOMDBL * 4
  NER_TOP_CENTRAL_CUBE_ICB = NER_TOP_CENTRAL_CUBE_ICB * 4

  NER_ICB_CMB = NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL + NER_TOPDBL_CMB
  NER_DOUBLING_OUTER_CORE = NER_TOP_CENTRAL_CUBE_ICB + NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL

  end subroutine derive_parameters

