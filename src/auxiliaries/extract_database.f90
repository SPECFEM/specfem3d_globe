!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


  program extract_databases

! right now just extract the mantle databases
!
! Usage: extract_databases proc000000_reg1_solver_data.bin ireg n  proc000000_reg1_kappa.bin
! for example: n = 10 (rho), 11 (kappav), 12 (muv)
!
!  Qinya Liu, Caltech, May 2007

  use constants, only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,MAX_STRING_LEN

  use shared_parameters, only: &
    LOCAL_PATH

  use constants_solver, only: &
    NGLOB_CRUST_MANTLE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_INNER_CORE,NSPEC_OUTER_CORE

  implicit none

  character(len=MAX_STRING_LEN) :: infile, s_num, outfile, s_ireg
  integer :: num, i, nspec, nglob, ireg

  ! uses nspec_crust_mantle for allocation as this is the maximum possible size
  integer,dimension(:,:,:,:),allocatable :: idummy_sem ! ibool
  integer,dimension(:),allocatable :: idummy_arr  !idoubling
  logical,dimension(:),allocatable :: ldummy_arr  ! ispec_is_tiso

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: junk
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: junk1
  !real(kind=CUSTOM_REAL) :: junk2

  integer :: ier
  integer :: idummy

  ! gets arguments
  call get_command_argument(1,infile)
  call get_command_argument(2,s_ireg)
  call get_command_argument(3,s_num)
  call get_command_argument(4,outfile)

  if (len_trim(infile) == 0 .or. len_trim(s_ireg) == 0 .or. len_trim(s_num) == 0 .or. len_trim(outfile) == 0) then
     print *, 'Usage: extract_databases infile ireg num outfile'
     print *, '  ireg = 1, 2, 3'
     print *, '  num = 10 for rho,  11 for kappav, 12 for muv '
     stop ' Reenter command line options'
  endif

  read(s_ireg,*) ireg
  read(s_num,*) num

  ! user output
  print *,"Extracting database values:"
  print *,"  database file : ",trim(infile)
  print *,"  region code   : ",ireg
  print *,"  array number  : ",num
  print *

  ! reads mesh parameters
  ! index of last occurrence of '/' to get directory from name (DATABASES_MPI/proc****_reg1_solver_data.bin)
  LOCAL_PATH = ''
  i = index(infile,'/',.true.)
  if (i > 1) then
    LOCAL_PATH(1:i-1) = infile(1:i-1)
  else
    LOCAL_PATH = 'DATABASES_MPI'        ! mesh_parameters.bin file in DATABASES_MPI/
  endif
  call read_mesh_parameters()

  ! user output
  print *,'mesh parameters (from input directory):'
  print *,'  input directory    = ',trim(LOCAL_PATH)
  print *,'  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
  print *,'  NSPEC_OUTER_CORE   = ',NSPEC_OUTER_CORE
  print *,'  NSPEC_INNER_CORE   = ',NSPEC_INNER_CORE
  print *

  if (ireg == 1) then
    nspec = NSPEC_CRUST_MANTLE
    nglob = NGLOB_CRUST_MANTLE
  else if (ireg == 2) then
    nspec = NSPEC_OUTER_CORE
    nglob = NGLOB_OUTER_CORE
  else if (ireg == 3) then
    nspec = NSPEC_INNER_CORE
    nglob = NGLOB_INNER_CORE
  else
    stop 'Error: ireg has to be 1, 2 or 3'
  endif

  if (num < 10 .or. num > 12) then
    stop 'Error: num has to be 10, 11 or 12'
  endif

  if (num == 10) then
    print *,"  extracting for: ","rhostore"
  else if (num == 11) then
    print *,"  extracting for: ","kappavstore"
  else if (num == 12) then
    print *,"  extracting for: ","muvstore"
  else
    print *,"  extracting for: ","unknown store"
  endif
  print *
  print *,"  output to file: ",trim(outfile)
  print *

  allocate(idummy_sem(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE)) ! ibool
  allocate(idummy_arr(NSPEC_CRUST_MANTLE))  !idoubling
  allocate(ldummy_arr(NSPEC_CRUST_MANTLE))  ! ispec_is_tiso

  allocate(junk(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           junk1(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating dummy arrays'

  ! opens solver_data file
  open(11,file=trim(infile),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening input file'

  ! new solver_data.bin: ordering see save_arrays_solver.f90
  ! 0.a nspec
  ! 0.b nglob
  ! 1. xstore
  ! 2. ystore
  ! 3. zstore
  ! 4. ibool - integer array
  ! 5. idoubling - integer array
  ! 6. ispec_is_tiso - logical array
  ! 7. xixstore
  ! 8. xiystore
  ! 9. xizstore
  ! 10. etaxstore
  ! 11. etaystore
  ! 12. etazstore
  ! 13. gammaxstore
  ! 14. gammaystore
  ! 15. gammazstore
  ! 16. rhostore
  ! 17. kappavstore

  ! skipps nspec
  read(11) idummy
  if (idummy /= nspec) stop 'Error invalid nspec in input file (solver_data.bin)'
  ! skipps nglob
  read(11) idummy
  if (idummy /= nglob) stop 'Error invalid nglob in input file (solver_data.bin)'

  ! user output
  print *,"solver data:"
  print *,"  number of spectral elements = ",nspec
  print *,"  number of global points = ",nglob
  print *

  ! xstore,ystore,zstore
  read(11) junk1(1:nglob)
  print *,"  xstore = ",junk1(1:5),".."

  read(11) junk1(1:nglob)
  print *,"  ystore = ",junk1(1:5),".."

  read(11) junk1(1:nglob)
  print *,"  zstore = ",junk1(1:5),".."

  ! ibool
  read(11) idummy_sem(:,:,:,1:nspec)

  ! idoubling
  read(11) idummy_arr(1:nspec)

  ! ispec_is_tiso
  read(11) ldummy_arr(1:nspec)

  ! skipps arrays, i.e. xixstore,..
  do i = 1, num-1
    read(11) junk(:,:,:,1:nspec)
  enddo

  ! reads desired data array, e.g. rhostore,..
  read(11) junk(:,:,:,1:nspec)
  print *
  print *,"  desired data = ",junk(:,1,1,1),".."

  ! closes file
  close(11)

  ! writes out data
  open(12,file=trim(outfile),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening output file'
  write(12) junk(:,:,:,1:nspec)
  close(12)

  ! user output
  print *
  print *,"done extracting, see file: ",trim(outfile)
  print *

  ! frees arrays
  deallocate(idummy_sem,idummy_arr,ldummy_arr)
  deallocate(junk,junk1)

  end program extract_databases

