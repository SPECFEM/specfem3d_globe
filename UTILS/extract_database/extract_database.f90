!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, November 2007
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
! Usage: extract_databases proc000000_reg1_solver_data_1.bin ireg n  proc000000_reg1_kappa.bin
! for example: n = 10 (rho), 11 (kappav), 12 (muv)

  implicit none

  include 'constants.h'
  include 'OUTPUT_FILES/values_from_mesher.h'

  character(len=150) :: infile, s_num, outfile, s_ireg
  integer :: num, i, nspec, ireg

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: junk
  real(kind=CUSTOM_REAL) :: junk2

  call getarg(1,infile)
  call getarg(2,s_ireg)
  call getarg(3,s_num)
  call getarg(4,outfile)
  if (trim(infile) == '' .or. trim(s_ireg) == '' .or. trim(s_num) == '' &
   .or. trim(outfile) == '') then 
     print *, 'Usage: extract_databases infile ireg num outfile'
     print *, '  ireg = 1, 2, 3'
     print *, '  num = 10 for rho,  11 for kappav, 12 for muv '
     stop ' Reenter command line options'
  endif

  read(s_ireg,*) ireg
  if (ireg == 1) then
    nspec=NSPEC_CRUST_MANTLE
  else if (ireg == 2) then
    nspec=NSPEC_OUTER_CORE
  else if (ireg == 3) then
    nspec=NSPEC_INNER_CORE
  else
    stop 'Error: ireg has to be 1, 2 or 3'
  endif

  read(s_num,*) num

  open(11,file=infile,status='old',form='unformatted')
  do i = 1, num-1
    read(11) junk2
  enddo
  read(11) junk(:,:,:,1:nspec)
  close(11)
  
  open(12,file=outfile,status='unknown',form='unformatted')
  write(12) junk(:,:,:,1:nspec)
  close(12)

end program extract_databases

