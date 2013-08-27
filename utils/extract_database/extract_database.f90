program extract_databases

! right now just extract the mantle databases
! Usage: extract_databases proc000000_reg1_solver_data_1.bin ireg n  proc000000_reg1_kappa.bin
! for example: n = 10 (rho), 11 (kappav), 12 (muv)
!  Qinya Liu, Caltech, May 2007

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  character(len=150) :: infile, s_num, outfile, s_ireg
  integer :: num, i, nspec, ireg

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: junk
  real(kind=CUSTOM_REAL) :: junk2

  call get_command_argument(1,infile)
  call get_command_argument(2,s_ireg)
  call get_command_argument(3,s_num)
  call get_command_argument(4,outfile)
  if (len_trim(infile) == 0 .or. len_trim(s_ireg) == 0 .or. len_trim(s_num) == 0 &
   .or. len_trim(outfile) == 0) then
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

