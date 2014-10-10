!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

subroutine read_model_iso()

! reads in current isotropic model: vp & vs & rho

  use tomography_model_iso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  integer :: ival,ier
  character(len=150) :: m_file

  ! reads in current vp & vs & rho model
  ! vp model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vp.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vp(:,:,:,1:nspec)
  close(IIN)

  ! vs model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vs.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vs(:,:,:,1:nspec)
  close(IIN)

  ! rho model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_rho.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_rho(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(model_vp),min_vp)
  call max_all_cr(maxval(model_vp),max_vp)

  call min_all_cr(minval(model_vs),min_vs)
  call max_all_cr(maxval(model_vs),max_vs)

  call min_all_cr(minval(model_rho),min_rho)
  call max_all_cr(maxval(model_rho),max_rho)

  if (myrank == 0) then
    print*,'initial models:'
    print*,'  vs min/max: ',min_vs,max_vs
    print*,'  vp min/max: ',min_vp,max_vp
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

  ! global addressing
  write(m_file,'(a,i6.6,a)') 'topo/proc',myrank,'_reg1_solver_data.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif

  read(IIN) ival !nspec
  if (ival /= nspec) call exit_mpi(myrank,'Error invalid nspec value in solver_data.bin')
  read(IIN) ival !nglob
  if (ival /= nglob) call exit_mpi(myrank,'Error invalid nspec value in solver_data.bin')

  read(IIN) x(1:nglob)
  read(IIN) y(1:nglob)
  read(IIN) z(1:nglob)
  read(IIN) ibool(:,:,:,1:nspec)
  close(IIN)

end subroutine read_model_iso

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_model_tiso()

! reads in current transverse isotropic model: vpv.. & vsv.. & eta & rho

  use tomography_model_tiso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta,min_rho,max_rho
  integer :: ival,ier
  character(len=150) :: m_file

  ! vpv model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vpv.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vpv(:,:,:,1:nspec)
  close(IIN)

  ! vph model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vph.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vph(:,:,:,1:nspec)
  close(IIN)

  ! vsv model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vsv.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vsv(:,:,:,1:nspec)
  close(IIN)

  ! vsh model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_vsh.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_vsh(:,:,:,1:nspec)
  close(IIN)

  ! eta model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_eta.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_eta(:,:,:,1:nspec)
  close(IIN)

  ! rho model
  write(m_file,'(a,i6.6,a)') 'INPUT_MODEL/proc',myrank,'_reg1_rho.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) model_rho(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(model_vpv),min_vpv)
  call max_all_cr(maxval(model_vpv),max_vpv)

  call min_all_cr(minval(model_vph),min_vph)
  call max_all_cr(maxval(model_vph),max_vph)

  call min_all_cr(minval(model_vsv),min_vsv)
  call max_all_cr(maxval(model_vsv),max_vsv)

  call min_all_cr(minval(model_vsh),min_vsh)
  call max_all_cr(maxval(model_vsh),max_vsh)

  call min_all_cr(minval(model_eta),min_eta)
  call max_all_cr(maxval(model_eta),max_eta)

  call min_all_cr(minval(model_rho),min_rho)
  call max_all_cr(maxval(model_rho),max_rho)

  if (myrank == 0) then
    print*,'initial models:'
    print*,'  vpv min/max: ',min_vpv,max_vpv
    print*,'  vph min/max: ',min_vph,max_vph
    print*,'  vsv min/max: ',min_vsv,max_vsv
    print*,'  vsh min/max: ',min_vsh,max_vsh
    print*,'  eta min/max: ',min_eta,max_eta
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

  ! global addressing
  write(m_file,'(a,i6.6,a)') 'topo/proc',myrank,'_reg1_solver_data.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif

  read(IIN) ival ! nspec
  if (ival /= nspec) call exit_mpi(myrank,'Error invalid nspec value in solver_data.bin')
  read(IIN) ival ! nglob
  if (ival /= nglob) call exit_mpi(myrank,'Error invalid nglob value in solver_data.bin')

  read(IIN) x(1:nglob)
  read(IIN) y(1:nglob)
  read(IIN) z(1:nglob)
  read(IIN) ibool(:,:,:,1:nspec)
  read(IIN) idoubling(1:nspec)
  read(IIN) ispec_is_tiso(1:nspec)
  close(IIN)

end subroutine read_model_tiso

