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


subroutine write_new_model_perturbations_iso()

! file output for new model perturbations

  use tomography_model_iso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: total_model
  character(len=150) :: m_file

  ! vp relative perturbations
  ! logarithmic perturbation: log( v_new) - log( v_old) = log( v_new / v_old )
  total_model = 0.0_CUSTOM_REAL
  where( model_vp /= 0.0 ) total_model = log( model_vp_new / model_vp)
  ! or
  ! linear approximation: (v_new - v_old) / v_old
  !where( model_vp /= 0.0 ) total_model = ( model_vp_new - model_vp) / model_vp

  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvpvp.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vp)
  call min_all_cr(minval(total_model),min_vp)

  ! vs relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vs /= 0.0 ) total_model = log( model_vs_new / model_vs)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvsvs.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vs)
  call min_all_cr(minval(total_model),min_vs)

  ! rho relative model perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_rho /= 0.0 ) total_model = log( model_rho_new / model_rho)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_drhorho.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_rho)
  call min_all_cr(minval(total_model),min_rho)

  if (myrank == 0) then
    print*,'relative update:'
    print*,'  dvp/vp min/max: ',min_vp,max_vp
    print*,'  dvs/vs min/max: ',min_vs,max_vs
    print*,'  drho/rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine write_new_model_perturbations_iso

!
!-------------------------------------------------------------------------------------------------
!

subroutine write_new_model_perturbations_tiso()

! file output for TI new model perturbations

  use tomography_model_tiso
  implicit none
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta,min_rho,max_rho
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: total_model

  character(len=150) :: m_file

  ! vpv relative perturbations
  ! logarithmic perturbation: log( v_new) - log( v_old) = log( v_new / v_old )
  total_model = 0.0_CUSTOM_REAL
  where( model_vpv /= 0.0 ) total_model = log( model_vpv_new / model_vpv)
  ! or
  ! linear approximation: (v_new - v_old) / v_old
  !where( model_vpv /= 0.0 ) total_model = ( model_vpv_new - model_vpv) / model_vpv

  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvpvvpv.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vpv)
  call min_all_cr(minval(total_model),min_vpv)

  ! vph relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vph /= 0.0 ) total_model = log( model_vph_new / model_vph)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvphvph.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vph)
  call min_all_cr(minval(total_model),min_vph)

  ! vsv relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsv /= 0.0 ) total_model = log( model_vsv_new / model_vsv)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvsvvsv.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vsv)
  call min_all_cr(minval(total_model),min_vsv)

  ! vsh relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_vsh /= 0.0 ) total_model = log( model_vsh_new / model_vsh)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_dvshvsh.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_vsh)
  call min_all_cr(minval(total_model),min_vsh)

  ! eta relative perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_eta /= 0.0 ) total_model = log( model_eta_new / model_eta)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_detaeta.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_eta)
  call min_all_cr(minval(total_model),min_eta)

  ! rho relative model perturbations
  total_model = 0.0_CUSTOM_REAL
  where( model_rho /= 0.0 ) total_model = log( model_rho_new / model_rho)
  write(m_file,'(a,i6.6,a)') 'OUTPUT_MODEL/proc',myrank,'_reg1_drhorho.bin'
  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) total_model
  close(IOUT)
  call max_all_cr(maxval(total_model),max_rho)
  call min_all_cr(minval(total_model),min_rho)

  if (myrank == 0) then
    print*,'relative update:'
    print*,'  dvpv/vpv min/max: ',min_vpv,max_vpv
    print*,'  dvph/vph min/max: ',min_vph,max_vph
    print*,'  dvsv/vsv min/max: ',min_vsv,max_vsv
    print*,'  dvsh/vsh min/max: ',min_vsh,max_vsh
    print*,'  deta/eta min/max: ',min_eta,max_eta
    print*,'  drho/rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine write_new_model_perturbations_tiso

