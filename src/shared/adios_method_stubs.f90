!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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


! placeholders for non-adios compilation

! for xmeshfem3D compilation

  subroutine adios_cleanup()
  end subroutine

  subroutine adios_setup()
  end subroutine

  subroutine crm_save_mesh_files_adios()
  end subroutine

  subroutine get_absorb_adios()
  end subroutine

  subroutine save_arrays_solver_adios()
  end subroutine

  subroutine save_mpi_arrays_adios()
  end subroutine


! for xspecfem3D compilation

  subroutine read_arrays_solver_adios()
  end subroutine

  subroutine read_attenuation_adios()
  end subroutine

  subroutine read_forward_arrays_adios()
  end subroutine

  subroutine read_intermediate_forward_arrays_adios()
  end subroutine

  subroutine read_mesh_databases_coupling_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_cm_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_ic_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_oc_adios()
  end subroutine

  subroutine read_mesh_databases_stacey_adios()
  end subroutine

  subroutine save_forward_arrays_adios()
  end subroutine

  subroutine save_intermediate_forward_arrays_adios()
  end subroutine

