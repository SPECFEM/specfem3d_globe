program test_rotate

  implicit none

  double precision, parameter :: PI = 3.141592653589793d0
  integer,parameter :: CUSTOM_REAL = 4

  ! local parameters
  double precision :: A,C,L,N,F,Gc,Gs
  double precision :: Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es
  double precision :: theta,phi
  double precision :: A_out,C_out,L_out,N_out,F_out

  double precision :: vpv,vph,vsv,vsh,rho,eta_aniso

  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  real(kind=CUSTOM_REAL) :: A_cr,C_cr,L_cr,N_cr,eta_cr,theta_cr,phi_cr
  real(kind=CUSTOM_REAL) :: cr11,cr12,cr13,cr14,cr15,cr16,cr22,cr23,cr24,cr25,cr26, &
                            cr33,cr34,cr35,cr36,cr44,cr45,cr46,cr55,cr56,cr66
  real(kind=CUSTOM_REAL), dimension(21) :: cij_kl
  real(kind=CUSTOM_REAL), dimension(21) :: cij_kl_spherical

  print *,'program: test_rotate'

  ! non-dimensionalized values
  vph = 1.2021723985671997
  vpv = 1.1566073894500732
  vsh = 0.68188071250915527
  vsv = 0.66188687086105347
  eta_aniso = 0.91191035509109497
  rho = 0.60045558214187622

  print *
  print *,'Checks:'
  print *

  ! Love parameterization
  A = rho * vph**2
  C = rho * vpv**2
  L = rho * vsv**2
  N = rho * vsh**2
  F = eta_aniso * (A - 2.d0 * L)

  ! angle in radians
  ! pointing in z-direction, rotation along phi should be invariant for tiso
  theta = 0.d0
  phi = PI / 6.d0

  print *,'test invariance: Love to global'
  print *,'  theta/phi = ',theta*180.0/PI,phi*180.0/PI

  ! radial axis symmetry:
  ! C11 = A = rho * vph**2
  ! C33 = C = rho * vpv**2
  ! C44 = L = rho * vsv**2
  ! C13 = F = eta * (A - 2*L)
  ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
  ! C22 = C11
  ! C23 = C13
  ! C55 = C44
  ! C66 = N = rho * vsh**2 = (C11-C12)/2

  ! local (radial) coordinate system to global SPECFEM reference
  call rotate_tensor_Love_to_global(theta,phi, &
                                    A,C,N,L,F, &
                                    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
  ! rotates back
  call rotate_tensor_global_to_radial(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! check
  print *,'check: Love to global'
  call check_par(A,c11,"c11")
  call check_par(C,c33,"c33")
  call check_par(L,c44,"c44")
  call check_par(N,c66,"c66")
  call check_par(F,c13,"c13")
  call check_par(A-2.d0*N,c12,"c12")
  call check_par(c11,c22,"c22")
  call check_par(c13,c23,"c23")
  call check_par(c44,c55,"c55")
  print *,'  A,C,L,N,F all good'
  print *
  print *,'check: Love global to radial back-rotation'
  call check_par(A,d11,"d11")
  call check_par(C,d33,"d33")
  call check_par(L,d44,"d44")
  call check_par(N,d66,"d66")
  call check_par(F,d13,"d13")
  call check_par(A-2.d0*N,d12,"d12")
  call check_par(c11,d22,"d22")
  call check_par(c13,d23,"d23")
  call check_par(c44,d55,"d55")
  print *,'  back A,C,L,N,F all good'
  print *
  print *,'check: Love global to radial back-rotation invariance'
  call check_par(c11,d11,"d11")
  call check_par(c12,d12,"d12")
  call check_par(c13,d13,"d13")
  call check_par(c14,d14,"d14")
  call check_par(c15,d15,"d15")
  call check_par(c16,d16,"d16")
  call check_par(c22,d22,"d22")
  call check_par(c23,d23,"d23")
  call check_par(c24,d24,"d24")
  call check_par(c25,d25,"d25")
  call check_par(c26,d26,"d26")
  call check_par(c33,d33,"d33")
  call check_par(c34,d34,"d34")
  call check_par(c35,d35,"d35")
  call check_par(c36,d36,"d36")
  call check_par(c44,d44,"d44")
  call check_par(c45,d45,"d45")
  call check_par(c46,d46,"d46")
  call check_par(c55,d55,"d55")
  call check_par(c56,d56,"d56")
  call check_par(c66,d66,"d66")
  print *,'  back all good'
  print *

  ! oblique angle in radians
  ! should lead to different, full cij
  theta = PI / 6.d0
  phi = PI / 6.d0

  print *,'test oblique: Love to global - global to Love'
  print *,'  theta/phi = ',theta*180.0/PI,phi*180.0/PI

  call rotate_tensor_Love_to_global(theta,phi, &
                                    A,C,N,L,F, &
                                    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  call rotate_tensor_global_to_Love(theta,phi, &
                                    A_out,C_out,N_out,L_out,F_out, &
                                    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  print *,'check: Love to global to Love back-rotation'
  call check_par(A,A_out,"A")
  call check_par(C,C_out,"C")
  call check_par(L,L_out,"L")
  call check_par(N,N_out,"N")
  call check_par(F,F_out,"F")
  print *,'  back A,C,L,N,F all good'
  print *

  call rotate_tensor_global_to_radial(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  print *,'check: Love to global to radial back-rotation'
  call check_par(A,d11,"d11")
  call check_par(C,d33,"d33")
  call check_par(L,d44,"d44")
  call check_par(N,d66,"d66")
  call check_par(F,d13,"d13")
  call check_par(A-2.d0*N,d12,"d12")
  call check_par(d11,d22,"d22")
  call check_par(d13,d23,"d23")
  call check_par(d44,d55,"d55")
  print *,'  back A,C,L,N,F all good'
  print *


  ! tiso routine
  ! CUSTOM_REAL parameters
  theta_cr = theta
  phi_cr = phi
  eta_cr = eta_aniso
  A_cr = A
  C_cr = C
  L_cr = L
  N_cr = N
  call rotate_tensor_tiso_to_cij(theta_cr,phi_cr, &
                                 A_cr,C_cr,L_cr,N_cr,eta_cr, &
                                 cr11,cr12,cr13,cr14,cr15,cr16,cr22,cr23,cr24,cr25,cr26, &
                                 cr33,cr34,cr35,cr36,cr44,cr45,cr46,cr55,cr56,cr66)
  ! converts back to double
  print *,'check: rotate_tensor_tiso_to_cij() routine'
  call check_par_cr(c11,cr11,"cr11")
  call check_par_cr(c12,cr12,"cr12")
  call check_par_cr(c13,cr13,"cr13")
  call check_par_cr(c14,cr14,"cr14")
  call check_par_cr(c15,cr15,"cr15")
  call check_par_cr(c16,cr16,"cr16")
  call check_par_cr(c22,cr22,"cr22")
  call check_par_cr(c23,cr23,"cr23")
  call check_par_cr(c24,cr24,"cr24")
  call check_par_cr(c25,cr25,"cr25")
  call check_par_cr(c26,cr26,"cr26")
  call check_par_cr(c33,cr33,"cr33")
  call check_par_cr(c34,cr34,"cr34")
  call check_par_cr(c35,cr35,"cr35")
  call check_par_cr(c36,cr36,"cr36")
  call check_par_cr(c44,cr44,"cr44")
  call check_par_cr(c45,cr45,"cr45")
  call check_par_cr(c46,cr46,"cr46")
  call check_par_cr(c55,cr55,"cr55")
  call check_par_cr(c56,cr56,"cr56")
  call check_par_cr(c66,cr66,"cr66")
  print *,'  all good'
  print *

  ! vector routine for kernels
  cij_kl(1) = c11
  cij_kl(2) = c12
  cij_kl(3) = c13
  cij_kl(4) = c14
  cij_kl(5) = c15
  cij_kl(6) = c16
  cij_kl(7) = c22
  cij_kl(8) = c23
  cij_kl(9) = c24
  cij_kl(10) = c25
  cij_kl(11) = c26
  cij_kl(12) = c33
  cij_kl(13) = c34
  cij_kl(14) = c35
  cij_kl(15) = c36
  cij_kl(16) = c44
  cij_kl(17) = c45
  cij_kl(18) = c46
  cij_kl(19) = c55
  cij_kl(20) = c56
  cij_kl(21) = c66

  call rotate_tensor_global_to_radial_vector(cij_kl,cij_kl_spherical,theta_cr,phi_cr)

  print *,'check: rotate_tensor_global_to_radial_vector() routine'
  call check_par_cr(d11,cij_kl_spherical(1),"cij_kl_spherical(1)")
  call check_par_cr(d12,cij_kl_spherical(2),"cij_kl_spherical(2)")
  call check_par_cr(d13,cij_kl_spherical(3),"cij_kl_spherical(3)")
  call check_par_cr(d14,cij_kl_spherical(4),"cij_kl_spherical(4)")
  call check_par_cr(d15,cij_kl_spherical(5),"cij_kl_spherical(5)")
  call check_par_cr(d16,cij_kl_spherical(6),"cij_kl_spherical(6)")
  call check_par_cr(d22,cij_kl_spherical(7),"cij_kl_spherical(7)")
  call check_par_cr(d23,cij_kl_spherical(8),"cij_kl_spherical(8)")
  call check_par_cr(d24,cij_kl_spherical(9),"cij_kl_spherical(9)")
  call check_par_cr(d25,cij_kl_spherical(10),"cij_kl_spherical(10)")
  call check_par_cr(d26,cij_kl_spherical(11),"cij_kl_spherical(11)")
  call check_par_cr(d33,cij_kl_spherical(12),"cij_kl_spherical(12)")
  call check_par_cr(d34,cij_kl_spherical(13),"cij_kl_spherical(13)")
  call check_par_cr(d35,cij_kl_spherical(14),"cij_kl_spherical(14)")
  call check_par_cr(d36,cij_kl_spherical(15),"cij_kl_spherical(15)")
  call check_par_cr(d44,cij_kl_spherical(16),"cij_kl_spherical(16)")
  call check_par_cr(d45,cij_kl_spherical(17),"cij_kl_spherical(17)")
  call check_par_cr(d46,cij_kl_spherical(18),"cij_kl_spherical(18)")
  call check_par_cr(d55,cij_kl_spherical(19),"cij_kl_spherical(19)")
  call check_par_cr(d56,cij_kl_spherical(20),"cij_kl_spherical(20)")
  call check_par_cr(d66,cij_kl_spherical(21),"cij_kl_spherical(21)")
  print *,'  all good'
  print *


  ! local (azimuthal) coordinate system to global SPECFEM reference
  Gc = 0.d0
  Gs = 0.d0
  call rotate_tensor_azimuthal_to_global(theta,phi, &
                                         A,C,N,L,F, &
                                         Gc,Gs, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                         c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  call rotate_tensor_global_to_Love(theta,phi, &
                                    A_out,C_out,N_out,L_out,F_out, &
                                    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  print *,'check: rotate_tensor_azimuthal_to_global() routine'
  call check_par(A,A_out,"A")
  call check_par(C,C_out,"C")
  call check_par(L,L_out,"L")
  call check_par(N,N_out,"N")
  call check_par(F,F_out,"F")
  print *,'  back A,C,L,N,F all good'
  print *


  ! full azimuthal to global reference
  Jc = 0.d0
  Js = 0.d0
  Kc = 0.d0
  Ks = 0.d0
  Mc = 0.d0
  Ms = 0.d0
  Bc = 0.d0
  Bs = 0.d0
  Hc = 0.d0
  Hs = 0.d0
  Dc = 0.d0
  Ds = 0.d0
  Ec = 0.d0
  Es = 0.d0
  call rotate_tensor_aniso_to_global(theta,phi, &
                                     A,C,N,L,F, &
                                     Gc,Gs, &
                                     Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es, &
                                     c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                     c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  call rotate_tensor_global_to_Love(theta,phi, &
                                    A_out,C_out,N_out,L_out,F_out, &
                                    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  print *,'check: rotate_tensor_aniso_to_global() routine'
  call check_par(A,A_out,"A")
  call check_par(C,C_out,"C")
  call check_par(L,L_out,"L")
  call check_par(N,N_out,"N")
  call check_par(F,F_out,"F")
  print *,'  back A,C,L,N,F all good'
  print *

  ! done
  print *,'test_rotate done successfully'

end program test_rotate

!----------------------------------------------------------------------------------------

  subroutine check_par(par_in,par_test,name)

  implicit none
  double precision :: par_in,par_test
  character(len=*) :: name

  double precision, parameter :: THRES = 1.d-15

  print *,'  ',trim(name),' = ',par_in,'/',par_test

  if (abs(par_in - par_test) > THRES) then
    print *,'Error: parameters differ!'
    print *,'  parameter  name = ',trim(name)
    print *,'            value = ',par_test,' should be equal to ',par_in
    print *,'exiting...'
    stop 1
  endif

  end subroutine check_par


!----------------------------------------------------------------------------------------

  subroutine check_par_cr(par_in,par_test,name)

  implicit none
  integer, parameter :: CUSTOM_REAL = 4

  double precision :: par_in
  real(kind=CUSTOM_REAL) :: par_test     ! custom real
  character(len=*) :: name

  double precision, parameter :: THRES = 1.d-6   ! lower accuracy for custom real

  print *,'  ',trim(name),' = ',par_in,'/',par_test

  if (abs(par_in - par_test) > THRES) then
    print *,'Error: parameters differ!'
    print *,'  parameter  name = ',trim(name)
    print *,'            value = ',par_test,' should be equal to ',par_in
    print *,'exiting...'
    stop 1
  endif

  end subroutine check_par_cr

