
  program convolve_jeroen

!! January 2001, by Jeroen Tromp

!! Dimitri Komatitsch: added option to convolve with a triangle
!!                     also now reads hdur and number of lines from external file

!! suppressed last "hdur" seconds in seismograms
!! to remove artefact from the convolution at the end of the seismogram
!! which was then creating some problems in the "remove_offset" routine

  implicit none

  double precision, parameter :: DECAY_RATE = 2.628d0

  double precision, parameter :: PI = 3.141592653589793d0

! number of samples in the ASCII seismogram files (i.e. it is the value of NSTEP in the SPECFEM3D code)
  integer, parameter :: nlines = 15500

  integer i,j,N_j
  integer number_remove
  real, dimension(:), allocatable :: time,sem
  double precision, dimension(:), allocatable :: sem_fil
  double precision alpha,dt,tau_j,source
  double precision t1,t2,displ1,displ2,gamma,height

  double precision hdur
  logical triangle

! read file with number of lines in input
! open(unit=33,file='input_convolve_code.txt',status='old')
! read(33,*) nlines
! read(33,*) hdur
! read(33,*) triangle
! close(33)
  hdur = 25
  triangle = .false.

! allocate arrays
  allocate(time(nlines),sem(nlines),sem_fil(nlines))

! read the seismogram from the standard input
  do i = 1,nlines
    read(*,*) time(i),sem(i)
  enddo

  alpha = DECAY_RATE/hdur
  dt = dble(time(2))-dble(time(1))
  N_j=int(hdur/dt)
  do i = 1,nlines
    sem_fil(i)=0.0d0
    do j=-N_j,N_j
      tau_j=dble(j)*dt

! convolve with a triangle
    if (triangle) then
       height = 1. / hdur
       if (abs(tau_j) > hdur) then
         source = 0.
       else if (tau_j < 0) then
         t1 = - N_j * dt
         displ1 = 0.
         t2 = 0
         displ2 = height
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1. - gamma) * displ1 + gamma * displ2
       else
         t1 = 0
         displ1 = height
         t2 = + N_j * dt
         displ2 = 0.
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1. - gamma) * displ1 + gamma * displ2
       endif

      else

! convolve with a Gaussian
        source = alpha*dexp(-alpha*alpha*tau_j*tau_j)/dsqrt(PI)

      endif

      if (i > j .and. i-j <= nlines) &
              sem_fil(i)=sem_fil(i)+dble(sem(i-j))*source*dt
    enddo
  enddo

! compute number of samples to remove from end of seismograms
  number_remove = int(hdur / dt) + 1
! write the convolved seismogram to the standard output
  do i = 1,nlines - number_remove
    write(*,*) time(i),sngl(sem_fil(i))
  enddo

  end program convolve_jeroen

