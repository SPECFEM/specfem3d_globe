
  program convolve_source_time_function

!
! convolve seismograms computed for a Heaviside with given source time function
!

  implicit none

  double precision, parameter :: DECAY_RATE=2.628d0

  double precision, parameter :: PI=3.141592653589793d0

  integer i,j,N_j
  integer number_remove
  double precision, dimension(:), allocatable :: time,sem
  double precision, dimension(:), allocatable :: sem_fil
  double precision alpha,dt,tau_j,source
  double precision t1,t2,displ1,displ2,gamma,height,exponent

  integer nlines
  double precision hdur

! read file with number of lines in input
  open(unit=33,file='input_convolve_code.txt',status='old')
  read(33,*) nlines
  read(33,*) hdur
  close(33)

! for Gaussian use 1.66667*hdur to get roughly a triangle with half-duration hdur
  hdur = hdur*1.6667

! allocate arrays
  allocate(time(nlines),sem(nlines),sem_fil(nlines))

  do i=1,nlines
    read(5,*) time(i),sem(i)
  enddo

  alpha=DECAY_RATE/hdur
  dt=time(2)-time(1)
  N_j=int(hdur/dt)
  do i=1,nlines
    sem_fil(i)=0.
    do j=-N_j,N_j
      tau_j=dble(j)*dt

   exponent = alpha*alpha*tau_j*tau_j

   if(exponent < 100.) then

! convolve with a Gaussian
!       source = alpha*exp(-exponent)/sqrt(PI)
!using a normalized Ricker
        source = 1. - 2.*exponent*exp(-exponent)

      if(i > j .and. i-j <= nlines) sem_fil(i) = sem_fil(i)+sem(i-j)*source*dt

    endif

    enddo
  enddo

! compute number of samples to remove from end of seismograms
  number_remove = int(hdur / dt) + 1
  do i=1,nlines - number_remove
    write(*,*) sngl(time(i)),' ',sngl(sem_fil(i))
  enddo

  end program convolve_source_time_function

