program cut_velocity

! this program cuts certain portion of the seismograms and convert them into 
! the adjoints sources for generating banana-dougnut kernels.

  implicit none

  integer :: i, nfile, is, ie, nstep, j, itime ,ifile,ios, i1, i2, nstep_old
  character(len=100) :: arg(100), file(100)
  integer,parameter :: NMAX = 30000
  real*8, parameter :: EPS = 1.0d-17
  real*8, parameter :: PI = 3.1415926d0
  real*8 :: ts, te, data(5,NMAX), out(NMAX), adj(NMAX), tw(NMAX), norm
  real*8 :: dt, t0, t0_old, dt_old, costh, sinth, th, baz
  logical :: lrot

  i = 1
  lrot = .false.

  do while (1 == 1) 
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') stop 'cut_velocity t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    if (trim(arg(i)) == '') exit
    if (i == 1) then 
      read(arg(i),*,iostat=ios) ts
      if (ios /= 0) stop 'Error reading ts'
    else if (i == 2) then
      read(arg(i),*,iostat=ios) te
      if (ios /= 0) stop 'Error reading te'
    else if (i == 3) then
      read(arg(i),*) ifile
      if (ios /= 0) stop 'Error reading ifile'
    else if (i == 4 .or. i == 5 .or. i == 6) then
      file(i-3) = trim(arg(i))
    else if (i == 7) then
      read(arg(i),*,iostat=ios) baz
      if (ios /= 0) stop 'Error reading baz'
      lrot = .true.
    else if (i > 7) then
      stop 'Error: cut_velocity t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    endif
    i = i + 1
  enddo

  i = i - 1
  if (lrot) then
    if (i /= 7) stop 'cut_velocity t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
    if (ifile /= 4 .and. ifile /= 5) stop 'ifile = 4 or 5 when baz is present'
    th = (baz - 180.0) / 180.0 * PI
    costh = cos(th)
    sinth = sin(th)
  else
    if (ifile > 3 .or. ifile < 0) stop 'Error ifile should be between 0 - 3 when baz is not present'
    if (i /= 6) stop 'cut_velocity t1 t2 ifile[0-5] E/N/Z-ascii-files [baz]'
  endif
  
  print *, 'ifile = ', ifile, '  lrot = ', lrot
  print *, ' '

 
  do i = 1, 3
    print *, 'reading asc file '//trim(file(i))//' ...'
    call dread_ascfile_f(file(i),t0,dt,nstep,data(i,:))
    if (nstep > NMAX) stop 'Change the data array range limit'
    if (i == 1) then
      t0_old = t0; dt_old = dt; nstep_old = nstep
    else
      if (i > 1 .and. abs(t0_old - t0) > EPS .and. abs(dt_old - dt) > EPS .and. nstep_old /= nstep) &
                 stop 'Error different t0, dt, nstep'
    endif
  enddo
 
  if (lrot) then
    data(4,:) = costh * data(1,:) - sinth * data(2,:)
    data(5,:) = sinth * data(1,:) + costh * data(2,:)
    call dwrite_ascfile_f('t.txt',t0,dt,nstep,data(4,:))
    call dwrite_ascfile_f('r.txt',t0,dt,nstep,data(5,:))
    i1 = 3; i2 = 5
  else
    i1 = 1; i2 = 3
  endif
    
  
  do i = i1, i2
    is = (ts - t0) / dt + 1
    ie = (te - t0) / dt + 1
    if (is < 1 .or. ie <= is .or. ie > nstep) then
      print *, 'Error in ts, te'; stop
    endif
    tw(1:nstep) = 0.
    do j = is, ie
      tw(j) = 1 - (2 * (dble(j) - is)/(ie - is) - 1) ** 2
    enddo
    do itime = 2, nstep-1
       out(itime) =  (data(i,itime+1) - data(i,itime-1)) / (2 * dt)
    enddo
    out(1) = (data(i,2) - data(i,1)) / dt
    out(nstep) = (data(i,nstep) - data(i,nstep-1)) /dt
    norm = dt * sum( tw(1:nstep) * out(1:nstep) * out(1:nstep))
    print *, 'i = ', i, 'norm = ', norm
    if (ifile /= 0 .and. ifile /= i) norm = 0.0
    
    if (abs(norm) > EPS) then
      adj(1:nstep) = - out(1:nstep) * tw(1:nstep) / norm
    else
      print *, 'norm < EPS for file '//trim(file(i))
      adj(:) = 0.
    endif
    data(i,:) = adj(:)
  enddo
  
  if (lrot) then
    call dwrite_ascfile_f('t-cut.txt',t0,dt,nstep,data(4,:))
    call dwrite_ascfile_f('r-cut.txt',t0,dt,nstep,data(5,:))
    data(1,:) = costh * data(4,:) + sinth * data(5,:)
    data(2,:) = -sinth * data(4,:) + costh * data(5,:)
  endif

  do i = 1, 3
    print *, 'write to asc file '//trim(file(i))//'.ad'
    call dwrite_ascfile_f(trim(file(i))//'.ad',t0,dt,nstep,data(i,:))
  enddo

end program cut_velocity
