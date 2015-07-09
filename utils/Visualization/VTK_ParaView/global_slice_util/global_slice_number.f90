program global_slice_number

! this program computes the slice numbers the source-receiver lin goes through
! for given source receiver location and simulation parameters
! Qinya Liu, Caltech, May 2007

  use sub_slice_number

  implicit none

  character(len=200) :: ch_elon, ch_elat, ch_slon, ch_slat, ch_nproc, ch_narc, ch_usage, ch_arg7, ch_arg8, &
             ch_eta_width, ch_clon, ch_clat, ch_grot, ch_arg13
  logical :: belt_region
  integer :: nchunks, narc, nproc, np, nt, i, ii, k, j, nproc_xi, nproc_eta, np_short
  integer, dimension(NPTS) :: chunk, slices, slices_short

  real(kind=CUSTOM_REAL) :: elon, elat, slon, slat, xi_width, eta_width, lat0, th0, clon, clat, grot, &
             ths, phs, thr, phr, thn, phn, ths_new, phs_new, thr_new, phr_new, delta, rotation_matrix(3,3), &
             x, y, z, xn, yn, zn, xik, etak, xi1, eta1
  real(kind=CUSTOM_REAL), dimension(NPTS) :: th, ph, thp, php, xi,eta

  !============================

  ! set up the input arguments

  call getarg(1,ch_elon)
  call getarg(2,ch_elat)
  call getarg(3,ch_slon)
  call getarg(4,ch_slat)
  call getarg(5,ch_nproc)
  call getarg(6,ch_narc)

  ch_usage = 'Usage: global_slice_number elon elat slon slat nproc narc [nchunks xi_width eta_width clon clat grot] [lat0]'

  if (trim(ch_elon) == '' .or. trim(ch_elat) == '' .or. trim(ch_slon) == '' .or. &
             trim(ch_slat) == '' .or. trim(ch_nproc) == '' .or. trim(ch_narc) == '')  then
    print *, trim(ch_usage)
    stop
  endif

  read(ch_elon,*) elon
  read(ch_elat,*) elat
  read(ch_slon,*) slon
  read(ch_slat,*) slat
  read(ch_nproc,*) nproc
  read(ch_narc,*) narc

  call getarg(7,ch_arg7)
  if (trim(ch_arg7) == '') then
    belt_region = .false.
    nchunks = 6
    xi_width = 90.
    eta_width = 90.
  else
    call getarg(8,ch_arg8)
    if (trim(ch_arg8) == '') then
      belt_region = .true.
      nchunks = 6
      xi_width = 90.
      eta_width = 90.
      read(ch_arg7,*) lat0
    else
      read(ch_arg7,*) nchunks
      read(ch_arg8,*) xi_width
      call getarg(9,ch_eta_width)
      call getarg(10,ch_clon)
      call getarg(11,ch_clat)
      call getarg(12,ch_grot)
      if (trim(ch_eta_width) == '' .or. trim(ch_clon) == '' .or. trim(ch_clat) == '' .or. trim(ch_grot) == '') then
        print *, trim(ch_usage)
        stop
      else
        read(ch_eta_width,*) eta_width
        read(ch_clon,*) clon
        read(ch_clat,*) clat
        read(ch_grot,*) grot

        call getarg(13,ch_arg13)
        if (trim(ch_arg13) == '') then
          belt_region = .false.
        else
          belt_region = .true.
          read(ch_arg13,*) lat0
        endif
      endif
    endif
  endif

  ! check the input arguements
  if (nchunks /= 6 .and. nchunks /= 1 .and. nchunks /= 2) stop 'nchunks = 1, 2, 6'
  if (narc /= 1 .and. narc /= 0) stop 'narc = 0,1'
  if (xi_width > 90.0 .or. eta_width > 90.0) stop 'xi(eta)_width <= 90.0'

  ! convert to degrees to radians
  ths = (90-elat)/180. * pi
  phs = elon/180. * pi
  thr = (90-slat)/180. * pi
  phr = slon/180. * pi
  xi_width = xi_width/180. * pi
  eta_width = eta_width/180. * pi
  if (nchunks /= 6) then
    clon = clon/180. * pi
    clat = clat/180. * pi
    grot = grot/180. * pi
  endif

  if (belt_region) th0 = lat0/180. * pi

  ! figure out the normal of source-receiver plane, and rotate them to the equator
  call tp2norm(ths,phs,thr,phr,thn,phn)
  call norm_rot(thn,phn,ths,phs,ths_new,phs_new)
  call norm_rot(thn,phn,thr,phr,thr_new,phr_new)
  if ((ths_new-pi/2) > EPS .or. (thr_new-pi/2) > EPS) stop 'New lat of source/receiver is not 0'

  ! discretize the source receiver garc, and figure out the slice number for individual point
  np = 180
  if (phr_new < phs_new) phr_new = phr_new + 2 * pi
  if (((phr_new - phs_new) <= pi .and. narc == 0) .or. ((phr_new-phs_new) > pi .and. narc == 1)) then
    call linarray(ph,np,phs_new,phr_new)
    delta = phr_new - phs_new
  else
    call linarray(ph,np,phr_new-2*pi,phs_new)
    delta = phs_new - phr_new + 2*pi
  endif
  !print *, 'Delta = ', delta

  if (nchunks /= 6) call rotmat(clon,clat,grot,rotation_matrix)

  chunk(1:np) = 10.

  if (belt_region) then
    nt = 10
    call linarray(th,nt,-th0,th0)
  else
    nt = 1
    th(1) = 0.
    open(11,file='garc_station.txt',status='unknown')
    open(12,file='xsection_translate.txt',status='unknown')
  endif

  j = 0
  do ii = 1, nt
    do i = 1, np
      call norm_rot_back(thn,phn,pi/2+th(ii),ph(i),thp(i),php(i))
      call tp2xyz(thp(i),php(i),x,y,z)

      if (.not. belt_region) then
        write(11,*) thp(i) * 180./pi, php(i) * 180./pi
        if (i == np/2) write(12,*) x, y, z
      endif

      if (nchunks /= 6) then
        xn = x * rotation_matrix(1,1) + y * rotation_matrix(2,1) + z * rotation_matrix(3,1)
        yn = x * rotation_matrix(1,2) + y * rotation_matrix(2,2) + z * rotation_matrix(3,2)
        zn = x * rotation_matrix(1,3) + y * rotation_matrix(2,3) + z * rotation_matrix(3,3)
        x = xn; y = yn; z = zn
      endif

      do k = 1, nchunks
        call chunk_map(k,x,y,z,xik,etak)
        if (abs(xik) <= pi/4 .and. abs(etak) <= pi/4) then
          chunk(i) = k; xi(i) = xik; eta(i) = etak
        endif
      enddo

      if (nchunks /= 6 .and. (chunk(i) > nchunks .or. abs(xi(i)) > pi/2 .or. abs(eta(i)) > pi/2)) then
        write(ch_usage,*) 'Check if point ', x, ',', y, ',', z, ' is within the region or not'
        print *, trim(ch_usage)
        stop
      endif

      xi1 = xi(i)/xi_width * 2; eta1 = eta(i) /eta_width * 2
      nproc_xi = floor((xi1+1)/2 * nproc);  nproc_eta = floor((eta1+1)/2 * nproc)
      j = j + 1; slices(j) = nproc * nproc * (chunk(i)-1) + nproc * nproc_eta + nproc_xi
    enddo
  enddo

  ! compact the slice numbers into an array
  np = np * nt
  call compact_int_array(slices,np,slices_short,np_short)

  ! write the result out
  do i = 1, np_short
    print *, slices_short(i)
  enddo

  if (.not. belt_region) close(11); close(12)


end program global_slice_number
