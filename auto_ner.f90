

  subroutine auto_ner(WIDTH, NEX_MAX, &
       NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
       NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
       NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
       R_CENTRAL_CUBE, CASE_3D)
    
    implicit none

    include 'constants.h'

    double precision WIDTH
    integer NEX_MAX
    integer NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
         NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
         NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB
    double precision R_CENTRAL_CUBE
    logical CASE_3D

    integer,          parameter                :: NUM_REGIONS = 14
    integer,          dimension(NUM_REGIONS)   :: scaling
    double precision, dimension(NUM_REGIONS)   :: radius
    double precision, dimension(NUM_REGIONS-1) :: ratio_top
    double precision, dimension(NUM_REGIONS-1) :: ratio_bottom
    integer,          dimension(NUM_REGIONS-1) :: NER

    ! This is PREM in Kilometers, well ... kinda, not really ....
    radius(1)  = 6371.00d0 ! Surface
    radius(2)  = 6346.60d0 !    Moho - 1st Mesh Doubling Interface
    radius(3)  = 6291.60d0 !      80
    radius(4)  = 6151.00d0 !     220
    radius(5)  = 5971.00d0 !     400
    radius(6)  = 5771.00d0 !     600
    radius(7)  = 5701.00d0 !     670
    radius(8)  = 5600.00d0 !     771
    radius(9)  = 4712.00d0 !    1650 - 2nd Mesh Doubling: Geochemical Layering; Kellogg et al. 1999, Science
    radius(10) = 3630.00d0 !     D''
    radius(11) = 3480.00d0 !     CMB
    radius(12) = 2511.00d0 !    3860 - 3rd Mesh Doubling Interface
    radius(13) = 1371.00d0 !    5000 - 4th Mesh Doubling Interface
    radius(14) =  982.00d0 ! Top Central Cube


    ! Mesh Doubling
    scaling(1)     = 1  ! SURFACE TO MOHO
    scaling(2:8)   = 2  ! MOHO    TO G'' (Geochemical Mantle 1650)
    scaling(9:11)  = 4  ! G''     TO MIC (Middle Inner Core)
    scaling(12)    = 8  ! MIC     TO MIC-II
    scaling(13:14) = 16 ! MIC-II  TO Central Cube -> Center of the Earth

    ! Minimum Number of Elements a Region must have
    NER(:)    = 1
    NER(3:5)  = 2
    if(CASE_3D) then
       NER(1) = 2
    endif

    ! Find the Number of Radial Elements in a region based upon
    ! the aspect ratio of the elements
    call auto_optimal_ner(NUM_REGIONS, WIDTH, NEX_MAX, radius, scaling, NER, ratio_top, ratio_bottom)

    ! Set Output arguments
    NER_CRUST                = NER(1)
    NER_80_MOHO              = NER(2)
    NER_220_80               = NER(3)
    NER_400_220              = NER(4)
    NER_600_400              = NER(5)
    NER_670_600              = NER(6)
    NER_771_670              = NER(7)
    NER_TOPDDOUBLEPRIME_771  = NER(8) + NER(9)
    NER_CMB_TOPDDOUBLEPRIME  = NER(10)
    NER_OUTER_CORE           = NER(11) + NER(12)
    NER_TOP_CENTRAL_CUBE_ICB = NER(13)
    R_CENTRAL_CUBE           = 950000.d0

  end subroutine auto_ner

  subroutine auto_optimal_ner(NUM_REGIONS, width, NEX, r, scaling, NER, rt, rb)

    implicit none

    include 'constants.h'

    integer NUM_REGIONS
    integer NEX
    double precision  width                                ! Width of the Chunk in Degrees
    integer,          dimension(NUM_REGIONS-1) :: NER      ! Elements per Region    - IN-N-OUT - Yummy !
    integer,          dimension(NUM_REGIONS)   :: scaling  ! Element Doubling       - INPUT
    double precision, dimension(NUM_REGIONS)   :: r        ! Radius                 - INPUT
    double precision, dimension(NUM_REGIONS-1) :: rt       ! Ratio at Top           - OUTPUT
    double precision, dimension(NUM_REGIONS-1) :: rb       ! Ratio at Bottom        - OUTPUT

    double precision dr, w, ratio, xi, ximin, wt, wb
    integer ner_test
    integer i

    ! Find optimal elements per region
    do i = 1,NUM_REGIONS-1
       dr = r(i) - r(i+1)              ! Radial Length of Ragion
       wt = width * PI/180.0d0 * r(i)   / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Top
       wb = width * PI/180.0d0 * r(i+1) / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Bottom
       w  = (wt + wb) * 0.5d0          ! Average Width of Region
       ner_test = NER(i)               ! Initial solution
       ratio = (dr / ner_test) / w     ! Aspect Ratio of Element
       xi = dabs(ratio - 1.0d0)        ! Aspect Ratio should be near 1.0
       ximin = 1e7                     ! Initial Minimum

       do while(xi <= ximin)
          NER(i) = ner_test            ! Found a better solution
          ximin = xi                   !
          ner_test = ner_test + 1      ! Increment ner_test and
          ratio = (dr / ner_test) / w  ! look for a better
          xi = dabs(ratio - 1.0d0)     ! solution
       end do
       rt(i) = dr / NER(i) / wt        ! Find the Ratio of Top
       rb(i) = dr / NER(i) / wb        ! and Bottom for completeness
    end do

  end subroutine auto_optimal_ner

  subroutine find_r_central_cube(nex_xi_in, rcube, nex_eta_in)
    implicit none

    integer, parameter :: NBNODE = 8
    double precision, parameter :: alpha = 0.41d0

    integer i, npts
    integer nex_xi, nex_eta_in, nex_xi_in
    integer nex_eta
    double precision rcube, rcubestep, rcube_test, rcubemin, rcubemax
    double precision xi, ximin
    double precision , allocatable, dimension(:,:) :: points
    double precision elem(NBNODE+1, 2)
    integer nspec_cube, nspec_chunks, ispec, nspec
    double precision edgemax, edgemin
    double precision max_edgemax, min_edgemin
    double precision aspect_ratio, max_aspect_ratio

    nex_xi = nex_xi_in / 16


    rcubestep    = 1.0d0
    rcube_test   = 930.0d0
    rcubemax     = 1100.0d0
    nex_eta_in   = -1
    ximin        = 1e7
    rcube        = rcube_test
    
    do while(rcube_test <= rcubemax)
       max_edgemax = -1e7
       min_edgemin = 1e7
       max_aspect_ratio = 0.0d0
       call compute_nex(nex_xi, rcube_test, alpha, nex_eta)
       npts = (4 * nex_xi * nex_eta * NBNODE) + (nex_xi * nex_xi * NBNODE)
       allocate(points(npts, 2))
       call compute_IC_mesh(rcube_test, points, npts, nspec_cube, nspec_chunks, nex_xi, nex_eta)
       nspec = nspec_cube + nspec_chunks
       do ispec = 1,nspec
          call get_element(points, ispec, npts, elem)
          call get_size_min_max(elem, edgemax, edgemin)
          aspect_ratio = edgemax / edgemin
          max_edgemax = MAX(max_edgemax, edgemax)
          min_edgemin = MIN(min_edgemin, edgemin)
          max_aspect_ratio = MAX(max_aspect_ratio, aspect_ratio)
       end do
       xi = (max_edgemax / min_edgemin) 
!       xi = abs(rcube_test - 981.0d0) / 45.0d0
!       write(*,'(a,5(f14.4,2x))')'rcube, xi, ximin:-',rcube_test, xi, min_edgemin,max_edgemax,max_aspect_ratio
       deallocate(points)
       if(xi < ximin) then
          ximin      = xi
          rcube      = rcube_test
          nex_eta_in = nex_eta
       endif
       rcube_test = rcube_test + rcubestep
    enddo

  end subroutine find_r_central_cube

  subroutine compute_nex(nex_xi, rcube, alpha, ner) 
    implicit none

    double precision, parameter :: RICB_KM = 1221.0d0
    double precision, parameter :: PI = 3.1415

    integer nex_xi, ner
    double precision rcube, alpha
    integer ix
    double precision ratio_x, factx, xi
    double precision x, y
    double precision surfx, surfy
    double precision dist_cc_icb, somme, dist_moy

    somme = 0.0d0

    do ix = 0,nex_xi/2,1
       ratio_x = (ix * 1.0d0) / ( nex_xi * 1.0d0)
       factx = 2.0d0 * ratio_x - 1.0d0
       xi = (PI / 2.0d0) * factx
       x = (rcube / sqrt(2.0d0)) * factx
       y = (rcube / sqrt(2.0d0)) * (1 + cos(xi) * alpha / (PI / 2.0d0))

       surfx = RICB_KM * cos(3 * (PI/4.0d0) - ratio_x * (PI/2.0d0))
       surfy = RICB_KM * sin(3 * (PI/4.0d0) - ratio_x * (PI/2.0d0))

       dist_cc_icb = sqrt((surfx -x)**2 + (surfy - y)**2)
       if(ix /= nex_xi/2) then
          dist_cc_icb = dist_cc_icb * 2
       endif
       somme = somme + dist_cc_icb
    end do
    dist_moy = somme / (nex_xi + 1)
    ner = nint(dist_moy / ((PI * RICB_KM) / (2*nex_xi)))
  end subroutine compute_nex

  subroutine get_element(points, ispec, npts, pts)
    implicit none
    integer i, npts, ispec
    integer, parameter :: NBNODE = 8
    double precision pts(NBNODE+1,2), points(npts,2)
    pts(1:8,:) = points( ( (ispec-1) * NBNODE)+1 : ( (ispec) * NBNODE )+1, : )
    pts(NBNODE+1,:) = pts(1,:)  ! Use first point as the last point
  end subroutine get_element

  subroutine get_size_min_max(pts, edgemax, edgemin)
    implicit none
    integer ie, ix1,ix2,ix3
    integer, parameter :: NBNODE = 8
    double precision edgemax, edgemin, edge
    double precision pts(NBNODE+1, 2)


    edgemax = -1e7
    edgemin = -edgemax
    do ie = 1,NBNODE/2,1
        ix1 = (ie * 2) - 1
        ix2 = ix1 + 1
        ix3 = ix1 + 2
        edge = sqrt( (pts(ix1,1) - pts(ix2,1))**2 + (pts(ix1,2) - pts(ix2,2))**2 ) + &
               sqrt( (pts(ix2,1) - pts(ix3,1))**2 + (pts(ix2,2) - pts(ix3,2))**2 )
        edgemax = MAX(edgemax, edge)
        edgemin = MIN(edgemin, edge)
    end do
  end subroutine get_size_min_max

  subroutine compute_IC_mesh(rcube, points, npts, nspec_cube, nspec_chunks, nex_xi, nex_eta)
    implicit none

    integer, parameter :: NBNODE = 8
    integer npts
    integer nspec_chunks, nspec_cube
    double precision rcube
    double precision alpha 
    double precision points(npts, 2)
    double precision x, y
    integer kmax
    integer nex_eta, nex_xi
    integer ic, ix, iy, in
    integer, parameter, dimension(NBNODE) :: iaddx(NBNODE) = (/0,1,2,2,2,1,0,0/)
    integer, parameter, dimension(NBNODE) :: iaddy(NBNODE) = (/0,0,0,1,2,2,2,1/)
    integer k

    k = 1
    alpha = 0.41d0
    nspec_chunks = 0
    do ic = 0,3
       do ix = 0,(nex_xi-1)*2,2
          do iy = 0,(nex_eta-1)*2,2
             do in = 1,NBNODE
                call compute_coordinate(ix+iaddx(in), iy+iaddy(in), nex_xi*2, nex_eta*2, rcube, ic, alpha, x,y)
                points(k,1) = x
                points(k,2) = y
                k = k + 1
             end do
             nspec_chunks = nspec_chunks + 1
          end do
       end do
    end do

    nspec_cube = 0
    do ix = 0,(nex_xi-1)*2,2
       do iy = 0,(nex_xi-1)*2,2
          do in = 1,NBNODE
             call compute_coordinate_central_cube(ix+iaddx(in), iy+iaddy(in), nex_xi*2, nex_xi*2, rcube, alpha,x,y)
             points(k,1) = x
             points(k,2) = y
             k = k + 1
          end do
          nspec_cube = nspec_cube + 1
       end do
    end do
    
  end subroutine compute_IC_mesh

  subroutine compute_coordinate_central_cube(ix,iy,nbx,nby,radius, alpha, x, y)
    implicit none

    double precision, parameter :: PI = 3.1415d0

    integer ix, iy, nbx, nby
    double precision radius, alpha
    double precision x, y

    double precision ratio_x, ratio_y
    double precision factx, facty
    double precision xi, eta

    ratio_x = (ix * 1.0d0) / (nbx * 1.0d0)
    ratio_y = (iy * 1.0d0) / (nby * 1.0d0)

    factx = 2.0d0 * ratio_x - 1.0d0
    facty = 2.0d0 * ratio_y - 1.0d0
    
    xi  = (PI / 2.0d0) * factx
    eta = (PI / 2.0d0) * facty

    x = (radius / sqrt(2.0d0)) * factx * ( 1 + cos(eta) * alpha / (PI / 2.0d0))
    y = (radius / sqrt(2.0d0)) * facty * ( 1 + cos(xi)  * alpha / (PI / 2.0d0))

  end subroutine compute_coordinate_central_cube

  subroutine compute_coordinate(ix,iy,nbx, nby, rcube, ic, alpha, x, y)
    implicit none

    double precision, parameter :: PI      = 3.1415d0
    double precision, parameter :: RICB_KM = 1221.0d0

    integer ix, iy, nbx, nby, ic
    double precision rcube, alpha
    double precision x, y

    double precision ratio_x, ratio_y
    double precision factx, xi
    double precision xcc, ycc
    double precision xsurf, ysurf
    double precision deltax, deltay
    double precision temp

    ratio_x = (ix * 1.0d0) / (nbx * 1.0d0)
    ratio_y = (iy * 1.0d0) / (nby * 1.0d0)
    
    factx = 2.0d0 * ratio_x - 1.0d0
    xi = (PI/2.0d0) * factx

    xcc = (rcube / sqrt(2.0d0)) * factx
    ycc = (rcube / sqrt(2.0d0)) * (1 + cos(xi) * alpha / (PI/2.0d0))
    
    xsurf = RICB_KM * cos(3.0d0 * (PI/4.0d0) - ratio_x * (PI/2.0d0))
    ysurf = RICB_KM * sin(3.0d0 * (PI/4.0d0) - ratio_x * (PI/2.0d0))

    deltax = xsurf - xcc
    deltay = ysurf - ycc

    x = xsurf - ratio_y * deltax
    y = ysurf - ratio_y * deltay

    if(ic == 1) then
       temp = x
       x    = y
       y    = temp
    else if (ic == 2) then
       x = -x
       y = -y
    else if (ic == 3) then
       temp = x
       x    = -y
       y    = temp
    end if
  end subroutine compute_coordinate
