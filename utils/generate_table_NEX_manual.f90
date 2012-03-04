
  program generate_table_NEX_manual

! generate Table showing possible values of NEX for the SPECFEM3D_GLOBE manual

! Author: Dimitri Komatitsch, University of Pau, France, September 2007

  implicit none

  integer NPROC_XI,c,NEX_XI,counter

  integer, parameter :: NB_COLUMNS_TABLE = 10
  integer, dimension(NB_COLUMNS_TABLE) :: store_NEX_XI

! maximum total number of processors we want to see in the table
  integer, parameter :: MAX_NUMBER_OF_PROCS = 100000
! integer, parameter :: MAX_NUMBER_OF_PROCS = 212992 ! current maximum on BlueGene at LLNL

! base value depends if we implement three or four doublings (default is three)
  integer, parameter :: NB_DOUBLING = 3
  integer :: BASE_VALUE

! output in LaTeX format or in regular ASCII format
  logical, parameter :: OUTPUT_LATEX_FORMAT = .false.

  BASE_VALUE = 2**NB_DOUBLING

! total number of processors is 6 * NPROC_XI^2
  do NPROC_XI = 1,int(sqrt(MAX_NUMBER_OF_PROCS / 6.d0))

    counter = 1
    c = 0

    do while (counter <= NB_COLUMNS_TABLE)
      c = c + 1
      NEX_XI = BASE_VALUE * c * NPROC_XI
      if(NEX_XI >= 64 .and. mod(NEX_XI,2*BASE_VALUE) == 0) then
        store_NEX_XI(counter) = NEX_XI
        counter = counter + 1
      endif
    enddo

    if(OUTPUT_LATEX_FORMAT) then
      write(*,"(i6,' &')") NPROC_XI
      write(*,"(i6,' &')") 6*NPROC_XI**2
      do counter = 1,NB_COLUMNS_TABLE-1
        write(*,"(i6,' &')") store_NEX_XI(counter)
      enddo
      write(*,*) store_NEX_XI(NB_COLUMNS_TABLE),' ',achar(92),'tabularnewline'
    else
      write(*,"(i6,' | ',i6,' | ',i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
        NPROC_XI, 6*NPROC_XI**2,(store_NEX_XI(counter), counter = 1,NB_COLUMNS_TABLE)
    endif
  enddo

  end program generate_table_NEX_manual

