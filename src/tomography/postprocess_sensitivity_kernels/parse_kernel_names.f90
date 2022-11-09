!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  use postprocess_par, only: MAX_STRING_LEN, MAX_KERNEL_NAMES

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN),dimension(MAX_KERNEL_NAMES),intent(out) :: kernel_names
  integer,intent(out) :: nker

  ! local parameters
  integer :: iker
  character(len=MAX_STRING_LEN) :: tmp,kernel_name
  character(len=1),parameter :: delimiter = ','

  ! safety check
  if (len_trim(kernel_names_comma_delimited) == 0) then
    print *,'Error: empty kernel name as first entry ',kernel_names_comma_delimited
    print *,'Please provide a kernel name'
    stop 'Error empty kernel name string'
  endif

  ! gets first name/token
  iker = 1
  call tokenize_string(kernel_names_comma_delimited, delimiter, kernel_name)
  kernel_names(iker) = kernel_name

  ! null-string as first argument for successive tokenize_string-calls
  tmp(1:1) = char(0)

  ! gets next names/tokens (tokenize_string will return null-terminated token when finished)
  do while (kernel_names(iker)(1:1) /= char(0))
    ! increases name/token number
    iker = iker + 1
    if (iker > MAX_KERNEL_NAMES) stop 'number of tokens exceeds MAX_KERNEL_NAMES'

    ! gets next successive token (with null-terminated string as first argument)
    call tokenize_string(tmp, delimiter, kernel_name)
    kernel_names(iker) = kernel_name
  enddo

  ! number of kernel names
  nker = iker-1

  ! checks name lengths (e.g. if kernel_name argument is "vsv,vsh," we will have a 3. kernel name with empty string)
  do iker = 1,nker
    if (len_trim(kernel_names(iker)) == 0) then
      print *,'Error encountered kernel name with zero length: kernel name number ',iker,' out of ',nker,' is empty'
      print *,'Please check your kernel_names argument...'
      stop 'Error kernel name with zero length'
    endif
  enddo

  end subroutine parse_kernel_names


!
!-------------------------------------------------------------------------------------------------
!

! The following utility function was modified from http://Fortranwiki.org/Fortran/show/strtok
!

! note: IBM's xlf compiler seems to have an issue with a subroutine name "strtok".
!    we're using a more Fortran-like name.
  subroutine tokenize_string (source_string, delimiter, token)

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in delimiter.
!
!             then, if the returned value is not equal to char(0), keep calling until it is
!             with SOURCE_STRING set to char(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning char(0).
!
!     Input:  source_string =   Source string to tokenize.
!             delimiter    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures

  use postprocess_par, only: MAX_STRING_LEN

  !     PARAMETERS:
  character(len=MAX_STRING_LEN), intent(in)  :: source_string
  character(len=1), intent(in)  :: delimiter
  character(len=MAX_STRING_LEN), intent(out) :: token

  !     SAVED VALUES:
  character(len=MAX_STRING_LEN),save :: saved_string
  integer,save :: isaved_start  ! points to beginning of unprocessed data
  integer,save :: isource_len   ! length of original input string

  !     LOCAL VALUES:
  integer :: ibegin        ! beginning of token to return
  integer :: ifinish       ! end of token to return
  integer :: length

  ! initialize stored copy of input string and pointer into input string on first call
  if (source_string(1:1) /= char(0)) then
    isaved_start = 1                      ! beginning of unprocessed data
    saved_string = trim(source_string)    ! save input string from first call in series
    isource_len = len_trim(saved_string)  ! length of input string from first call
  endif

  ! token = '' is not allowed by IBM's xlf compiler. assigning empty strings is non-standard Fortran.
  ! as a work-around, we set a null character '' here to initialize the string
  token = char(0)
  ibegin = isaved_start

  ! sets first index ibegin to beginning of (next) token
  do while (.true.)
    if ( (ibegin <= isource_len) .and. (index(delimiter,saved_string(ibegin:ibegin)) /= 0)) then
      ! delimiter is encountered, starts with next index (next token)
      ibegin = ibegin + 1
    else
      ! exits do-loop
      exit
    endif
  enddo

  if (ibegin > isource_len) then
    token(1:1) = char(0)
    return
  endif

  ! sets second index ifinish to end of token (including delimiter)
  ifinish = ibegin

  do while (.true.)
    if ((ifinish <= isource_len) .and. (index(delimiter,saved_string(ifinish:ifinish)) == 0)) then
      ! delimiter is not encountered yet, increases finish index
      ifinish = ifinish + 1
    else
      ! exits do-loop
      exit
    endif
    ! avoids accessing saved_string in next if-statement out of bounds
    if (ifinish > isource_len) exit
  enddo

  ! sets token string
  !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
  ! IBM's xlf compiler doesn't like assigining strings of zero lengths?
  length = ifinish-1 - ibegin + 1
  if (length > 0 .and. length < MAX_STRING_LEN) then
    token = saved_string(ibegin:ifinish-1)
  else
    stop 'invalid token, either too short or too long'
  endif

  isaved_start = ifinish

  end subroutine tokenize_string

