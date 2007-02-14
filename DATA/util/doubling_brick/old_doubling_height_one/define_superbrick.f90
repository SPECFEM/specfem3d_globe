
  subroutine define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick)

! define the symmetric mesh doubling superbrick composed of 28 elements and 58 points

  include "constants_modified.h"

  integer, dimension(NGNOD_DOUBLING_SUPERBRICK,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

  x_superbrick( 1) = 3.d0 / 2.d0
  y_superbrick( 1) = 1.d0
  z_superbrick( 1) = 1.d0

  x_superbrick( 2) = 3.d0 / 2.d0
  y_superbrick( 2) = 1.d0
  z_superbrick( 2) = 2.d0 / 3.d0

  x_superbrick( 3) = 3.d0 / 2.d0
  y_superbrick( 3) = 3.d0 / 2.d0
  z_superbrick( 3) = 2.d0 / 3.d0

  x_superbrick( 4) = 3.d0 / 2.d0
  y_superbrick( 4) = 3.d0 / 2.d0
  z_superbrick( 4) = 1.d0

  x_superbrick( 5) = 2.d0
  y_superbrick( 5) = 1.d0
  z_superbrick( 5) = 1.d0

  x_superbrick( 6) = 2.d0
  y_superbrick( 6) = 1.d0
  z_superbrick( 6) = 1.d0 / 3.d0

  x_superbrick( 7) = 2.d0
  y_superbrick( 7) = 3.d0 / 2.d0
  z_superbrick( 7) = 1.d0 / 3.d0

  x_superbrick( 8) = 2.d0
  y_superbrick( 8) = 3.d0 / 2.d0
  z_superbrick( 8) = 1.d0

  x_superbrick( 9) = 3.d0 / 2.d0
  y_superbrick( 9) = 2.d0
  z_superbrick( 9) = 1.d0 / 3.d0

  x_superbrick(10) = 3.d0 / 2.d0
  y_superbrick(10) = 2.d0
  z_superbrick(10) = 1.d0

  x_superbrick(11) = 2.d0
  y_superbrick(11) = 2.d0
  z_superbrick(11) = 0.d0

  x_superbrick(12) = 2.d0
  y_superbrick(12) = 2.d0
  z_superbrick(12) = 1.d0

  x_superbrick(13) = 1.d0
  y_superbrick(13) = 1.d0
  z_superbrick(13) = 1.d0 / 3.d0

  x_superbrick(14) = 1.d0
  y_superbrick(14) = 1.d0
  z_superbrick(14) = 0.d0

  x_superbrick(15) = 1.d0
  y_superbrick(15) = 2.d0
  z_superbrick(15) = 0.d0

  x_superbrick(16) = 1.d0
  y_superbrick(16) = 2.d0
  z_superbrick(16) = 1.d0 / 3.d0

  x_superbrick(17) = 3.d0 / 2.d0
  y_superbrick(17) = 1.d0
  z_superbrick(17) = 1.d0 / 3.d0

  x_superbrick(18) = 2.d0
  y_superbrick(18) = 1.d0
  z_superbrick(18) = 0.d0

  x_superbrick(19) = 1.d0
  y_superbrick(19) = 1.d0
  z_superbrick(19) = 2.d0 / 3.d0

  x_superbrick(20) = 1.d0
  y_superbrick(20) = 1.d0
  z_superbrick(20) = 1.d0

  x_superbrick(21) = 1.d0
  y_superbrick(21) = 3.d0 / 2.d0
  z_superbrick(21) = 2.d0 / 3.d0

  x_superbrick(22) = 1.d0
  y_superbrick(22) = 3.d0 / 2.d0
  z_superbrick(22) = 1.d0

  x_superbrick(23) = 1.d0
  y_superbrick(23) = 2.d0
  z_superbrick(23) = 1.d0

  x_superbrick(24) = 3.d0 / 2.d0
  y_superbrick(24) = 1.d0 / 2.d0
  z_superbrick(24) = 2.d0 / 3.d0

  x_superbrick(25) = 3.d0 / 2.d0
  y_superbrick(25) = 1.d0 / 2.d0
  z_superbrick(25) = 1.d0

  x_superbrick(26) = 2.d0
  y_superbrick(26) = 1.d0 / 2.d0
  z_superbrick(26) = 1.d0 / 3.d0

  x_superbrick(27) = 2.d0
  y_superbrick(27) = 1.d0 / 2.d0
  z_superbrick(27) = 1.d0

  x_superbrick(28) = 3.d0 / 2.d0
  y_superbrick(28) = 0.d0
  z_superbrick(28) = 1.d0 / 3.d0

  x_superbrick(29) = 3.d0 / 2.d0
  y_superbrick(29) = 0.d0
  z_superbrick(29) = 1.d0

  x_superbrick(30) = 2.d0
  y_superbrick(30) = 0.d0
  z_superbrick(30) = 0.d0

  x_superbrick(31) = 2.d0
  y_superbrick(31) = 0.d0
  z_superbrick(31) = 1.d0

  x_superbrick(32) = 1.d0
  y_superbrick(32) = 0.d0
  z_superbrick(32) = 0.d0

  x_superbrick(33) = 1.d0
  y_superbrick(33) = 0.d0
  z_superbrick(33) = 1.d0 / 3.d0

  x_superbrick(34) = 1.d0
  y_superbrick(34) = 1.d0 / 2.d0
  z_superbrick(34) = 2.d0 / 3.d0

  x_superbrick(35) = 1.d0
  y_superbrick(35) = 1.d0 / 2.d0
  z_superbrick(35) = 1.d0

  x_superbrick(36) = 1.d0
  y_superbrick(36) = 0.d0
  z_superbrick(36) = 1.d0

  x_superbrick(37) = 1.d0 / 2.d0
  y_superbrick(37) = 1.d0
  z_superbrick(37) = 1.d0

  x_superbrick(38) = 1.d0 / 2.d0
  y_superbrick(38) = 1.d0
  z_superbrick(38) = 2.d0 / 3.d0

  x_superbrick(39) = 1.d0 / 2.d0
  y_superbrick(39) = 3.d0 / 2.d0
  z_superbrick(39) = 2.d0 / 3.d0

  x_superbrick(40) = 1.d0 / 2.d0
  y_superbrick(40) = 3.d0 / 2.d0
  z_superbrick(40) = 1.d0

  x_superbrick(41) = 0.d0
  y_superbrick(41) = 1.d0
  z_superbrick(41) = 1.d0

  x_superbrick(42) = 0.d0
  y_superbrick(42) = 1.d0
  z_superbrick(42) = 1.d0 / 3.d0

  x_superbrick(43) = 0.d0
  y_superbrick(43) = 3.d0 / 2.d0
  z_superbrick(43) = 1.d0 / 3.d0

  x_superbrick(44) = 0.d0
  y_superbrick(44) = 3.d0 / 2.d0
  z_superbrick(44) = 1.d0

  x_superbrick(45) = 1.d0 / 2.d0
  y_superbrick(45) = 2.d0
  z_superbrick(45) = 1.d0 / 3.d0

  x_superbrick(46) = 1.d0 / 2.d0
  y_superbrick(46) = 2.d0
  z_superbrick(46) = 1.d0

  x_superbrick(47) = 0.d0
  y_superbrick(47) = 2.d0
  z_superbrick(47) = 0.d0

  x_superbrick(48) = 0.d0
  y_superbrick(48) = 2.d0
  z_superbrick(48) = 1.d0

  x_superbrick(49) = 1.d0 / 2.d0
  y_superbrick(49) = 1.d0
  z_superbrick(49) = 1.d0 / 3.d0

  x_superbrick(50) = 0.d0
  y_superbrick(50) = 1.d0
  z_superbrick(50) = 0.d0

  x_superbrick(51) = 1.d0 / 2.d0
  y_superbrick(51) = 1.d0 / 2.d0
  z_superbrick(51) = 2.d0 / 3.d0

  x_superbrick(52) = 1.d0 / 2.d0
  y_superbrick(52) = 1.d0 / 2.d0
  z_superbrick(52) = 1.d0

  x_superbrick(53) = 0.d0
  y_superbrick(53) = 1.d0 / 2.d0
  z_superbrick(53) = 1.d0 / 3.d0

  x_superbrick(54) = 0.d0
  y_superbrick(54) = 1.d0 / 2.d0
  z_superbrick(54) = 1.d0

  x_superbrick(55) = 1.d0 / 2.d0
  y_superbrick(55) = 0.d0
  z_superbrick(55) = 1.d0 / 3.d0

  x_superbrick(56) = 1.d0 / 2.d0
  y_superbrick(56) = 0.d0
  z_superbrick(56) = 1.d0

  x_superbrick(57) = 0.d0
  y_superbrick(57) = 0.d0
  z_superbrick(57) = 0.d0

  x_superbrick(58) = 0.d0
  y_superbrick(58) = 0.d0
  z_superbrick(58) = 1.d0

  ibool_superbrick(1, 1) =  1
  ibool_superbrick(2, 1) =  2
  ibool_superbrick(3, 1) =  3
  ibool_superbrick(4, 1) =  4
  ibool_superbrick(5, 1) =  5
  ibool_superbrick(6, 1) =  6
  ibool_superbrick(7, 1) =  7
  ibool_superbrick(8, 1) =  8

  ibool_superbrick(1, 2) =  4
  ibool_superbrick(2, 2) =  3
  ibool_superbrick(3, 2) =  9
  ibool_superbrick(4, 2) = 10
  ibool_superbrick(5, 2) =  8
  ibool_superbrick(6, 2) =  7
  ibool_superbrick(7, 2) = 11
  ibool_superbrick(8, 2) = 12

  ibool_superbrick(1, 3) = 13
  ibool_superbrick(2, 3) = 14
  ibool_superbrick(3, 3) = 15
  ibool_superbrick(4, 3) = 16
  ibool_superbrick(5, 3) = 17
  ibool_superbrick(6, 3) = 18
  ibool_superbrick(7, 3) = 11
  ibool_superbrick(8, 3) =  9

  ibool_superbrick(1, 4) = 20
  ibool_superbrick(2, 4) = 19
  ibool_superbrick(3, 4) = 21
  ibool_superbrick(4, 4) = 22
  ibool_superbrick(5, 4) =  1
  ibool_superbrick(6, 4) =  2
  ibool_superbrick(7, 4) =  3
  ibool_superbrick(8, 4) =  4

  ibool_superbrick(1, 5) =  2
  ibool_superbrick(2, 5) = 17
  ibool_superbrick(3, 5) =  9
  ibool_superbrick(4, 5) =  3
  ibool_superbrick(5, 5) =  6
  ibool_superbrick(6, 5) = 18
  ibool_superbrick(7, 5) = 11
  ibool_superbrick(8, 5) =  7

  ibool_superbrick(1, 6) = 22
  ibool_superbrick(2, 6) = 21
  ibool_superbrick(3, 6) = 16
  ibool_superbrick(4, 6) = 23
  ibool_superbrick(5, 6) =  4
  ibool_superbrick(6, 6) =  3
  ibool_superbrick(7, 6) =  9
  ibool_superbrick(8, 6) = 10

  ibool_superbrick(1, 7) = 19
  ibool_superbrick(2, 7) = 13
  ibool_superbrick(3, 7) = 16
  ibool_superbrick(4, 7) = 21
  ibool_superbrick(5, 7) =  2
  ibool_superbrick(6, 7) = 17
  ibool_superbrick(7, 7) =  9
  ibool_superbrick(8, 7) =  3

  ibool_superbrick(1, 8) =  1
  ibool_superbrick(2, 8) =  2
  ibool_superbrick(3, 8) = 24
  ibool_superbrick(4, 8) = 25
  ibool_superbrick(5, 8) =  5
  ibool_superbrick(6, 8) =  6
  ibool_superbrick(7, 8) = 26
  ibool_superbrick(8, 8) = 27

  ibool_superbrick(1, 9) = 25
  ibool_superbrick(2, 9) = 24
  ibool_superbrick(3, 9) = 28
  ibool_superbrick(4, 9) = 29
  ibool_superbrick(5, 9) = 27
  ibool_superbrick(6, 9) = 26
  ibool_superbrick(7, 9) = 30
  ibool_superbrick(8, 9) = 31

  ibool_superbrick(1,10) = 13
  ibool_superbrick(2,10) = 14
  ibool_superbrick(3,10) = 32
  ibool_superbrick(4,10) = 33
  ibool_superbrick(5,10) = 17
  ibool_superbrick(6,10) = 18
  ibool_superbrick(7,10) = 30
  ibool_superbrick(8,10) = 28

  ibool_superbrick(1,11) = 20
  ibool_superbrick(2,11) = 19
  ibool_superbrick(3,11) = 34
  ibool_superbrick(4,11) = 35
  ibool_superbrick(5,11) =  1
  ibool_superbrick(6,11) =  2
  ibool_superbrick(7,11) = 24
  ibool_superbrick(8,11) = 25

  ibool_superbrick(1,12) =  2
  ibool_superbrick(2,12) = 17
  ibool_superbrick(3,12) = 28
  ibool_superbrick(4,12) = 24
  ibool_superbrick(5,12) =  6
  ibool_superbrick(6,12) = 18
  ibool_superbrick(7,12) = 30
  ibool_superbrick(8,12) = 26

  ibool_superbrick(1,13) = 35
  ibool_superbrick(2,13) = 34
  ibool_superbrick(3,13) = 33
  ibool_superbrick(4,13) = 36
  ibool_superbrick(5,13) = 25
  ibool_superbrick(6,13) = 24
  ibool_superbrick(7,13) = 28
  ibool_superbrick(8,13) = 29

  ibool_superbrick(1,14) = 19
  ibool_superbrick(2,14) = 13
  ibool_superbrick(3,14) = 33
  ibool_superbrick(4,14) = 34
  ibool_superbrick(5,14) =  2
  ibool_superbrick(6,14) = 17
  ibool_superbrick(7,14) = 28
  ibool_superbrick(8,14) = 24

  ibool_superbrick(1,15) = 37
  ibool_superbrick(2,15) = 38
  ibool_superbrick(3,15) = 39
  ibool_superbrick(4,15) = 40
  ibool_superbrick(5,15) = 41
  ibool_superbrick(6,15) = 42
  ibool_superbrick(7,15) = 43
  ibool_superbrick(8,15) = 44

  ibool_superbrick(1,16) = 40
  ibool_superbrick(2,16) = 39
  ibool_superbrick(3,16) = 45
  ibool_superbrick(4,16) = 46
  ibool_superbrick(5,16) = 44
  ibool_superbrick(6,16) = 43
  ibool_superbrick(7,16) = 47
  ibool_superbrick(8,16) = 48

  ibool_superbrick(1,17) = 13
  ibool_superbrick(2,17) = 14
  ibool_superbrick(3,17) = 15
  ibool_superbrick(4,17) = 16
  ibool_superbrick(5,17) = 49
  ibool_superbrick(6,17) = 50
  ibool_superbrick(7,17) = 47
  ibool_superbrick(8,17) = 45

  ibool_superbrick(1,18) = 20
  ibool_superbrick(2,18) = 19
  ibool_superbrick(3,18) = 21
  ibool_superbrick(4,18) = 22
  ibool_superbrick(5,18) = 37
  ibool_superbrick(6,18) = 38
  ibool_superbrick(7,18) = 39
  ibool_superbrick(8,18) = 40

  ibool_superbrick(1,19) = 38
  ibool_superbrick(2,19) = 49
  ibool_superbrick(3,19) = 45
  ibool_superbrick(4,19) = 39
  ibool_superbrick(5,19) = 42
  ibool_superbrick(6,19) = 50
  ibool_superbrick(7,19) = 47
  ibool_superbrick(8,19) = 43

  ibool_superbrick(1,20) = 22
  ibool_superbrick(2,20) = 21
  ibool_superbrick(3,20) = 16
  ibool_superbrick(4,20) = 23
  ibool_superbrick(5,20) = 40
  ibool_superbrick(6,20) = 39
  ibool_superbrick(7,20) = 45
  ibool_superbrick(8,20) = 46

  ibool_superbrick(1,21) = 19
  ibool_superbrick(2,21) = 13
  ibool_superbrick(3,21) = 16
  ibool_superbrick(4,21) = 21
  ibool_superbrick(5,21) = 38
  ibool_superbrick(6,21) = 49
  ibool_superbrick(7,21) = 45
  ibool_superbrick(8,21) = 39

  ibool_superbrick(1,22) = 37
  ibool_superbrick(2,22) = 38
  ibool_superbrick(3,22) = 51
  ibool_superbrick(4,22) = 52
  ibool_superbrick(5,22) = 41
  ibool_superbrick(6,22) = 42
  ibool_superbrick(7,22) = 53
  ibool_superbrick(8,22) = 54

  ibool_superbrick(1,23) = 52
  ibool_superbrick(2,23) = 51
  ibool_superbrick(3,23) = 55
  ibool_superbrick(4,23) = 56
  ibool_superbrick(5,23) = 54
  ibool_superbrick(6,23) = 53
  ibool_superbrick(7,23) = 57
  ibool_superbrick(8,23) = 58

  ibool_superbrick(1,24) = 13
  ibool_superbrick(2,24) = 14
  ibool_superbrick(3,24) = 32
  ibool_superbrick(4,24) = 33
  ibool_superbrick(5,24) = 49
  ibool_superbrick(6,24) = 50
  ibool_superbrick(7,24) = 57
  ibool_superbrick(8,24) = 55

  ibool_superbrick(1,25) = 20
  ibool_superbrick(2,25) = 19
  ibool_superbrick(3,25) = 34
  ibool_superbrick(4,25) = 35
  ibool_superbrick(5,25) = 37
  ibool_superbrick(6,25) = 38
  ibool_superbrick(7,25) = 51
  ibool_superbrick(8,25) = 52

  ibool_superbrick(1,26) = 38
  ibool_superbrick(2,26) = 49
  ibool_superbrick(3,26) = 55
  ibool_superbrick(4,26) = 51
  ibool_superbrick(5,26) = 42
  ibool_superbrick(6,26) = 50
  ibool_superbrick(7,26) = 57
  ibool_superbrick(8,26) = 53

  ibool_superbrick(1,27) = 35
  ibool_superbrick(2,27) = 34
  ibool_superbrick(3,27) = 33
  ibool_superbrick(4,27) = 36
  ibool_superbrick(5,27) = 52
  ibool_superbrick(6,27) = 51
  ibool_superbrick(7,27) = 55
  ibool_superbrick(8,27) = 56

  ibool_superbrick(1,28) = 19
  ibool_superbrick(2,28) = 13
  ibool_superbrick(3,28) = 33
  ibool_superbrick(4,28) = 34
  ibool_superbrick(5,28) = 38
  ibool_superbrick(6,28) = 49
  ibool_superbrick(7,28) = 55
  ibool_superbrick(8,28) = 51

  end subroutine define_superbrick

