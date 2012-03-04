
  subroutine define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick)

  include "constants_modified_v40.h"

  integer, dimension(NGNOD_DOUBLING_SUPERBRICK,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

  x_superbrick( 1) = 3.d0 / 2.d0
  y_superbrick( 1) = 1.d0
  z_superbrick( 1) = 2.d0

  x_superbrick( 2) = 3.d0 / 2.d0
  y_superbrick( 2) = 1.d0
  z_superbrick( 2) = 3.d0 / 2.d0

  x_superbrick( 3) = 3.d0 / 2.d0
  y_superbrick( 3) = 3.d0 / 2.d0
  z_superbrick( 3) = 3.d0 / 2.d0

  x_superbrick( 4) = 3.d0 / 2.d0
  y_superbrick( 4) = 3.d0 / 2.d0
  z_superbrick( 4) = 2.d0

  x_superbrick( 5) = 2.d0
  y_superbrick( 5) = 1.d0
  z_superbrick( 5) = 2.d0

  x_superbrick( 6) = 2.d0
  y_superbrick( 6) = 1.d0
  z_superbrick( 6) = 1.d0

  x_superbrick( 7) = 2.d0
  y_superbrick( 7) = 3.d0 / 2.d0
  z_superbrick( 7) = 1.d0

  x_superbrick( 8) = 2.d0
  y_superbrick( 8) = 3.d0 / 2.d0
  z_superbrick( 8) = 2.d0

  x_superbrick( 9) = 3.d0 / 2.d0
  y_superbrick( 9) = 2.d0
  z_superbrick( 9) = 1.d0

  x_superbrick(10) = 3.d0 / 2.d0
  y_superbrick(10) = 2.d0
  z_superbrick(10) = 2.d0

  x_superbrick(11) = 2.d0
  y_superbrick(11) = 2.d0
  z_superbrick(11) = 1.d0 / 2.d0

  x_superbrick(12) = 2.d0
  y_superbrick(12) = 2.d0
  z_superbrick(12) = 2.d0

  x_superbrick(13) = 1.d0
  y_superbrick(13) = 1.d0
  z_superbrick(13) = 1.d0

  x_superbrick(14) = 1.d0
  y_superbrick(14) = 1.d0
  z_superbrick(14) = 1.d0 / 2.d0

  x_superbrick(15) = 1.d0
  y_superbrick(15) = 2.d0
  z_superbrick(15) = 1.d0 / 2.d0

  x_superbrick(16) = 1.d0
  y_superbrick(16) = 2.d0
  z_superbrick(16) = 1.d0

  x_superbrick(17) = 3.d0 / 2.d0
  y_superbrick(17) = 1.d0
  z_superbrick(17) = 1.d0

  x_superbrick(18) = 2.d0
  y_superbrick(18) = 1.d0
  z_superbrick(18) = 1.d0 / 2.d0

  x_superbrick(19) = 1.d0
  y_superbrick(19) = 1.d0
  z_superbrick(19) = 3.d0 / 2.d0

  x_superbrick(20) = 1.d0
  y_superbrick(20) = 1.d0
  z_superbrick(20) = 2.d0

  x_superbrick(21) = 1.d0
  y_superbrick(21) = 3.d0 / 2.d0
  z_superbrick(21) = 3.d0 / 2.d0

  x_superbrick(22) = 1.d0
  y_superbrick(22) = 3.d0 / 2.d0
  z_superbrick(22) = 2.d0

  x_superbrick(23) = 1.d0
  y_superbrick(23) = 2.d0
  z_superbrick(23) = 2.d0

  x_superbrick(24) = 1.d0
  y_superbrick(24) = 1.d0
  z_superbrick(24) = 0.d0

  x_superbrick(25) = 2.d0
  y_superbrick(25) = 1.d0
  z_superbrick(25) = 0.d0

  x_superbrick(26) = 2.d0
  y_superbrick(26) = 2.d0
  z_superbrick(26) = 0.d0

  x_superbrick(27) = 1.d0
  y_superbrick(27) = 2.d0
  z_superbrick(27) = 0.d0

  x_superbrick(28) = 3.d0 / 2.d0
  y_superbrick(28) = 1.d0 / 2.d0
  z_superbrick(28) = 3.d0 / 2.d0

  x_superbrick(29) = 3.d0 / 2.d0
  y_superbrick(29) = 1.d0 / 2.d0
  z_superbrick(29) = 2.d0

  x_superbrick(30) = 2.d0
  y_superbrick(30) = 1.d0 / 2.d0
  z_superbrick(30) = 1.d0

  x_superbrick(31) = 2.d0
  y_superbrick(31) = 1.d0 / 2.d0
  z_superbrick(31) = 2.d0

  x_superbrick(32) = 3.d0 / 2.d0
  y_superbrick(32) = 0.d0
  z_superbrick(32) = 1.d0

  x_superbrick(33) = 3.d0 / 2.d0
  y_superbrick(33) = 0.d0
  z_superbrick(33) = 2.d0

  x_superbrick(34) = 2.d0
  y_superbrick(34) = 0.d0
  z_superbrick(34) = 1.d0 / 2.d0

  x_superbrick(35) = 2.d0
  y_superbrick(35) = 0.d0
  z_superbrick(35) = 2.d0

  x_superbrick(36) = 1.d0
  y_superbrick(36) = 0.d0
  z_superbrick(36) = 1.d0 / 2.d0

  x_superbrick(37) = 1.d0
  y_superbrick(37) = 0.d0
  z_superbrick(37) = 1.d0

  x_superbrick(38) = 1.d0
  y_superbrick(38) = 1.d0 / 2.d0
  z_superbrick(38) = 3.d0 / 2.d0

  x_superbrick(39) = 1.d0
  y_superbrick(39) = 1.d0 / 2.d0
  z_superbrick(39) = 2.d0

  x_superbrick(40) = 1.d0
  y_superbrick(40) = 0.d0
  z_superbrick(40) = 2.d0

  x_superbrick(41) = 2.d0
  y_superbrick(41) = 0.d0
  z_superbrick(41) = 0.d0

  x_superbrick(42) = 1.d0
  y_superbrick(42) = 0.d0
  z_superbrick(42) = 0.d0

  x_superbrick(43) = 1.d0 / 2.d0
  y_superbrick(43) = 1.d0
  z_superbrick(43) = 2.d0

  x_superbrick(44) = 1.d0 / 2.d0
  y_superbrick(44) = 1.d0
  z_superbrick(44) = 3.d0 / 2.d0

  x_superbrick(45) = 1.d0 / 2.d0
  y_superbrick(45) = 3.d0 / 2.d0
  z_superbrick(45) = 3.d0 / 2.d0

  x_superbrick(46) = 1.d0 / 2.d0
  y_superbrick(46) = 3.d0 / 2.d0
  z_superbrick(46) = 2.d0

  x_superbrick(47) = 0.d0
  y_superbrick(47) = 1.d0
  z_superbrick(47) = 2.d0

  x_superbrick(48) = 0.d0
  y_superbrick(48) = 1.d0
  z_superbrick(48) = 1.d0

  x_superbrick(49) = 0.d0
  y_superbrick(49) = 3.d0 / 2.d0
  z_superbrick(49) = 1.d0

  x_superbrick(50) = 0.d0
  y_superbrick(50) = 3.d0 / 2.d0
  z_superbrick(50) = 2.d0

  x_superbrick(51) = 1.d0 / 2.d0
  y_superbrick(51) = 2.d0
  z_superbrick(51) = 1.d0

  x_superbrick(52) = 1.d0 / 2.d0
  y_superbrick(52) = 2.d0
  z_superbrick(52) = 2.d0

  x_superbrick(53) = 0.d0
  y_superbrick(53) = 2.d0
  z_superbrick(53) = 1.d0 / 2.d0

  x_superbrick(54) = 0.d0
  y_superbrick(54) = 2.d0
  z_superbrick(54) = 2.d0

  x_superbrick(55) = 1.d0 / 2.d0
  y_superbrick(55) = 1.d0
  z_superbrick(55) = 1.d0

  x_superbrick(56) = 0.d0
  y_superbrick(56) = 1.d0
  z_superbrick(56) = 1.d0 / 2.d0

  x_superbrick(57) = 0.d0
  y_superbrick(57) = 1.d0
  z_superbrick(57) = 0.d0

  x_superbrick(58) = 0.d0
  y_superbrick(58) = 2.d0
  z_superbrick(58) = 0.d0

  x_superbrick(59) = 1.d0 / 2.d0
  y_superbrick(59) = 1.d0 / 2.d0
  z_superbrick(59) = 3.d0 / 2.d0

  x_superbrick(60) = 1.d0 / 2.d0
  y_superbrick(60) = 1.d0 / 2.d0
  z_superbrick(60) = 2.d0

  x_superbrick(61) = 0.d0
  y_superbrick(61) = 1.d0 / 2.d0
  z_superbrick(61) = 1.d0

  x_superbrick(62) = 0.d0
  y_superbrick(62) = 1.d0 / 2.d0
  z_superbrick(62) = 2.d0

  x_superbrick(63) = 1.d0 / 2.d0
  y_superbrick(63) = 0.d0
  z_superbrick(63) = 1.d0

  x_superbrick(64) = 1.d0 / 2.d0
  y_superbrick(64) = 0.d0
  z_superbrick(64) = 2.d0

  x_superbrick(65) = 0.d0
  y_superbrick(65) = 0.d0
  z_superbrick(65) = 1.d0 / 2.d0

  x_superbrick(66) = 0.d0
  y_superbrick(66) = 0.d0
  z_superbrick(66) = 2.d0

  x_superbrick(67) = 0.d0
  y_superbrick(67) = 0.d0
  z_superbrick(67) = 0.d0

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

  ibool_superbrick(1, 8) = 24
  ibool_superbrick(2, 8) = 25
  ibool_superbrick(3, 8) = 26
  ibool_superbrick(4, 8) = 27
  ibool_superbrick(5, 8) = 14
  ibool_superbrick(6, 8) = 18
  ibool_superbrick(7, 8) = 11
  ibool_superbrick(8, 8) = 15

  ibool_superbrick(1, 9) =  1
  ibool_superbrick(2, 9) =  2
  ibool_superbrick(3, 9) = 28
  ibool_superbrick(4, 9) = 29
  ibool_superbrick(5, 9) =  5
  ibool_superbrick(6, 9) =  6
  ibool_superbrick(7, 9) = 30
  ibool_superbrick(8, 9) = 31

  ibool_superbrick(1,10) = 29
  ibool_superbrick(2,10) = 28
  ibool_superbrick(3,10) = 32
  ibool_superbrick(4,10) = 33
  ibool_superbrick(5,10) = 31
  ibool_superbrick(6,10) = 30
  ibool_superbrick(7,10) = 34
  ibool_superbrick(8,10) = 35

  ibool_superbrick(1,11) = 13
  ibool_superbrick(2,11) = 14
  ibool_superbrick(3,11) = 36
  ibool_superbrick(4,11) = 37
  ibool_superbrick(5,11) = 17
  ibool_superbrick(6,11) = 18
  ibool_superbrick(7,11) = 34
  ibool_superbrick(8,11) = 32

  ibool_superbrick(1,12) = 20
  ibool_superbrick(2,12) = 19
  ibool_superbrick(3,12) = 38
  ibool_superbrick(4,12) = 39
  ibool_superbrick(5,12) =  1
  ibool_superbrick(6,12) =  2
  ibool_superbrick(7,12) = 28
  ibool_superbrick(8,12) = 29

  ibool_superbrick(1,13) =  2
  ibool_superbrick(2,13) = 17
  ibool_superbrick(3,13) = 32
  ibool_superbrick(4,13) = 28
  ibool_superbrick(5,13) =  6
  ibool_superbrick(6,13) = 18
  ibool_superbrick(7,13) = 34
  ibool_superbrick(8,13) = 30

  ibool_superbrick(1,14) = 39
  ibool_superbrick(2,14) = 38
  ibool_superbrick(3,14) = 37
  ibool_superbrick(4,14) = 40
  ibool_superbrick(5,14) = 29
  ibool_superbrick(6,14) = 28
  ibool_superbrick(7,14) = 32
  ibool_superbrick(8,14) = 33

  ibool_superbrick(1,15) = 19
  ibool_superbrick(2,15) = 13
  ibool_superbrick(3,15) = 37
  ibool_superbrick(4,15) = 38
  ibool_superbrick(5,15) =  2
  ibool_superbrick(6,15) = 17
  ibool_superbrick(7,15) = 32
  ibool_superbrick(8,15) = 28

  ibool_superbrick(1,16) = 24
  ibool_superbrick(2,16) = 25
  ibool_superbrick(3,16) = 41
  ibool_superbrick(4,16) = 42
  ibool_superbrick(5,16) = 14
  ibool_superbrick(6,16) = 18
  ibool_superbrick(7,16) = 34
  ibool_superbrick(8,16) = 36

  ibool_superbrick(1,17) = 43
  ibool_superbrick(2,17) = 44
  ibool_superbrick(3,17) = 45
  ibool_superbrick(4,17) = 46
  ibool_superbrick(5,17) = 47
  ibool_superbrick(6,17) = 48
  ibool_superbrick(7,17) = 49
  ibool_superbrick(8,17) = 50

  ibool_superbrick(1,18) = 46
  ibool_superbrick(2,18) = 45
  ibool_superbrick(3,18) = 51
  ibool_superbrick(4,18) = 52
  ibool_superbrick(5,18) = 50
  ibool_superbrick(6,18) = 49
  ibool_superbrick(7,18) = 53
  ibool_superbrick(8,18) = 54

  ibool_superbrick(1,19) = 13
  ibool_superbrick(2,19) = 14
  ibool_superbrick(3,19) = 15
  ibool_superbrick(4,19) = 16
  ibool_superbrick(5,19) = 55
  ibool_superbrick(6,19) = 56
  ibool_superbrick(7,19) = 53
  ibool_superbrick(8,19) = 51

  ibool_superbrick(1,20) = 20
  ibool_superbrick(2,20) = 19
  ibool_superbrick(3,20) = 21
  ibool_superbrick(4,20) = 22
  ibool_superbrick(5,20) = 43
  ibool_superbrick(6,20) = 44
  ibool_superbrick(7,20) = 45
  ibool_superbrick(8,20) = 46

  ibool_superbrick(1,21) = 44
  ibool_superbrick(2,21) = 55
  ibool_superbrick(3,21) = 51
  ibool_superbrick(4,21) = 45
  ibool_superbrick(5,21) = 48
  ibool_superbrick(6,21) = 56
  ibool_superbrick(7,21) = 53
  ibool_superbrick(8,21) = 49

  ibool_superbrick(1,22) = 22
  ibool_superbrick(2,22) = 21
  ibool_superbrick(3,22) = 16
  ibool_superbrick(4,22) = 23
  ibool_superbrick(5,22) = 46
  ibool_superbrick(6,22) = 45
  ibool_superbrick(7,22) = 51
  ibool_superbrick(8,22) = 52

  ibool_superbrick(1,23) = 19
  ibool_superbrick(2,23) = 13
  ibool_superbrick(3,23) = 16
  ibool_superbrick(4,23) = 21
  ibool_superbrick(5,23) = 44
  ibool_superbrick(6,23) = 55
  ibool_superbrick(7,23) = 51
  ibool_superbrick(8,23) = 45

  ibool_superbrick(1,24) = 24
  ibool_superbrick(2,24) = 57
  ibool_superbrick(3,24) = 58
  ibool_superbrick(4,24) = 27
  ibool_superbrick(5,24) = 14
  ibool_superbrick(6,24) = 56
  ibool_superbrick(7,24) = 53
  ibool_superbrick(8,24) = 15

  ibool_superbrick(1,25) = 43
  ibool_superbrick(2,25) = 44
  ibool_superbrick(3,25) = 59
  ibool_superbrick(4,25) = 60
  ibool_superbrick(5,25) = 47
  ibool_superbrick(6,25) = 48
  ibool_superbrick(7,25) = 61
  ibool_superbrick(8,25) = 62

  ibool_superbrick(1,26) = 60
  ibool_superbrick(2,26) = 59
  ibool_superbrick(3,26) = 63
  ibool_superbrick(4,26) = 64
  ibool_superbrick(5,26) = 62
  ibool_superbrick(6,26) = 61
  ibool_superbrick(7,26) = 65
  ibool_superbrick(8,26) = 66

  ibool_superbrick(1,27) = 13
  ibool_superbrick(2,27) = 14
  ibool_superbrick(3,27) = 36
  ibool_superbrick(4,27) = 37
  ibool_superbrick(5,27) = 55
  ibool_superbrick(6,27) = 56
  ibool_superbrick(7,27) = 65
  ibool_superbrick(8,27) = 63

  ibool_superbrick(1,28) = 20
  ibool_superbrick(2,28) = 19
  ibool_superbrick(3,28) = 38
  ibool_superbrick(4,28) = 39
  ibool_superbrick(5,28) = 43
  ibool_superbrick(6,28) = 44
  ibool_superbrick(7,28) = 59
  ibool_superbrick(8,28) = 60

  ibool_superbrick(1,29) = 44
  ibool_superbrick(2,29) = 55
  ibool_superbrick(3,29) = 63
  ibool_superbrick(4,29) = 59
  ibool_superbrick(5,29) = 48
  ibool_superbrick(6,29) = 56
  ibool_superbrick(7,29) = 65
  ibool_superbrick(8,29) = 61

  ibool_superbrick(1,30) = 39
  ibool_superbrick(2,30) = 38
  ibool_superbrick(3,30) = 37
  ibool_superbrick(4,30) = 40
  ibool_superbrick(5,30) = 60
  ibool_superbrick(6,30) = 59
  ibool_superbrick(7,30) = 63
  ibool_superbrick(8,30) = 64

  ibool_superbrick(1,31) = 19
  ibool_superbrick(2,31) = 13
  ibool_superbrick(3,31) = 37
  ibool_superbrick(4,31) = 38
  ibool_superbrick(5,31) = 44
  ibool_superbrick(6,31) = 55
  ibool_superbrick(7,31) = 63
  ibool_superbrick(8,31) = 59

  ibool_superbrick(1,32) = 24
  ibool_superbrick(2,32) = 57
  ibool_superbrick(3,32) = 67
  ibool_superbrick(4,32) = 42
  ibool_superbrick(5,32) = 14
  ibool_superbrick(6,32) = 56
  ibool_superbrick(7,32) = 65
  ibool_superbrick(8,32) = 36

  end subroutine define_superbrick

