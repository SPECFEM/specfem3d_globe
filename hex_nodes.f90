!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine usual_hex_nodes(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! check that the parameter file is correct
  if(NGNOD /= 27) stop 'elements should have 27 control nodes'

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=1
  iaddz(10)=0

  iaddx(11)=1
  iaddy(11)=2
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=1
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=2
  iaddy(15)=2
  iaddz(15)=1

  iaddx(16)=0
  iaddy(16)=2
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=1
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=2
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=1
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=1
  iaddz(21)=0

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=1
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=2
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=1
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=1
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=1
  iaddz(27)=1

  end subroutine usual_hex_nodes

  subroutine unusual_hex_nodes1(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=4
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=4
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=2
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=4
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=4
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=2
  iaddy(8)=4
  iaddz(8)=2

! midside nodes
  iaddx(9)=2
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=4
  iaddy(10)=2
  iaddz(10)=0

  iaddx(11)=2
  iaddy(11)=4
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=1
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=4
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=4
  iaddy(15)=4
  iaddz(15)=1

  iaddx(16)=1
  iaddy(16)=4
  iaddz(16)=1

  iaddx(17)=3
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=4
  iaddy(18)=2
  iaddz(18)=2

  iaddx(19)=3
  iaddy(19)=4
  iaddz(19)=2

  iaddx(20)=2
  iaddy(20)=2
  iaddz(20)=2

! side center nodes

  iaddx(21)=2
  iaddy(21)=2
  iaddz(21)=0

  iaddx(22)=3
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=4
  iaddy(23)=2
  iaddz(23)=1

  iaddx(24)=3
  iaddy(24)=4
  iaddz(24)=1

  iaddx(25)=1
  iaddy(25)=2
  iaddz(25)=1

  iaddx(26)=3
  iaddy(26)=2
  iaddz(26)=2

! center node

  iaddx(27)=2
  iaddy(27)=2
  iaddz(27)=1

  end subroutine unusual_hex_nodes1

  subroutine unusual_hex_nodes1p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=4
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=4
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

! midside nodes
  iaddx(9)=2
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=4
  iaddy(10)=2
  iaddz(10)=0

  iaddx(11)=2
  iaddy(11)=4
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=3
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=3
  iaddy(15)=4
  iaddz(15)=1

  iaddx(16)=0
  iaddy(16)=4
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=2
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=4
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=2
  iaddz(20)=2

! side center nodes

  iaddx(21)=2
  iaddy(21)=2
  iaddz(21)=0

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=3
  iaddy(23)=2
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=4
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=2
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=2
  iaddz(26)=2

! center node

  iaddx(27)=2
  iaddy(27)=2
  iaddz(27)=1

  end subroutine unusual_hex_nodes1p

  subroutine unusual_hex_nodes2(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=2

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=2

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=4

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=4

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=4

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=4

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=1

  iaddx(10)=2
  iaddy(10)=2
  iaddz(10)=2

  iaddx(11)=1
  iaddy(11)=4
  iaddz(11)=1

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=2

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=3

  iaddx(15)=2
  iaddy(15)=4
  iaddz(15)=3

  iaddx(16)=0
  iaddy(16)=4
  iaddz(16)=2

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=4

  iaddx(18)=2
  iaddy(18)=2
  iaddz(18)=4

  iaddx(19)=1
  iaddy(19)=4
  iaddz(19)=4

  iaddx(20)=0
  iaddy(20)=2
  iaddz(20)=4

! side center nodes

  iaddx(21)=1
  iaddy(21)=2
  iaddz(21)=1

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=3

  iaddx(23)=2
  iaddy(23)=2
  iaddz(23)=3

  iaddx(24)=1
  iaddy(24)=4
  iaddz(24)=3

  iaddx(25)=0
  iaddy(25)=2
  iaddz(25)=2

  iaddx(26)=1
  iaddy(26)=2
  iaddz(26)=4

! center node

  iaddx(27)=1
  iaddy(27)=2
  iaddz(27)=3

  end subroutine unusual_hex_nodes2

  subroutine unusual_hex_nodes2p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=-2

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=-2

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=-1

  iaddx(10)=2
  iaddy(10)=2
  iaddz(10)=-2

  iaddx(11)=1
  iaddy(11)=4
  iaddz(11)=-1

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=0

  iaddx(15)=2
  iaddy(15)=4
  iaddz(15)=0

  iaddx(16)=0
  iaddy(16)=4
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=2
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=4
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=2
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=2
  iaddz(21)=-1

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=2
  iaddz(23)=0

  iaddx(24)=1
  iaddy(24)=4
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=2
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=2
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=2
  iaddz(27)=1

  end subroutine unusual_hex_nodes2p

  subroutine unusual_hex_nodes3(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=2
  iaddz(10)=0

  iaddx(11)=1
  iaddy(11)=4
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=2
  iaddy(15)=4
  iaddz(15)=1

  iaddx(16)=0
  iaddy(16)=4
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=2
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=4
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=2
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=2
  iaddz(21)=0

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=2
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=4
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=2
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=2
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=2
  iaddz(27)=1

  end subroutine unusual_hex_nodes3

  subroutine unusual_hex_nodes4(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=2
  iaddz(10)=0

  iaddx(11)=1
  iaddy(11)=4
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=2
  iaddy(15)=3
  iaddz(15)=1

  iaddx(16)=0
  iaddy(16)=3
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=1
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=2
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=1
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=2
  iaddz(21)=0

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=1
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=3
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=1
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=1
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=1
  iaddz(27)=1

  end subroutine unusual_hex_nodes4

  subroutine unusual_hex_nodes4p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=4
  iaddz(3)=0

  iaddx(4)=0
  iaddy(4)=4
  iaddz(4)=0

  iaddx(5)=0
  iaddy(5)=2
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=2
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=4
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=4
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=2
  iaddz(10)=0

  iaddx(11)=1
  iaddy(11)=4
  iaddz(11)=0

  iaddx(12)=0
  iaddy(12)=2
  iaddz(12)=0

  iaddx(13)=0
  iaddy(13)=1
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=1
  iaddz(14)=1

  iaddx(15)=2
  iaddy(15)=4
  iaddz(15)=1

  iaddx(16)=0
  iaddy(16)=4
  iaddz(16)=1

  iaddx(17)=1
  iaddy(17)=2
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=3
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=4
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=3
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=2
  iaddz(21)=0

  iaddx(22)=1
  iaddy(22)=1
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=3
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=4
  iaddz(24)=1

  iaddx(25)=0
  iaddy(25)=3
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=3
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=3
  iaddz(27)=1

  end subroutine unusual_hex_nodes4p

  subroutine unusual_hex_nodes6(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=-2

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=-2

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=2

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=2

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=2

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=2

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=1
  iaddz(10)=-1

  iaddx(11)=1
  iaddy(11)=2
  iaddz(11)=-2

  iaddx(12)=0
  iaddy(12)=1
  iaddz(12)=-1

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=1

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=1

  iaddx(15)=2
  iaddy(15)=2
  iaddz(15)=0

  iaddx(16)=0
  iaddy(16)=2
  iaddz(16)=0

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=2

  iaddx(18)=2
  iaddy(18)=1
  iaddz(18)=2

  iaddx(19)=1
  iaddy(19)=2
  iaddz(19)=2

  iaddx(20)=0
  iaddy(20)=1
  iaddz(20)=2

! side center nodes

  iaddx(21)=1
  iaddy(21)=1
  iaddz(21)=-1

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=1

  iaddx(23)=2
  iaddy(23)=1
  iaddz(23)=1

  iaddx(24)=1
  iaddy(24)=2
  iaddz(24)=0

  iaddx(25)=0
  iaddy(25)=1
  iaddz(25)=1

  iaddx(26)=1
  iaddy(26)=1
  iaddz(26)=2

! center node

  iaddx(27)=1
  iaddy(27)=1
  iaddz(27)=1

  end subroutine unusual_hex_nodes6

  subroutine unusual_hex_nodes6p(iaddx,iaddy,iaddz)

  implicit none

  include "constants.h"

  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! define the topology of the hexahedral elements

! corner nodes
  iaddx(1)=0
  iaddy(1)=0
  iaddz(1)=0

  iaddx(2)=2
  iaddy(2)=0
  iaddz(2)=0

  iaddx(3)=2
  iaddy(3)=2
  iaddz(3)=2

  iaddx(4)=0
  iaddy(4)=2
  iaddz(4)=2

  iaddx(5)=0
  iaddy(5)=0
  iaddz(5)=4

  iaddx(6)=2
  iaddy(6)=0
  iaddz(6)=4

  iaddx(7)=2
  iaddy(7)=2
  iaddz(7)=4

  iaddx(8)=0
  iaddy(8)=2
  iaddz(8)=4

! midside nodes
  iaddx(9)=1
  iaddy(9)=0
  iaddz(9)=0

  iaddx(10)=2
  iaddy(10)=1
  iaddz(10)=1

  iaddx(11)=1
  iaddy(11)=2
  iaddz(11)=2

  iaddx(12)=0
  iaddy(12)=1
  iaddz(12)=1

  iaddx(13)=0
  iaddy(13)=0
  iaddz(13)=2

  iaddx(14)=2
  iaddy(14)=0
  iaddz(14)=2

  iaddx(15)=2
  iaddy(15)=2
  iaddz(15)=3

  iaddx(16)=0
  iaddy(16)=2
  iaddz(16)=3

  iaddx(17)=1
  iaddy(17)=0
  iaddz(17)=4

  iaddx(18)=2
  iaddy(18)=1
  iaddz(18)=4

  iaddx(19)=1
  iaddy(19)=2
  iaddz(19)=4

  iaddx(20)=0
  iaddy(20)=1
  iaddz(20)=4

! side center nodes

  iaddx(21)=1
  iaddy(21)=1
  iaddz(21)=1

  iaddx(22)=1
  iaddy(22)=0
  iaddz(22)=2

  iaddx(23)=2
  iaddy(23)=1
  iaddz(23)=3

  iaddx(24)=1
  iaddy(24)=2
  iaddz(24)=3

  iaddx(25)=0
  iaddy(25)=1
  iaddz(25)=3

  iaddx(26)=1
  iaddy(26)=1
  iaddz(26)=4

! center node

  iaddx(27)=1
  iaddy(27)=1
  iaddz(27)=3

  end subroutine unusual_hex_nodes6p
