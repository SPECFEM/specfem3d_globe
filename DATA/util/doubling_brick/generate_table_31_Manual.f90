
program generate

! generate Table 3.1 of the SPECFEM3D_GLOBE manual

implicit none

integer NPROC_XI,c,NEX_XI,compteur

integer, parameter :: NB_COLONNES = 7
integer, dimension(NB_COLONNES) :: store_NEX_XI

! base value depends if we implement three or four doublings
!!!integer, parameter :: BASE_VALUE = 16
integer, parameter :: BASE_VALUE = 32

do NPROC_XI = 1,26

  compteur = 1

  do c = 1,20
    NEX_XI = BASE_VALUE * c * NPROC_XI
    if(NEX_XI >= 64 .and. compteur <= NB_COLONNES) then
      store_NEX_XI(compteur) = NEX_XI
      compteur = compteur + 1
    endif
  enddo

  write(*,"(i6,' | ',i6,' | ',i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
     NPROC_XI, 6*NPROC_XI**2,(store_NEX_XI(compteur), compteur = 1,NB_COLONNES)

enddo

end

