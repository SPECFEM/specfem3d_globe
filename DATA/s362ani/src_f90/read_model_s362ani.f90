
   subroutine read_model_s362ani(THREE_D_MODEL, &
              THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
              THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)

  implicit none

  integer THREE_D_MODEL,THREE_D_MODEL_S362ANI
  integer THREE_D_MODEL_S362WMANI
  integer THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA

  integer lu
  character(len=128) modeldef
  logical exists
  integer numvar
  integer ierror

  include 'mod.h'

! -------------------------------------

  lu=1                    ! --- log unit: input 3-D model
  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI) then
    modeldef='DATA/s362ani/S362ANI'
  else if (THREE_D_MODEL == THREE_D_MODEL_S362WMANI) then
    modeldef='DATA/s362ani/S362WMANI'
  else if (THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM) then
    modeldef='DATA/s362ani/S362ANI_PREM'
  else if (THREE_D_MODEL == THREE_D_MODEL_S29EA) then
    modeldef='DATA/s362ani/S2.9EA'
  else
    stop 'unknown 3D model in read_model_s362ani'
  endif
  inquire(file=modeldef,exist=exists)
  if (exists) then
    call gt3dmodl(lu,modeldef, &
        maxhpa,maxker,maxcoe, &
        numhpa,numker,numcoe,lmxhpa, &
        ihpakern,itypehpa,coe, &
        itpspl,xlaspl,xlospl,radspl, &
        numvar,ivarkern,varstr, &
        refmdl,kerstr,hsplfl,dskker,ierror)
  else
    write(*,"('the model ',a,' does not exits')") &
             modeldef(1:len_trim(modeldef))
  endif

!         --- check arrays

  if (numker > maxker) stop 'numker > maxker'
  do ihpa=1,numhpa
    if (itypehpa(ihpa) == 1) then
      if (lmxhpa(ihpa) > maxl) stop 'lmxhpa(ihpa) > maxl'
    else if (itypehpa(ihpa) == 2) then
      if (numcoe(ihpa) > maxcoe) stop 'numcoe(ihpa) > maxcoe'
    else
      stop 'problem with itypehpa'
    endif
  enddo

  end subroutine read_model_s362ani

