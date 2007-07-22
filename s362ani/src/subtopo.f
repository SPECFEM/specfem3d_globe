c --- evaluate depressions of the 410- and 650-km discontinuities in km

	subroutine subtopo(xcolat,xlon,topo410out,topo650out,ifknowmodel,
     #         THREE_D_MODEL,THREE_D_MODEL_S362ANI,
     #         THREE_D_MODEL_S362WMANI,
     #         THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)

	implicit none

	real*4 xcolat,xlon
	real*4 topo410out,topo650out
	integer ifknowmodel
	integer THREE_D_MODEL,THREE_D_MODEL_S362ANI
	integer THREE_D_MODEL_S362WMANI
	integer THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA

c --- model evaluation

	integer ieval	! --- 1 for velocity, 2 for anisotropy
	real*4 valu(2)	! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
	real*4 value	! --- used in single evaluation of perturbation
	integer isel	! --- if variable should be included
	real*4 x,y	! --- lat lon

c --- 
	character*128 modeldef
	integer iker,i
	character*40 vstr
	logical exists
	integer lu,numvar,lstr
	integer lnblnk,ierror

	include 'mod.h'

c -------------------------------------

	lu=1 			! --- log unit: input 3-D model 

c --- read the model if necessary

	if(ifknowmodel.ne.0) then
	else
	  if(THREE_D_MODEL .eq. THREE_D_MODEL_S362ANI) then
            modeldef='DATA/s362ani/S362ANI'
          elseif(THREE_D_MODEL .eq. THREE_D_MODEL_S362WMANI) then
            modeldef='DATA/s362ani/S362WMANI'
          elseif(THREE_D_MODEL .eq. THREE_D_MODEL_S362ANI_PREM) then
            modeldef='DATA/s362ani/S362ANI_PREM'
          elseif(THREE_D_MODEL .eq. THREE_D_MODEL_S29EA) then
            modeldef='DATA/s362ani/S2.9EA'
          else
            stop 'unknown 3D model in subshsv'
          endif
          inquire(file=modeldef,exist=exists)
          if(exists) then
            call gt3dmodl(lu,modeldef,
     #        maxhpa,maxker,maxcoe,
     #        numhpa,numker,numcoe,lmxhpa,
     #        ihpakern,itypehpa,coe,
     #        itpspl,xlaspl,xlospl,radspl,
     #        numvar,ivarkern,varstr,
     #        refmdl,kerstr,hsplfl,dskker,ierror) 
	  else
	    write(6,"('the model ',a,' does not exits')")
     #              modeldef(1:lnblnk(modeldef))
	  endif
	  ifknowmodel=1

c 	  --- check arrays

	  if(numker.gt.maxker) stop 'numker.gt.maxker'
          do ihpa=1,numhpa
           if(itypehpa(ihpa).eq.1) then
	    if(lmxhpa(ihpa).gt.maxl) stop 'lmxhpa(ihpa).gt.maxl'
	   else if(itypehpa(ihpa).eq.2) then
 	    if(numcoe(ihpa).gt.maxcoe) stop 'numcoe(ihpa).gt.maxcoe'
	   else 
	    stop 'problem with itypehpa'
	   endif
	  enddo
	endif		

c       --- contributing horizontal basis functions at xlat,xlon
	  
	y=90.0-xcolat
	x=xlon
	do ihpa=1,numhpa
            if(itypehpa(ihpa).eq.1) then
              lmax=lmxhpa(ihpa)
              call ylm(y,x,lmax,ylmcof(1,ihpa),wk1,wk2,wk3)
            else if(itypehpa(ihpa).eq.2) then
              numcof=numcoe(ihpa)
              call splcon(y,x,numcof,xlaspl(1,ihpa),
     #            xlospl(1,ihpa),radspl(1,ihpa),
     #            nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
            else
              write(6,"('problem 1')") 
            endif
	enddo

c         --- evaluate topography (depression) in km 
          
	valu(1)=0. ! --- 410  
	valu(2)=0. ! --- 650 

	do ieval=1,2
	    value=0.
	    do iker=1,numker
	      isel=0
	      lstr=lnblnk(varstr(ivarkern(iker)))
	      vstr=(varstr(ivarkern(iker)))
	      if(ieval.eq.1) then
	        if(vstr(1:lstr).eq.'Topo 400,') then
	          isel=1
		endif
	      else if(ieval.eq.2) then
	        if(vstr(1:lstr).eq.'Topo 670,') then
	          isel=1
	        endif
	      endif

	      if(isel.eq.1) then
                  if(itypehpa(ihpakern(iker)).eq.1) then
		    ihpa=ihpakern(iker)
                    nylm=(lmxhpa(ihpakern(iker))+1)**2
                    do i=1,nylm
                      value=value+ylmcof(i,ihpa)*coe(i,iker)
                    enddo
                  else if(itypehpa(ihpakern(iker)).eq.2) then
		    ihpa=ihpakern(iker)
                    do i=1,nconpt(ihpa)
                      iver=iconpt(i,ihpa)
                      value=value+conpt(i,ihpa)*coe(iver,iker)
                    enddo
                  else
                    write(6,"('problem 2')")
                    stop
                  endif ! --- itypehpa
	      endif ! --- isel.eq.1
	    enddo ! --- end of do iker=1,numker
	   
	    valu(ieval)=value
	enddo ! --- ieval

	topo410out=valu(1) 	
	topo650out=valu(2) 	
	end


	

