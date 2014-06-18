c --- evaluate depressions of the 410- and 650-km discontinuities in km

	subroutine subtopo(xcolat,xlon,topo410,topo650)

	implicit none

	real*4 xcolat,xlon
	real*4 topo410,topo650

c --- model evaluation

	integer ieval	! --- 1 for velocity, 2 for anisotropy
	real*4 valu(2)	! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
	real*4 value	! --- used in single evaluation of perturbation
	integer isel	! --- if variable should be included
	real*4 x,y	! --- lat lon

c --- 
	integer iker,i
	character*40 vstr
	integer lstr
	integer lnblnk

	include 'mod.h'

c -------------------------------------

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

	topo410=valu(1) 	
	topo650=valu(2) 	
	end


	

