c --- evaluate perturbations in per cent 

	subroutine subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv)

	implicit none

	real*4 xcolat,xlon,xrad
	real*4 dvsh,dvsv,dvph,dvpv
 
c --- model evaluation

	integer ish	! --- 0 if SV, 1 if SH
	integer ieval	! --- 1 for velocity, 2 for anisotropy
	real*4 valu(2)	! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
	real*4 value	! --- used in single evaluation of perturbation
	integer isel	! --- if variable should be included
	real*4 depth	! --- depth
	real*4 x,y	! --- lat lon
	real*4 vsh3drel 	! --- relative perturbation
	real*4 vsv3drel 	! --- relative perturbation

c --- 
	
	integer iker,i
	character*40 vstr
	integer lstr
	integer lnblnk,ierror

	include 'mod.h'

c -------------------------------------

	depth=6371.0-xrad
	call evradker (depth,kerstr,numker,vercof,vercofd,ierror)
	if(ierror.ne.0) stop 'ierror evradker'

c --- loop over sv and sh (sv=0,sh=1)

	do ish=0,1			

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

c         --- evaluate 3-D perturbations in velocity and anisotropy 
          
	  valu(1)=0. ! --- velocity 
	  valu(2)=0. ! --- anisotropy

	  do ieval=1,2
	    value=0.
	    do iker=1,numker
	      isel=0
	      lstr=lnblnk(varstr(ivarkern(iker)))
	      vstr=(varstr(ivarkern(iker)))
	      if(ieval.eq.1) then
	        if(vstr(1:lstr).eq.'UM (SH+SV)*0.5,'.or.
     #		   vstr(1:lstr).eq.'LM (SH+SV)*0.5,'.or.
     #		   vstr(1:lstr).eq.'EA (SH+SV)*0.5,') then
	          isel=1
		endif
	      else if(ieval.eq.2) then
	        if(vstr(1:lstr).eq.'UM SH-SV,'.or.
     #	       	   vstr(1:lstr).eq.'LM SH-SV,'.or.
     #	       	   vstr(1:lstr).eq.'EA SH-SV,') then
	          isel=1
	        endif
	      endif

	      if(isel.eq.1) then
	        if(vercof(iker).ne.0.) then 
                  if(itypehpa(ihpakern(iker)).eq.1) then
		    ihpa=ihpakern(iker)
                    nylm=(lmxhpa(ihpakern(iker))+1)**2
                    do i=1,nylm
                      value=value+vercof(iker)*ylmcof(i,ihpa)
     #				*coe(i,iker)
                    enddo
                  else if(itypehpa(ihpakern(iker)).eq.2) then
		    ihpa=ihpakern(iker)
                    do i=1,nconpt(ihpa)
                      iver=iconpt(i,ihpa)
                      value=value+vercof(iker)*conpt(i,ihpa)
     #				*coe(iver,iker)
                    enddo
                  else
                    write(6,"('problem 2')")
                    stop
                  endif ! --- itypehpa
	        endif ! --- vercof(iker).ne.0.
	      endif ! --- isel.eq.1
	    enddo ! --- end of do iker=1,numker
	   
	    valu(ieval)=value
	  enddo ! --- ieval

c 	  --- evaluate perturbations in vsh and vsv 

	  if(ish.eq.1) then 
	    vsh3drel=valu(1)+0.5*valu(2) 	
	  else if(ish.eq.0) then
	    vsv3drel=valu(1)-0.5*valu(2) 	
	  else
	    stop 'something wrong'
	  endif

	enddo ! --- by ish

c --- evaluate perturbations in per cent

	dvsh=vsh3drel
	dvsv=vsv3drel
	dvph=0.55*dvsh	! --- scaling used in the inversion
	dvpv=0.55*dvsv	! --- scaling used in the inversion
	end


	

