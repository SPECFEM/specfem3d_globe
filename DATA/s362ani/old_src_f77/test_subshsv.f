c --- subshsv returns velocity and density perturbations in percent
c --- subtopo returns topography of the 410- and 650-km discontinuities in km (depression)

	implicit none

	real*4 xlat,xcolat,xlon,xdep,xrad
	real*4 vshout,vsvout,vphout,vpvout,etaout,rhoout
	real*4 topo410out,topo650out
	integer ifknowmodel 	
	
c ---
	
	write(6,"('xlat=',$)")
	read(5,*)xlat
	write(6,"('xlon=',$)")
	read(5,*)xlon
	write(6,"('xdep=',$)")
	read(5,*)xdep
	
	xcolat=90.0-xlat
	xrad=6371.0-xdep
	ifknowmodel=0	

	call subshsv(xcolat,xlon,xrad,
     #		     vshout,vsvout,vphout,vpvout,etaout,rhoout,
     #	             ifknowmodel)
	write(6,"('    vsh       vsv       vph       vpv       eta
     #       rho    ')") 
	write(6,"(6f10.5)") vshout,vsvout,vphout,vpvout,etaout,rhoout

	call subtopo(xcolat,xlon,topo410out,topo650out,ifknowmodel)
	write(6,"('   topo410    topo650 ')") 
	write(6,"(2f11.5)") topo410out,topo650out
	end
