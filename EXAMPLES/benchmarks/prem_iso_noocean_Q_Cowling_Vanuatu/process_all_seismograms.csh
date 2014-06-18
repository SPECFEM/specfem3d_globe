#!/bin/csh

# you will find the Perl and "sac2asc" processing scripts below in directory utils/seis_process
# and you will need the Seismic Analysis Code (SAC) software http://www.iris.edu/dms/nodes/dmc/software/downloads/sac
# (i.e. the "sac" and "saclst" executable commands) to be installed and executable on your machine
# for this processing script to run

# Dimitri Komatitsch, CNRS, Marseille, France, September 2013

########### SEM spectral-element results ###########

./process_syn.pl -m CMTSOLUTION_to_filter_nonzero_hdur -a STATIONS -l 0/6000 -t 150/500 SEMD/*.sem.ascii

./rotate.pl -l 0 -L 6000 SEMD/*.MXE.sem.ascii.sac

foreach file ( SEMD/*.MXR.sem.ascii.sac SEMD/*.MXT.sem.ascii.sac SEMD/*.MXZ.sem.ascii.sac )
  echo "converting " $file " back to final result in ASCII..."
  ./sac2asc $file > $file.asciinew
end

########### QMXD normal-mode results ###########

./process_syn.pl -m CMTSOLUTION_to_filter_nonzero_hdur -a STATIONS -l 0/6000 -t 150/500 QMXD/*.qmxd

./rotate.pl -l 0 -L 6000 QMXD/*.MXE.qmxd.sac

foreach file ( QMXD/*.MXR.qmxd.sac QMXD/*.MXT.qmxd.sac QMXD/*.MXZ.qmxd.sac )
  echo "converting " $file " back to final result in ASCII..."
  ./sac2asc $file > $file.asciinew
end

