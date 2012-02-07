#!/bin/csh -f
#
# Vala Hjorleifsdottir and Qinya Liu, Caltech, Jan 2007
#
# uses `asc2sac` binary which can be compiled from UTILS/lib/asc2sac.c
# (requires SAC libraries for compilation)

#############################################################
## modify to match your location 

asc2sac_bin=/opt/seismo-util/bin/asc2sac

#############################################################

foreach file ($*)
	echo $file
    
  # call depending on which version you use:
  #
  # seismo-util version by Brian
  #set nlines = `wc -l $file | awk '{print $1}'`
  #$asc2sac_bin $file $nlines $file.sac
  
  # UTILS/lib version by Eh
  $asc2sac_bin $file
  
end
