#!/bin/csh -f

# Vala Hjorleifsdottir and Qinya Liu, Caltech, Jan 2007

foreach file ($*)
	echo $file
  set nlines = `wc -l $file | awk '{print $1}'`
  /opt/seismo-util/bin/asc2sac $file $nlines $file.sac
end
