#! /usr/bin/env pvpython

import os
import sys
import fileinput
import string
from paraview import servermanager


# loop
# be sure to have called : 
# > ls -1 AVS_movie_0*.inp | gawk '{print "sed -e \"s/AVS_movie_0001000.inp/"$1"/g\"< paraview_movie.pvsm > "$1".pvsm"}' | sh 
# > ls -1 AVS_movie_0*.inp.pvsm > tmp.in

f = open("tmp.in", "r")

receiverInput = open("tmp.in","r")
line = receiverInput.readline()
icounter = 0

while line :
  length= len(line)
  filename = line[0:length-1]
  print filename
  icounter = icounter + 1
  str_i = "%6.6d " % ( icounter )
  
  os.system("./paraviewpython-example.py "+filename+" "+str_i)

  #raw_input("Enter to continue")

  line = receiverInput.readline()
receiverInput.close()
