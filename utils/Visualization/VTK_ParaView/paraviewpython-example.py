#! /usr/bin/env pvpython

import os
import sys
import fileinput
import string
from paraview import servermanager

## input: 
if len(sys.argv) == 3:
  filename = str(sys.argv[1])  
  icounter = str(sys.argv[2])
else:
  print "usage: ./paraviewpython-example.py state-file counter[e.g.=0]"
  sys.exit()
  
outfile = "paraview_movie."
print "file root: ",outfile

## paraview
servermanager.Connect()
view = servermanager.CreateRenderView()
servermanager.LoadState(filename)
view = servermanager.GetRenderView()
view.UseOffscreenRenderingForScreenshots = 1

## save as jpeg
jpegfilename = outfile+str(icounter)+".jpg"
print "plot to: " + jpegfilename
view.WriteImage(jpegfilename, "vtkJPEGWriter", 1)
print

