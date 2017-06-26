#! /usr/bin/env pvpython
#
#
# usage: ./paraviewpython-example.py alpha_kernel.pvsm [2]
#
# creates jpg: image*.jpg

import os
import sys
import fileinput
import string
from paraview import servermanager

number = ""

## input:
if len(sys.argv) == 2:
  filename = str(sys.argv[1])
else :
  if len(sys.argv) == 3:
    filename = str(sys.argv[1])
    number = str(sys.argv[2])
  else :
    print "usage: /paraviewpython-example.py alpha_kernel.pvsm "
    sys.exit()

outfile = "image"+number
print "file root: ",outfile

## paraview
servermanager.Connect()
view = servermanager.CreateRenderView()
servermanager.LoadState(filename)
view = servermanager.GetRenderView()

# to avoid segmentation fault
view.UseOffscreenRenderingForScreenshots = 0

## save as jpeg
jpegfilename = outfile+".jpg"
print "plot to: " + jpegfilename
view.WriteImage(jpegfilename, "vtkJPEGWriter", 1)
print

## save as png
#print "plot to: " + "image.png"
#view.WriteImage("image.png", "vtkPNGWriter",1)

## save as bmp
#print "plot to: " + "image.bmp"
#view.WriteImage("image.bmp", "vtkBMPWriter",1)

