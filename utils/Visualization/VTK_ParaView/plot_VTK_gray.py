#!/usr/bin/env python
#
# creates a (gray-scale) PNG file from a polygon data file *.vtk
#
from vtk import *
import sys

# input: 
if len(sys.argv) == 2:
  modelfile = str(sys.argv[1])  
else :
  print "Usage: python plot_VTK_gray.py OUTPUT_FILES/bin_movie_009000.d.vtk"
  sys.exit(1)
  
print "modelfile: ",modelfile


# creates 2D projection view
# graticule
latLevel = 3
lngLevel = 3
pname = "eqc" # equidistant cylindrical (plate caree) 
pcs = vtkGeoProjection()
pcs.SetName( pname )
pcs.SetCentralMeridian( 0. )
gcs = vtkGeoProjection()
xfm = vtkGeoTransform()
xfm.SetSourceProjection( gcs )
xfm.SetDestinationProjection( pcs )

# reads in polygon data
model = vtkPolyDataReader()
model.SetFileName(modelfile)
model.Update()

# Delaunay triangulation on data
delaunay2D = vtkDelaunay2D()
delaunay2D.SetInput( model.GetOutput() )
delaunay2D.Update()
print "delaunay : "
print "   points: ",delaunay2D.GetOutput().GetNumberOfPoints()
print "   polygons: ",delaunay2D.GetOutput().GetNumberOfPolys()

xf3 = vtkTransformFilter()
xf3.SetTransform( xfm )
xf3.SetInputConnection( delaunay2D.GetOutputPort() )


# coloring
#output from: ./convert_cpt_lookuptable_python.bash gray_pyramid_inv.cpt 
colortable = vtkLookupTable()
colortable.SetNumberOfTableValues(25)
colortable.SetTableValue(  0 , 1.0 , 1.0 , 1.0 , 1 )
colortable.SetTableValue(  1 , 0.878431 , 0.878431 , 0.878431 , 1 )
colortable.SetTableValue(  2 , 0.8 , 0.8 , 0.8 , 1 )
colortable.SetTableValue(  3 , 0.721569 , 0.721569 , 0.721569 , 1 )
colortable.SetTableValue(  4 , 0.639216 , 0.639216 , 0.639216 , 1 )
colortable.SetTableValue(  5 , 0.560784 , 0.560784 , 0.560784 , 1 )
colortable.SetTableValue(  6 , 0.47451 , 0.47451 , 0.47451 , 1 )
colortable.SetTableValue(  7 , 0.396078 , 0.396078 , 0.396078 , 1 )
colortable.SetTableValue(  8 , 0.317647 , 0.317647 , 0.317647 , 1 )
colortable.SetTableValue(  9 , 0.235294 , 0.235294 , 0.235294 , 1 )
colortable.SetTableValue(  10 , 0.156863 , 0.156863 , 0.156863 , 1 )
colortable.SetTableValue(  11 , 0.0745098 , 0.0745098 , 0.0745098 , 1 )
colortable.SetTableValue(  12 , 0.00392157 , 0.00392157 , 0.00392157 , 1 )
colortable.SetTableValue(  13 , 0.0823529 , 0.0823529 , 0.0823529 , 1 )
colortable.SetTableValue(  14 , 0.164706 , 0.164706 , 0.164706 , 1 )
colortable.SetTableValue(  15 , 0.243137 , 0.243137 , 0.243137 , 1 )
colortable.SetTableValue(  16 , 0.32549 , 0.32549 , 0.32549 , 1 )
colortable.SetTableValue(  17 , 0.403922 , 0.403922 , 0.403922 , 1 )
colortable.SetTableValue(  18 , 0.482353 , 0.482353 , 0.482353 , 1 )
colortable.SetTableValue(  19 , 0.560784 , 0.560784 , 0.560784 , 1 )
colortable.SetTableValue(  20 , 0.639216 , 0.639216 , 0.639216 , 1 )
colortable.SetTableValue(  21 , 0.721569 , 0.721569 , 0.721569 , 1 )
colortable.SetTableValue(  22 , 0.8 , 0.8 , 0.8 , 1 )
colortable.SetTableValue(  23 , 0.878431 , 0.878431 , 0.878431 , 1 )
colortable.SetTableValue(  24 , 1.0 , 1.0 , 1.0 , 1 )
colortable.SetTableRange(0.0,255.0)


# creates new actor based on PolyData mapper
mapper3 = vtkPolyDataMapper()
mapper3.SetInputConnection( xf3.GetOutputPort() )
mapper3.SetScalarRange(0.0,255.0)
mapper3.ColorArrayName= "displacement" 
mapper3.ColorAttributeType= 0 
mapper3.SetLookupTable( colortable )

actor3 = vtkActor()
actor3.SetMapper( mapper3 )


# view rendering
ren = vtkRenderer()
ren.AddActor( actor3 )
ren.SetBackground(0,0,0) # black

window = vtkRenderWindow()
window.SetMultiSamples(0)
window.AddRenderer( ren )
window.SetSize(1000, 500)
window.OffScreenRenderingOn()

ren.ResetCamera()
ren.GetActiveCamera().Zoom(2.2)
window.Render()

# writes png images
windowToImage = vtkWindowToImageFilter()
windowToImage.SetInput(window)
PNGWriter = vtkPNGWriter()
PNGWriter.SetInputConnection( windowToImage.GetOutputPort() )
PNGWriter.SetFileName("bin_mask.png")
PNGWriter.Write()

#window.GetInteractor().Initialize()
#window.GetInteractor().Start()
