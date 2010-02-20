#!/usr/bin/env python
#
# creates a PNG file from a polygon data file *.vtk
#
from vtk import *
import sys

# input: 
if len(sys.argv) == 2:
  modelfile = str(sys.argv[1])  
else :
  print "Usage: python plot_VTK.py OUTPUT_FILES/bin_movie_009000.d.vtk"
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
#output from: ./convert_cpt_lookuptable_python.bash blue_white_red.cpt 
colortable = vtkLookupTable()
colortable.SetNumberOfTableValues(25)
colortable.SetTableValue( 0 , 1 , 0.0392157 , 0.0392157 , 1.0 )
colortable.SetTableValue( 1 , 1 , 0.121569 , 0.121569 , 0.878431 )
colortable.SetTableValue( 2 , 1 , 0.2 , 0.2 , 0.8 )
colortable.SetTableValue( 3 , 1 , 0.278431 , 0.278431 , 0.721569 )
colortable.SetTableValue( 4 , 1 , 0.360784 , 0.360784 , 0.639216 )
colortable.SetTableValue( 5 , 1 , 0.439216 , 0.439216 , 0.560784 )
colortable.SetTableValue( 6 , 1 , 0.521569 , 0.521569 , 0.47451 )
colortable.SetTableValue( 7 , 1 , 0.6 , 0.6 , 0.396078 )
colortable.SetTableValue( 8 , 1 , 0.678431 , 0.678431 , 0.317647 )
colortable.SetTableValue( 9 , 1 , 0.760784 , 0.760784 , 0.235294 )
colortable.SetTableValue( 10 , 1 , 0.839216 , 0.839216 , 0.156863 )
colortable.SetTableValue( 11 , 1 , 0.921569 , 0.921569 , 0.0745098 )
colortable.SetTableValue( 12 , 0.996078 , 0.996078 , 0.996078 , 0.0 )
colortable.SetTableValue( 13 , 0.921569 , 0.921569 , 1 , 0.0823529 )
colortable.SetTableValue( 14 , 0.839216 , 0.839216 , 1 , 0.164706 )
colortable.SetTableValue( 15 , 0.760784 , 0.760784 , 1 , 0.243137 )
colortable.SetTableValue( 16 , 0.678431 , 0.678431 , 1 , 0.32549 )
colortable.SetTableValue( 17 , 0.6 , 0.6 , 1 , 0.403922 )
colortable.SetTableValue( 18 , 0.521569 , 0.521569 , 1 , 0.482353 )
colortable.SetTableValue( 19 , 0.439216 , 0.439216 , 1 , 0.560784 )
colortable.SetTableValue( 20 , 0.360784 , 0.360784 , 1 , 0.639216 )
colortable.SetTableValue( 21 , 0.278431 , 0.278431 , 1 , 0.721569 )
colortable.SetTableValue( 22 , 0.2 , 0.2 , 1 , 0.8 )
colortable.SetTableValue( 23 , 0.121569 , 0.121569 , 1 , 0.878431 )
colortable.SetTableValue( 24 , 0.0392157 , 0.0392157 , 1 , 1.0 )
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
PNGWriter.SetFileName("bin_color.png")
PNGWriter.Write()

#window.GetInteractor().Initialize()
#window.GetInteractor().Start()
