/*
 !=====================================================================
 !
 !          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
 !          --------------------------------------------------
 ! 
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, April 2014
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include "config.h"

#ifdef HAVE_VTK

#pragma message ("\nCompiling with: HAVE_VTK enabled\n")

#include <stdio.h>
#include <unistd.h>
using namespace std;

// VTK includes
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTimerLog.h>
#include <vtkPlane.h>
#include <vtkBox.h>
#include <vtkBoxClipDataSet.h>
#include <vtkBoundingBox.h>
#include <vtkClipPolyData.h>
#include <vtkClipVolume.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkLookupTable.h>
#include <vtkLegendBoxActor.h>
#include <vtkScalarBarActor.h>
#include <vtkLegendBoxActor.h>
#include <vtkScalarBarActor.h>
#include <vtkClipDataSet.h>
#include <vtkAppendFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSphereSource.h>
#include <vtkGeometryFilter.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkQuad.h>
#include <vtkContourFilter.h>

/* ----------------------------------------------------------------------------------------------- */

// USER PARAMETERS
int SHOW_FREESURFACE = 0;
int SHOW_VOLUMEDATA = 0;

// color range
double gcolor_min  = 0.0;
double gcolor_max  = 0.005;
double gcolor_incr = 0.001;

/* ----------------------------------------------------------------------------------------------- */


class Visualization {
  public:
    // 2D surface data
    vtkPoints* points2D;
    vtkDoubleArray* data_array2D;
    vtkCellArray* cells2D;
  
    vtkUnstructuredGrid* volume2D;

    vtkDataSetMapper* mapMesh2D;
    vtkActor* actor2D;
  
    // 3D volume data
    vtkPoints* points3D;
    vtkDoubleArray* data_array3D;
    vtkCellArray* cells;

    vtkUnstructuredGrid* volume3D;

    vtkDataSetMapper* mapMesh3D;
    vtkActor* actor3D;

    // clipping
    vtkPlane* clipPlane1;
    vtkPlane* clipPlane2;

    vtkTableBasedClipDataSet* clip1;
    vtkTableBasedClipDataSet* clip2;
    vtkTableBasedClipDataSet* clip1m;
  
    // visualizing
    vtkAppendFilter *merger;
    vtkDataSetSurfaceFilter* clippersurface;
  
    // colors
    vtkLookupTable* lut2D;
    vtkLookupTable* lut;
    int icolor;
    vtkScalarBarActor* legendcolor;

    // window text
    vtkTextActor* text;
    vtkTextActor* textGPU;
    vtkTextActor* help;

    // user interaction
    int helpOn;
    int jpegImageOn;
    int colorAdjustOn;
    int bgBlackOn;

    // rendering
    vtkRenderer* ren;
    vtkRenderWindow* renWin;
    vtkRenderWindowInteractor* iren;

    // camera
    vtkCamera* camera;
    double pcam[3]; // camera position
    double rclip[2]; // clipping range
    double pfocal[3]; // focal point
  
    // source sphere
    vtkSphereSource* sphere;
    vtkPolyDataMapper* mapperSphere;
    vtkActor* actorSphere;  
    double pos_source[3]; // source location

    // countours
    vtkContourFilter* contour;
    vtkPolyDataMapper* contourMapper;
    vtkActor* contourActor;
};

class VTKmesh{
  public:
    int np;
    int nspec;  
};

class VTKState {
  public:
    // VTK rendering
    Visualization vtk;
    
    // meshes
    VTKmesh freesurface;
    VTKmesh volume;

    // timing
    vtkTimerLog* timer;  
};

// global VTK state variable
static VTKState fs;

/* ----------------------------------------------------------------------------------------------- */

void interactor_usage() {
  cout << endl;
  cout << "interactor usage, press: " << endl;
  cout << "  a            - adjust color scale value (toggle on/off) " << endl;
  cout << "  b            - change background color (toggle white/black) " << endl;
  cout << "  c            - change color scale " << endl;
  cout << "  h            - help (toggle on/off) " << endl;
  cout << "  mouse click  - move/rotate view " << endl;
  cout << "  r            - reset view " << endl;
  cout << "  s            - save snapshot as vtu file " << endl;
  cout << "  i            - save snapshot as jpg image (toggle on/off) " << endl;
  cout << "  <space>,q    - continue simulation " << endl;
  cout << endl;
}

void save_snapshot_vtu() {
  // Write vtu file
  printf("writing unstructured grid data...\n");

  std::string filename = "test_snapshot.vtu";
  
  // creates writer
  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(fs.vtk.volume3D);
  writer->SetDataModeToAscii();
  writer->Write();
  
  //clean up
  writer->Delete();

  printf("snapshot written to file: %s\n\n",filename.c_str());
}


void save_snapshot_jpg(int it) {
  // Write vtu file
  printf("writing jpg image...\n");
  
  //std::string filename = "test_snapshot.jpg";
  
  char filename[180];
  if (it > 0 ) {
    sprintf(filename,"test_snapshot.%6.6d.jpg",it);
  }else{
    sprintf(filename,"test_snapshot.jpg");
  }
  
  // window filter
  vtkWindowToImageFilter* w2i = vtkWindowToImageFilter::New();
  w2i->SetInput(fs.vtk.renWin);
  w2i->Update();
  
  // creates writer
  vtkJPEGWriter* writer = vtkJPEGWriter::New();
  //writer->SetFileName(filename.c_str());
  writer->SetFileName(filename);
  writer->SetInputConnection(w2i->GetOutputPort());
  writer->Write();
  
  //clean up
  writer->Delete();
  w2i->Delete();
  
  printf("snapshot written to file: %s\n\n",filename);
}


void set_color_scale(int icolor) {
  if (icolor == 0 ) {
    // sets (default) rainbow color scale
    //# Set the hue range: from low to high the
    //# table passes through blue, green, yellow,
    //# orange, and red
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();
    
    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.66667, 0.0 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }else if (icolor == 1 ) {
    // red-blue color scale
    // Since the default lut is
    // a rainbow lut, we only have
    // to worry about the hue range.
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();

    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.0, 0.667 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }else if (icolor == 2 ) {
    // red color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();
    
    fs.vtk.lut->SetValueRange( 0.4, 1.0 );
    fs.vtk.lut->SetHueRange( 0.0, 0.4 );
    fs.vtk.lut->SetSaturationRange( 0.5, 0.0 );
  }else if (icolor == 3 ) {
    // blue color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();
    
    fs.vtk.lut->SetValueRange( 0.4, 1.0 );
    fs.vtk.lut->SetHueRange( 0.6, 0.6 );
    fs.vtk.lut->SetSaturationRange( 0.5, 0.0 );
  }else if (icolor == 4 ) {
    // black color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();
    
    fs.vtk.lut->SetValueRange( 0.4, 1.0 ); // from black to white
    fs.vtk.lut->SetHueRange( 0.0, 1.0 );
    fs.vtk.lut->SetSaturationRange( 0.0, 0.0 ); // no color saturation
  }else{
    // reset
    fs.vtk.icolor = 0;
    // sets (default) rainbow color scale
    //# Set the hue range: from low to high the
    //# table passes through blue, green, yellow,
    //# orange, and red
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();
    
    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.66667, 0.0 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }
  fs.vtk.lut->Build();
}

// A class not derived from vtkObjectBase
class MyInteractor{
  public:
    void KeypressCallbackFunction(vtkObject* caller,long unsigned int eventId,void* callData) {
      //cout << "Pressed: " << fs.vtk.iren->GetKeySym() << endl;
      // usage
      if (string(fs.vtk.iren->GetKeySym()) == "h") {
        interactor_usage();
        // displays help text
        if (fs.vtk.helpOn == 0 ) {
          fs.vtk.help = vtkTextActor::New();
          fs.vtk.help->GetTextProperty()->SetFontSize ( 16 );
          fs.vtk.help->SetPosition( 10, 80 );
          fs.vtk.help->GetTextProperty()->SetColor( 0.1, 0.1, 0.1 );
          // updates displayed text
          fs.vtk.help->SetInput( "interactor usage, press:\n  a            - adjust color scale value (toggle on/off)\n  b            - change background color (toggle white/black) \n  c            - change color scale \n  h            - help (toggle on/off) \n  mouse click  - move/rotate view \n  r            - reset view \n  s            - save snapshot as vtu file\n  i            - save snapshot as jpg image (toggle on/off)\n  <space>,q    - continue simulation \n" );
                    fs.vtk.ren->AddActor2D( fs.vtk.help );
          fs.vtk.renWin->Render();
          fs.vtk.helpOn = 1;
        }else{
          fs.vtk.ren->RemoveActor2D( fs.vtk.help );
          fs.vtk.renWin->Render();
          fs.vtk.help->Delete();
          fs.vtk.helpOn = 0;
        }                
      }
      // changes color scales
      if (string(fs.vtk.iren->GetKeySym()) == "c") {
        fs.vtk.icolor++;
        set_color_scale(fs.vtk.icolor);
        fs.vtk.renWin->Render();
      }
      // adjust color scaling
      if (string(fs.vtk.iren->GetKeySym()) == "a") {
        fs.vtk.colorAdjustOn++;
        fs.vtk.colorAdjustOn = fs.vtk.colorAdjustOn % 2 ;
        if (fs.vtk.colorAdjustOn) {
          printf("\ntoggle on: color adjusting\n");
        }else{
          printf("\ntoggle off: color adjusting off\n");
        }
        set_color_scale(fs.vtk.icolor);
        fs.vtk.renWin->Render();
      }      
      // saves vtu snapshot
      if (string(fs.vtk.iren->GetKeySym()) == "s") {
        save_snapshot_vtu();
      }
      // saves jpg snapshot
      if (string(fs.vtk.iren->GetKeySym()) == "i") {
        // toggles jpeg output flag
        if (fs.vtk.jpegImageOn == 0) {
          printf("\ntoggle on: save snapshot as jpeg-image\n");
          fs.vtk.jpegImageOn = 1;
        }else{
          printf("\ntoggle off: save snapshot as jpeg-image\n");
          fs.vtk.jpegImageOn = 0;
        }
      }
      // stops interaction, continues running simulation
      if (string(fs.vtk.iren->GetKeySym()) == "space" || string(fs.vtk.iren->GetKeySym()) == "q") {
        // updates text
        fs.vtk.text->SetInput( "...continue simulation " );
        fs.vtk.renWin->Render();        
        
        // stops interactive renderer
        fs.vtk.iren->TerminateApp();
      }
      // stops interaction, continues running simulation
      if (string(fs.vtk.iren->GetKeySym()) == "r") {
        // reposition the camera, so that actor can be fully seen
        //fs.vtk.ren->ResetCamera(); // only resets zoom effect

        // range
        fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
        // focal point
        fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
        // position
        fs.vtk.camera->SetPosition( fs.vtk.pcam );

        // view
        fs.vtk.camera->SetViewAngle( 30.0 );
        fs.vtk.camera->SetViewUp( 0.0, 0.0, 1.0 );
        
        fs.vtk.renWin->Render();        
      }
      // changes color scale maximum values
      if (string(fs.vtk.iren->GetKeySym()) == "Up") {
        // increases color max by 10%
        gcolor_max = gcolor_max + gcolor_incr;
        fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
        fs.vtk.lut->Build();
        fs.vtk.renWin->Render();
      }
      if (string(fs.vtk.iren->GetKeySym()) == "Down") {
        // decreases color max by 10%
        gcolor_max = gcolor_max - gcolor_incr;
        if (gcolor_max < gcolor_min ) gcolor_max = gcolor_min;
        fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
        fs.vtk.lut->Build();
        fs.vtk.renWin->Render();
      }
      // changes background
      if (string(fs.vtk.iren->GetKeySym()) == "b") {
        // toggles background
        if (fs.vtk.bgBlackOn == 0) {
          printf("\ntoggle on: background black/white\n");
          fs.vtk.bgBlackOn = 1;
          fs.vtk.ren->SetBackground(0,0,0); // Background color black
          fs.vtk.renWin->Render();
        }else{
          printf("\ntoggle off: background black/white\n");
          fs.vtk.bgBlackOn = 0;
          fs.vtk.ren->SetBackground(1,1,1); // Background color white
          fs.vtk.renWin->Render();
        }
      }
      // toggles freesurface visibility
      if (string(fs.vtk.iren->GetKeySym()) == "1") {
        if (SHOW_FREESURFACE == 1) {
          if (fs.vtk.actor2D->GetVisibility() == 1) {
            fs.vtk.contourActor->SetVisibility( 0 );
            fs.vtk.actor2D->SetVisibility( 0 );
          }else{
            fs.vtk.contourActor->SetVisibility( 1 );
            fs.vtk.actor2D->SetVisibility( 1 );
          }
          fs.vtk.renWin->Render();
        }
      }      
    }
};


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int* GPU_MODE) {

  // user output
  printf(" VTK window initialization\n");

  // initializes window
  //
  // camera position
  fs.vtk.camera = vtkCamera::New();

  // initial settings
  fs.vtk.rclip[0] = -1.5;
  fs.vtk.rclip[1] = 1.5;
  fs.vtk.pfocal[0] = 0.0;
  fs.vtk.pfocal[1] = 0.0;
  fs.vtk.pfocal[2] = 0.0;
  fs.vtk.pcam[0] = 1.0;
  fs.vtk.pcam[1] = 1.0;
  fs.vtk.pcam[2] = 0.3;

  // range
  fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
  // camer focal point
  fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
  // camer position
  fs.vtk.camera->SetPosition( fs.vtk.pcam );

  // view
  fs.vtk.camera->SetViewAngle( 30.0 );
  fs.vtk.camera->SetViewUp( 0.0, 0.0, 1.0 );

  // renderer
  fs.vtk.ren = vtkRenderer::New();
  fs.vtk.ren->SetActiveCamera(fs.vtk.camera);
  
  // Background color white
  fs.vtk.bgBlackOn = 0;
  fs.vtk.ren->SetBackground(1,1,1);
  
  // text actors
  // GPU flag
  fs.vtk.textGPU = vtkTextActor::New();
  fs.vtk.textGPU->GetTextProperty()->SetFontSize ( 32 );
  fs.vtk.textGPU->SetPosition( 10, 560 );
  if (*GPU_MODE) {
    fs.vtk.textGPU->SetInput( "GPU" );
    fs.vtk.textGPU->GetTextProperty()->SetColor( 0.2,0.8,0.2 );
  }else{
    fs.vtk.textGPU->SetInput( "CPU" );
    fs.vtk.textGPU->GetTextProperty()->SetColor( 0.8,0.2,0.2 );
  }
  fs.vtk.ren->AddActor2D( fs.vtk.textGPU );
  
  // progress text
  fs.vtk.text = vtkTextActor::New();
  fs.vtk.text->GetTextProperty()->SetFontSize ( 16 );
  fs.vtk.text->SetPosition( 10, 530 );
  fs.vtk.text->SetInput( "...initializing data " );
  fs.vtk.text->GetTextProperty()->SetColor( 0.5,0.5,0.5 );
  fs.vtk.ren->AddActor2D( fs.vtk.text );

  // color table
  int tableSize = 256;
  fs.vtk.icolor = 0; // from blue to red
  fs.vtk.colorAdjustOn = 1; // automatic adjust
  
  fs.vtk.lut = vtkLookupTable::New();
  fs.vtk.lut->SetNumberOfColors( tableSize );

  // sets (default) rainbow color scale
  set_color_scale(fs.vtk.icolor);

  // render window
  fs.vtk.renWin = vtkRenderWindow::New();
  fs.vtk.renWin->AddRenderer( fs.vtk.ren );

  fs.vtk.renWin->SetPosition(500,0);
  fs.vtk.renWin->SetSize(900,600);
  
  // window interactor
  fs.vtk.iren = vtkRenderWindowInteractor::New();
  fs.vtk.iren->SetRenderWindow(fs.vtk.renWin);
  fs.vtk.iren->Initialize();

  // callback function
  MyInteractor mykey;
  fs.vtk.iren->AddObserver(vtkCommand::KeyPressEvent, &mykey, &MyInteractor::KeypressCallbackFunction);

  // renders window
  fs.vtk.renWin->Render();

  // timer
  fs.timer = vtkTimerLog::New();
  
  //printf("      initialized VTK window successfully\n");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float* xs_x,float* xs_y, float* xs_z) {

  double xyz[3];
  
  // sets source location
  fs.vtk.pos_source[0] = *xs_x;
  fs.vtk.pos_source[1] = *xs_y;
  fs.vtk.pos_source[2] = *xs_z;
  printf("     sphere location: %f %f %f \n",fs.vtk.pos_source[0],fs.vtk.pos_source[1],fs.vtk.pos_source[2]);
  
  // creates sphere around source location
  fs.vtk.sphere = vtkSphereSource::New();
  fs.vtk.sphere->SetCenter(fs.vtk.pos_source);
  fs.vtk.sphere->SetRadius(0.02);

  fs.vtk.mapperSphere = vtkPolyDataMapper::New();
  fs.vtk.mapperSphere->SetInputConnection(fs.vtk.sphere->GetOutputPort());
 
  fs.vtk.actorSphere = vtkActor::New();
  fs.vtk.actorSphere->SetMapper(fs.vtk.mapperSphere);

  // adds actor
  fs.vtk.ren->AddActor(fs.vtk.actorSphere);

  // render window
  fs.vtk.renWin->Render();

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int* free_np,
                                                             float* free_x,
                                                             float* free_y,
                                                             float* free_z,
                                                             int* free_nspec,
                                                             int* free_conn) {
  double xyz[3];
  int id1,id2,id3,id4;
  
  // user output
  //printf("    prepare VTK freesurface...\n");

  // initializes
  SHOW_FREESURFACE = 1;
  
  // number of points
  fs.freesurface.np = *free_np;
  printf("     free surface points: %d\n", fs.freesurface.np);

  // checks
  if (fs.freesurface.np == 0) {
    fprintf(stderr,"Error: VTK_MODE without 2D freesurface points \n");
    exit(1);
  }

  // creates new points and data arrays
  fs.vtk.points2D = vtkPoints::New();
  fs.vtk.points2D->SetNumberOfPoints(fs.freesurface.np);

  fs.vtk.data_array2D = vtkDoubleArray::New();
  fs.vtk.data_array2D->SetNumberOfComponents(1);
  fs.vtk.data_array2D->SetNumberOfValues(fs.freesurface.np);
  fs.vtk.data_array2D->SetName("topo");

  for(int i = 0;i<fs.freesurface.np;i++) {
    xyz[0] = free_x[i];
    xyz[1] = free_y[i];
    xyz[2] = free_z[i];
    fs.vtk.points2D->SetPoint(i,xyz);
    fs.vtk.data_array2D->SetValue(i,6371000.0*(sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2])-1.0) );
  }
  // Find min and max
  double bounds[2];
  fs.vtk.data_array2D->GetValueRange(bounds);
  double min = bounds[0];
  double max = bounds[1];
  printf("     topo: min = %f max = %f\n", min,max);

  // black color scale
  int tableSize = 256;  
  fs.vtk.lut2D = vtkLookupTable::New();
  fs.vtk.lut2D->SetNumberOfColors(tableSize);
  fs.vtk.lut2D->SetValueRange( 0.4, 1.0 ); // from black to white
  fs.vtk.lut2D->SetHueRange( 0.0, 1.0 );
  fs.vtk.lut2D->SetSaturationRange( 0.0, 0.0 ); // no color saturation
  fs.vtk.lut2D->SetTableRange( min, max );
  fs.vtk.lut2D->Build();
  
  // creates cell connectivity
  fs.freesurface.nspec = *free_nspec;

  // cells
  fs.vtk.cells2D = vtkCellArray::New();
  vtkQuad* quad = vtkQuad::New();
  for(int ispec = 0;ispec<fs.freesurface.nspec;ispec++) {
    id1 = free_conn[0+ispec*4];
    id2 = free_conn[1+ispec*4];
    id3 = free_conn[2+ispec*4];
    id4 = free_conn[3+ispec*4];
    quad->GetPointIds()->SetId(0,id1);
    quad->GetPointIds()->SetId(1,id2);
    quad->GetPointIds()->SetId(2,id3);
    quad->GetPointIds()->SetId(3,id4);
    fs.vtk.cells2D->InsertNextCell(quad);
  }
  quad->Delete();

  fs.vtk.volume2D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume2D->SetPoints(fs.vtk.points2D);
  fs.vtk.volume2D->GetPointData()->SetScalars(fs.vtk.data_array2D);
  // cells
  fs.vtk.volume2D->SetCells(VTK_QUAD, fs.vtk.cells2D);

  // contour iso-surfacing
  fs.vtk.contour = vtkContourFilter::New();
  fs.vtk.contour->SetInput(  fs.vtk.volume2D );
  fs.vtk.contour->SetNumberOfContours(7);
  fs.vtk.contour->SetValue(0, 0.0);
  fs.vtk.contour->SetValue(1, min*0.25);
  fs.vtk.contour->SetValue(2, max*0.25);
  fs.vtk.contour->SetValue(3, min*0.5);
  fs.vtk.contour->SetValue(4, max*0.5);
  fs.vtk.contour->SetValue(5, min*0.75);
  fs.vtk.contour->SetValue(6, max*0.75);

  fs.vtk.contour->ComputeScalarsOff();
  fs.vtk.contour->ComputeGradientsOff();

  fs.vtk.contourMapper = vtkPolyDataMapper::New();
  fs.vtk.contourMapper->SetInputConnection(0, fs.vtk.contour->GetOutputPort() );
  fs.vtk.contourMapper->ScalarVisibilityOff();
  
  fs.vtk.contourActor = vtkActor::New();
  fs.vtk.contourActor->SetMapper( fs.vtk.contourMapper );
  fs.vtk.contourActor->GetProperty()->SetOpacity( 1.0 );
  fs.vtk.contourActor->SetScale( 1.001,1.001,1.001 );
  fs.vtk.contourActor->GetProperty()->SetColor( 0.5, 0.5, 0.5 );
  
  fs.vtk.ren->AddActor(fs.vtk.contourActor);
  
  // mapping
  fs.vtk.mapMesh2D = vtkDataSetMapper::New();
  fs.vtk.mapMesh2D->SetInput( fs.vtk.volume2D );
  
  fs.vtk.mapMesh2D->SetLookupTable(fs.vtk.lut2D);
  fs.vtk.mapMesh2D->SetColorModeToMapScalars();
  fs.vtk.mapMesh2D->SetScalarModeToUsePointData();
  fs.vtk.mapMesh2D->ScalarVisibilityOn();
  fs.vtk.mapMesh2D->UseLookupTableScalarRangeOn();

  // actor
  fs.vtk.actor2D = vtkActor::New();
  fs.vtk.actor2D->SetMapper(fs.vtk.mapMesh2D);
  fs.vtk.actor2D->GetProperty()->SetOpacity( 0.5 );
  fs.vtk.actor2D->SetScale( 1.001,1.001,1.001 );

  // 2D actor
  fs.vtk.ren->AddActor(fs.vtk.actor2D);

  // render window
  fs.vtk.renWin->Render();
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int* vol_np,
                                                 float* vol_x,
                                                 float* vol_y,
                                                 float* vol_z,
                                                 int* vol_nspec,
                                                 int* vol_conn) {

  double xyz[3];
  double data_bounds[2];
  double model_bounds[6];  
  int id1,id2,id3,id4,id5,id6,id7,id8;
  static int debug_file = 0;
  
  // user output
  //printf("    prepare VTK wavefield ...\n");

  // initializes
  SHOW_VOLUMEDATA = 1;
  
  // volumetric wavefield
  // text
  fs.vtk.text->SetInput( "...adding wavefield " );
  
  // update view
  fs.vtk.renWin->Render();
  
  // volume mesh
  fs.volume.np = *vol_np;
  printf("     volume points: %d\n", fs.volume.np);

  // checks
  if (fs.volume.np == 0) {
    fprintf(stderr,"Error: VTK_MODE without 3D volume points \n");
    exit(1);
  }
      
  // volume data
  // point locations
  fs.vtk.points3D = vtkPoints::New();
  fs.vtk.points3D->SetNumberOfPoints(fs.volume.np);

  // point data (wavefield)
  fs.vtk.data_array3D = vtkDoubleArray::New();
  fs.vtk.data_array3D->SetNumberOfComponents(1);
  fs.vtk.data_array3D->SetName("vnorm");
  fs.vtk.data_array3D->SetNumberOfValues(fs.volume.np);    

  for(int i = 0;i<fs.volume.np;i++) {
    xyz[0] = vol_x[i];
    xyz[1] = vol_y[i];
    xyz[2] = vol_z[i];
    fs.vtk.points3D->SetPoint(i,xyz);
    fs.vtk.data_array3D->SetValue(i,0.0);
    //fs.vtk.data_array3D->SetValue(i,i*1.0);
  }
  //fs.vtk.data_array3D->Update();
  // Find min and max
  fs.vtk.data_array3D->GetValueRange(data_bounds);
  double min = data_bounds[0];
  double max = data_bounds[1];
  printf("     data: min = %f max = %f\n",min,max);

  // Find min and max
  fs.vtk.points3D->GetBounds(model_bounds);
  double xmin = model_bounds[0];
  double xmax = model_bounds[1];
  double ymin = model_bounds[2];
  double ymax = model_bounds[3];
  double zmin = model_bounds[4];
  double zmax = model_bounds[5];
  printf("     model: xmin/xmax = %f / %f\n",xmin,xmax);
  printf("     model: ymin/ymax = %f / %f\n",ymin,ymax);
  printf("     model: zmin/zmax = %f / %f\n",zmin,zmax);
  
  // adjusts camera settings
  fs.vtk.rclip[0] = xmin - 0.1*(xmax-xmin);
  fs.vtk.rclip[1] = xmax + 0.1*(xmax-xmin);
  fs.vtk.pfocal[0] = fs.vtk.pos_source[0]; //(xmax-xmin)/2.0;
  fs.vtk.pfocal[1] = fs.vtk.pos_source[1]; //(ymax-ymin)/2.0;
  fs.vtk.pfocal[2] = fs.vtk.pos_source[2]; // (zmax-zmin)/2.0;
  fs.vtk.pcam[0] = fs.vtk.pos_source[0] + 0.5*(zmax-zmin); // (xmax-xmin)/2.0 - 0.1*(xmax-xmin);
  fs.vtk.pcam[1] = fs.vtk.pos_source[1] + 0.5*(zmax-zmin); // (ymax-ymin)/2.0 + 0.1*(ymax-ymin);
  fs.vtk.pcam[2] = fs.vtk.pos_source[2] + 0.5*(zmax-zmin); // (zmax-zmin)/2.0 + 0.5*(zmax-zmin);
  
  // range
  fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
  // camer focal point
  fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
  // camer position
  fs.vtk.camera->SetPosition( fs.vtk.pcam );
  
  //
  // unstructured grid
  //
  // creates cell connectivity
  fs.volume.nspec = *vol_nspec;

  // cells
  fs.vtk.cells = vtkCellArray::New();
  vtkHexahedron* hex = vtkHexahedron::New();
  for(int ispec = 0;ispec<fs.volume.nspec;ispec++) {
    id1 = vol_conn[0+ispec*8];
    id2 = vol_conn[1+ispec*8];
    id3 = vol_conn[2+ispec*8];
    id4 = vol_conn[3+ispec*8];
    id5 = vol_conn[4+ispec*8];
    id6 = vol_conn[5+ispec*8];
    id7 = vol_conn[6+ispec*8];
    id8 = vol_conn[7+ispec*8];
    hex->GetPointIds()->SetId(0,id1);
    hex->GetPointIds()->SetId(1,id2);
    hex->GetPointIds()->SetId(2,id3);
    hex->GetPointIds()->SetId(3,id4);
    hex->GetPointIds()->SetId(4,id5);
    hex->GetPointIds()->SetId(5,id6);
    hex->GetPointIds()->SetId(6,id7);
    hex->GetPointIds()->SetId(7,id8);
    fs.vtk.cells->InsertNextCell(hex);      
  }
  hex->Delete();
  
  fs.vtk.volume3D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume3D->SetPoints(fs.vtk.points3D);
  fs.vtk.volume3D->GetPointData()->SetScalars(fs.vtk.data_array3D);
  // cells
  fs.vtk.volume3D->SetCells(VTK_HEXAHEDRON, fs.vtk.cells);

  //
  // clip box
  //
  // bounds
  double bb_all[6] = { xmin,xmax,ymin,ymax,zmin,zmax };

  // normalized source location
  double norm  = sqrt( fs.vtk.pos_source[0]*fs.vtk.pos_source[0]
                   + fs.vtk.pos_source[1]*fs.vtk.pos_source[1]
                   + fs.vtk.pos_source[2]*fs.vtk.pos_source[2]);
  // normalized vector to source
  double v[3];
  if (norm > 0.0) {
    v[0] = fs.vtk.pos_source[0]/norm;
    v[1] = fs.vtk.pos_source[1]/norm;
    v[2] = fs.vtk.pos_source[2]/norm;
  }else{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 1.0;
  }
  
  // source location vector / perpendicular vectors
  //
  // vectorproduct(vector1, vector2, product)
  // calculates the vector product of vector1 and vector2
  //    product(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
  //    product(2) = - vector1(1)*vector2(3) + vector1(3)*vector2(1)
  //    product(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)
  double p1[3],p2[3];
  if (norm > 0) {
    p1[0] = 0.0;
    p1[1] = 0.0;
    p1[2] = 1.0;
    
    p2[0] = 0.0;
    p2[1] = 1.0;
    p2[2] = 0.0;
  }else{
    p1[0] = 0.0;
    p1[1] = 1.0;
    p1[2] = 0.0;
    
    p2[0] = 1.0;
    p2[1] = 0.0;
    p2[2] = 0.0;
  }
  double pn1[3],pn2[3];
  pn1[0] = v[1]*p1[2] - v[2]*p1[1];
  pn1[1] = -v[0]*p1[2] + v[2]*p1[0];
  pn1[2] = v[0]*p1[1] - v[1]*p1[0];

  pn2[0] = v[1]*p2[2] - v[2]*p2[1];
  pn2[1] = -v[0]*p2[2] + v[2]*p2[0];
  pn2[2] = v[0]*p2[1] - v[1]*p2[0];
  
  // flips normal depending on location of source
  if (fabs(v[0]-xmin) < fabs(v[0]-xmax) ) {
    pn1[0] *= -1.0;
    pn1[1] *= -1.0;
    pn1[2] *= -1.0;
  }
  if (fabs(v[1]-ymin) < fabs(v[1]-ymax) ) {
    pn2[0] *= -1.0;
    pn2[1] *= -1.0;
    pn2[2] *= -1.0;
  }
  
  // plane through source location
  fs.vtk.clipPlane1 = vtkPlane::New();
  fs.vtk.clipPlane1->SetNormal( pn1 );
  fs.vtk.clipPlane1->SetOrigin(0.0, 0.0, 0.0);

  fs.vtk.clipPlane2 = vtkPlane::New();
  fs.vtk.clipPlane2->SetNormal( pn2 );
  fs.vtk.clipPlane2->SetOrigin(0.0, 0.0, 0.0);

  // table based clipping
  fs.vtk.clip1 = vtkTableBasedClipDataSet::New();
  fs.vtk.clip1->SetInput( fs.vtk.volume3D );
  fs.vtk.clip1->SetClipFunction( fs.vtk.clipPlane1 );
  fs.vtk.clip1->Update();

  fs.vtk.clip1m = vtkTableBasedClipDataSet::New();
  fs.vtk.clip1m->SetInput( fs.vtk.volume3D );
  fs.vtk.clip1m->SetClipFunction( fs.vtk.clipPlane1 );
  fs.vtk.clip1m->InsideOutOn();
  fs.vtk.clip1m->Update();

  fs.vtk.clip2 = vtkTableBasedClipDataSet::New();
  fs.vtk.clip2->SetInput( fs.vtk.clip1m->GetOutput() );
  fs.vtk.clip2->SetClipFunction( fs.vtk.clipPlane2 );
  fs.vtk.clip2->Update();

  // merges mesh parts
  fs.vtk.merger = vtkAppendFilter::New();
  fs.vtk.merger->MergePointsOn();
  fs.vtk.merger->AddInput( fs.vtk.clip1->GetOutput() );
  fs.vtk.merger->AddInput( fs.vtk.clip2->GetOutput() );
  fs.vtk.merger->Update();
  
  // test file
  if (debug_file == 1) {
    vtkUnstructuredGrid* data = fs.vtk.merger->GetOutput();
    cout << "    merger: cells  = " << data->GetNumberOfCells() << endl;
    cout << "    merger: points = " << data->GetNumberOfPoints() << endl;
    // Write vtu file
    printf("\nwriting unstructured grid data...\n");
    std::string filename = "test_init_snapshot.vtu";
    // creates writer
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    writer->SetDataModeToAscii();
    writer->Write();
    //clean up
    writer->Delete();
    printf("snapshot written to file: %s\n\n",filename.c_str());
  }
  
  // clipper surface
  fs.vtk.clippersurface = vtkDataSetSurfaceFilter::New();
  fs.vtk.clippersurface->SetInputConnection(0, fs.vtk.merger->GetOutputPort());
  fs.vtk.clippersurface->Update();
  
  // cell connectivity mapping
  fs.vtk.mapMesh3D = vtkDataSetMapper::New();
  fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.clippersurface->GetOutputPort());

  fs.vtk.mapMesh3D->SetLookupTable(fs.vtk.lut);
  fs.vtk.mapMesh3D->SetColorModeToMapScalars();
  fs.vtk.mapMesh3D->SetScalarModeToUsePointData();
  fs.vtk.mapMesh3D->ScalarVisibilityOn();
  fs.vtk.mapMesh3D->UseLookupTableScalarRangeOn();
  
  //actor
  fs.vtk.actor3D = vtkActor::New();
  fs.vtk.actor3D->SetMapper(fs.vtk.mapMesh3D);
  //fs.vtk.actor3D->GetProperty()->SetRepresentationToSurface();
  //fs.vtk.actor3D->GetProperty()->SetEdgeVisibility(1);
  
  // 3D actor
  fs.vtk.ren->AddActor(fs.vtk.actor3D);
  
  // legend for colors
  fs.vtk.legendcolor = vtkScalarBarActor::New();
  fs.vtk.legendcolor->SetLookupTable(fs.vtk.lut);
  fs.vtk.legendcolor->SetTitle( "vnorm" );
  fs.vtk.legendcolor->SetWidth( 0.1 );
  //fs.vtk.legendcolor->SetFontSize ( 16 );
  fs.vtk.ren->AddActor(fs.vtk.legendcolor);

  // reposition the camera, so that actor can be fully seen
  fs.vtk.ren->ResetCamera();

  // text
  fs.vtk.text->SetInput( "...update view with mouse/keyboard, then press <space> to continue " );

  // renders window
  fs.vtk.renWin->Render();

  // window interactor
  // below call might crash if window is rendered before interactor is initialized
  // error in X11: GLXBadCurrentWindow
  interactor_usage();
  fs.vtk.iren->Start();
  
  // starts timer
  fs.timer->StartTimer();  
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int* it_h,float* time_h, float* data) {

  double xyz[3];
  double bounds[2];
  double min,max;
  
  static int VERBOSE = 1;
  
  int it = *it_h;
  float time = *time_h;
  
  if (VERBOSE) printf("     visual: it = %d time = %f \n",it,time);

  // time for calculating new wavefield
  fs.timer->StopTimer();
  double timeInSeconds = fs.timer->GetElapsedTime();
  fs.timer->StartTimer();
  
  // frames per second (based on computing new wavefield
  double fps = 1.0/timeInSeconds;
  // updates time string
  char inputString[180];
  sprintf(inputString,"time step: %d \t time: %6.3f \t fps: %f",it,time,fps);
  fs.vtk.text->SetInput(inputString);

  // updates data
  if (SHOW_VOLUMEDATA == 1) {
    for(int i = 0;i<fs.volume.np;i++) {
      fs.vtk.data_array3D->SetValue(i,data[i]);
    }
    // mark as modified to update rendering
    fs.vtk.data_array3D->Modified();

    // Find min and max
    fs.vtk.data_array3D->GetValueRange(bounds);
    min = bounds[0];
    max = bounds[1];
    
    // adjusts color maximum
    if (fs.vtk.colorAdjustOn) {
      if (gcolor_max < 0.0 ) gcolor_max = 1.e-10;
      gcolor_max = max;
      gcolor_incr = 0.1 * max;
      if (fs.vtk.icolor >= 2 ) gcolor_min = 1.e-3*gcolor_max;
      fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
      fs.vtk.lut->Build();
    }
  }

  // updates window
  fs.vtk.renWin->Render();

  // saves snapshot image
  if (fs.vtk.jpegImageOn) {
    save_snapshot_jpg( it );
  }

  // time for rendering
  fs.timer->StopTimer();
  double time_renderer = fs.timer->GetElapsedTime();
  fs.timer->StartTimer();

  if (VERBOSE) {
    printf("     visual: %s \n",inputString);
    printf("     timer : time for rendering = %f (s) \n", time_renderer);
    printf("     data  : min = %e max = %e \n\n",min,max);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)() {

  // text
  fs.vtk.text->SetInput( "...change view with mouse/keyboard, then press <space> to finish " );
  
  // render window
  fs.vtk.renWin->Render();

  // window interactor
  // below call might crash if window is rendered before interactor is initialized
  // error in X11: GLXBadCurrentWindow
  interactor_usage();
  fs.vtk.iren->Start();

  fs.timer->StopTimer();
  fs.timer->Delete();
  
  // free arrays
  if (SHOW_FREESURFACE == 1) {
    fs.vtk.points2D->Delete();
    fs.vtk.data_array2D->Delete();
    fs.vtk.cells2D->Delete();

    fs.vtk.lut2D->Delete();

    // contour
    fs.vtk.contourActor->Delete();
    fs.vtk.contourMapper->Delete();
    fs.vtk.contour->Delete();
    
    fs.vtk.actor2D->Delete();
    fs.vtk.mapMesh2D->Delete();    
  }

  if (SHOW_VOLUMEDATA == 1) {
    fs.vtk.points3D->Delete();
    fs.vtk.data_array3D->Delete();
    fs.vtk.cells->Delete();

    fs.vtk.legendcolor->Delete();
    fs.vtk.lut->Delete();

    fs.vtk.clipPlane1->Delete();
    fs.vtk.clipPlane2->Delete();
    
    fs.vtk.clippersurface->Delete();
    fs.vtk.merger->Delete();
    
    fs.vtk.clip1->Delete();
    fs.vtk.clip1m->Delete();
    fs.vtk.clip2->Delete();
    
    fs.vtk.actor3D->Delete();
    fs.vtk.mapMesh3D->Delete();
  }

  // source sphere
  fs.vtk.actorSphere->Delete();
  fs.vtk.mapperSphere->Delete();
  fs.vtk.sphere->Delete();

  fs.vtk.textGPU->Delete();
  fs.vtk.text->Delete();

  fs.vtk.iren->Delete();
  fs.vtk.renWin->Delete();
  fs.vtk.ren->Delete();
  fs.vtk.camera->Delete();  
}


#endif // HAVE_VTK

