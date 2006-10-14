/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: finance.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <vtk/vtkActor.h>
#include <vtk/vtkAxes.h>
#include <vtk/vtkContourFilter.h>
#include <vtk/vtkDataSet.h>
#include <vtk/vtkFloatArray.h>
#include <vtk/vtkGaussianSplatter.h>
#include <vtk/vtkImageData.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkPolyDataMapper.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkRenderWindow.h>
#include <vtk/vtkRenderWindowInteractor.h>
#include <vtk/vtkRenderer.h>
#include <vtk/vtkTubeFilter.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridReader.h>

int main( int argc, char *argv[] )
{
  double bounds[6];
  vtkDataSet *dataSet;
  
  vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
  reader -> SetFileName(argv[1]);
  reader -> Update();

  // construct pipeline for original population
  vtkGaussianSplatter *popSplatter = vtkGaussianSplatter::New();
    popSplatter->SetInput(reader -> GetOutput());
    popSplatter->SetSampleDimensions(50,50,50);
    popSplatter->SetRadius(0.05);
    popSplatter->ScalarWarpingOff();
  vtkContourFilter *popSurface = vtkContourFilter::New();
    popSurface->SetInput(reader->GetOutput());
    popSurface->SetValue(0,0.01);
  vtkPolyDataMapper *popMapper = vtkPolyDataMapper::New();
    popMapper->SetInput(popSurface->GetOutput());
    popMapper->ScalarVisibilityOff();
  vtkActor *popActor = vtkActor::New();
    popActor->SetMapper(popMapper);
    popActor->GetProperty()->SetOpacity(0.3);
    popActor->GetProperty()->SetColor(.9,.9,.9);

  // construct pipeline for delinquent population
  vtkGaussianSplatter *lateSplatter = vtkGaussianSplatter::New();
    lateSplatter->SetInput(reader -> GetOutput());
    lateSplatter->SetSampleDimensions(50,50,50);
    lateSplatter->SetRadius(0.05);
    lateSplatter->SetScaleFactor(0.005);
  vtkContourFilter *lateSurface = vtkContourFilter::New();
    lateSurface->SetInput(reader->GetOutput());
    lateSurface->SetValue(0,10.0);
  vtkPolyDataMapper *lateMapper = vtkPolyDataMapper::New();
    lateMapper->SetInput(lateSurface->GetOutput());
    lateMapper->ScalarVisibilityOn();
  vtkActor *lateActor = vtkActor::New();
    lateActor->SetMapper(lateMapper);
    lateActor->GetProperty()->SetColor(1.0,0.0,0.0);

    /*
  // create axes
  popSplatter->Update();
  popSplatter->GetOutput()->GetBounds(bounds);
  vtkAxes *axes = vtkAxes::New();
    axes->SetOrigin(bounds[0], bounds[2], bounds[4]);
    axes->SetScaleFactor(popSplatter->GetOutput()->GetLength()/5);
  vtkTubeFilter *axesTubes = vtkTubeFilter::New();
    axesTubes->SetInput(axes->GetOutput());
    axesTubes->SetRadius(axes->GetScaleFactor()/25.0);
    axesTubes->SetNumberOfSides(6);
  vtkPolyDataMapper *axesMapper = vtkPolyDataMapper::New();
    axesMapper->SetInput(axesTubes->GetOutput());
  vtkActor *axesActor = vtkActor::New();
    axesActor->SetMapper(axesMapper);
    */

  // graphics stuff
  vtkRenderer *renderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(renderer);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

  // read data  //set up renderer
  renderer->AddActor(lateActor);
  //  renderer->AddActor(axesActor);
  renderer->AddActor(popActor);
  renderer->SetBackground(1,1,1);
  renWin->SetSize(300,300);

  // interact with data
  iren->Initialize();

  renWin->Render();
  iren->Start();

  // Clean up
  renderer->Delete();
  renWin->Delete();
  iren->Delete();
  popSplatter->Delete();
  popSurface->Delete();
  popMapper->Delete();
  popActor->Delete();
  lateSplatter->Delete();
  lateSurface->Delete();
  lateMapper->Delete();
  lateActor->Delete();
  //  axes->Delete();
  //  axesTubes->Delete();
  //  axesMapper->Delete();
  //  axesActor->Delete();
  //dataSet->Delete();

  return 0;
}


