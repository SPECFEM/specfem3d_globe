//-----------------------------------------------------------------------------
// Program:     ugrid
// Description: Converts x, y, z, s points to VTK XML Unstructured Grid
//              via 3D Delaunay triangulation.
// File:        ugrid.cxx
// Author:      Nicholas Schwarz, schwarz@evl.uic.edu
//              Electronic Visualization Laboratory
//              University of Illinois at Chicago
// Date:        4 June 2004
//-----------------------------------------------------------------------------

#include <stdio.h>

#include <vtk/vtkFloatArray.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkXMLPolyDataWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkDelaunay3D.h>
#include <vtk/vtkPointSet.h>

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: ugrid input_file output_file\n");
    return 0;
  }

  printf("This may take a while...\n");

  float xyz[3];
  FILE *file;
  int i;
  int npts, ncells;
  int pid[8];
  
  if ((file = fopen(argv[1], "r")) == 0)
    {
      printf("Error opening file.\n");
      return 0;
    }

  fscanf(file, "%d", &npts);

  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  float *xV = new float[npts];
  float *yV = new float[npts];
  float *zV = new float[npts];
  float *sV = new float[npts];

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();

  for (i = 0 ; i < npts ; i++)
    {
      fscanf(file, "%f", &xV[i]);
      fscanf(file, "%f", &yV[i]);
      fscanf(file, "%f", &zV[i]);
      fscanf(file, "%f", &sV[i]);
      xyz[0] = xV[i]; 
      xyz[1] = yV[i]; 
      xyz[2] = zV[i];
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, sV[i]);
    }

  fclose(file);

  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);

  vtkDelaunay3D* del = vtkDelaunay3D::New();
  del -> SetInput((vtkPointSet*)dataSet);

  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  writer -> SetInput(del -> GetOutput());
  writer -> SetFileName(argv[2]);
  writer -> Write();

  del -> Delete();
  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();

  printf("Done.\n");
 
  return 0;

}
