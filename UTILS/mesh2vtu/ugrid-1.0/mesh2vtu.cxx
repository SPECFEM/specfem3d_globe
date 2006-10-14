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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <vtk/vtkFloatArray.h>
#include <vtk/vtkPoints.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkXMLPolyDataWriter.h>
#include <vtk/vtkUnstructuredGridToPolyDataFilter.h>
#include <vtk/vtkDelaunay3D.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkPointSet.h>
#include <vtk/vtkHexahedron.h>


int main(int argc, char** argv) {

  if (argc < 3) {
    printf("Usage: ugrid input_file output_file\n");
    return 0;
  }

  float xyz[3];
  int cell[8];
  FILE *file;
  int i, j;
  int npts, ncells;
  int pid[8];

  int fd;
  
  if((fd = open(argv[1], O_RDONLY)) == -1) {
    printf("Error opening file: %s.\n", argv[1]);
    return 0;
  }

  if(read(fd, &npts, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in points): %s\n", argv[1]);
  }
  
  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New();
  float *xV = new float[npts];
  float *yV = new float[npts];
  float *zV = new float[npts];
  float *sV = new float[npts];

  vtkPoints *newPts = vtkPoints::New();
  vtkFloatArray *newScalars = vtkFloatArray::New();
  printf("mesh2vtu: Reading in points: %d\n", npts);
  for (i = 0 ; i < npts ; i++)
    {
      read(fd, &xV[i], sizeof(float));
      read(fd, &yV[i], sizeof(float));
      read(fd, &zV[i], sizeof(float));
      read(fd, &sV[i], sizeof(float));
      xyz[0] = xV[i]; 
      xyz[1] = yV[i]; 
      xyz[2] = zV[i];
      newPts -> InsertPoint(i, xyz);
      newScalars -> InsertValue(i, sV[i]);
    }

  vtkCellArray *cells = vtkCellArray::New();
  if(read(fd, &ncells, sizeof(int)) != sizeof(int)) {
    printf("Bad read on file (in cells): %s\n", argv[1]);
  }
  printf("mesh2vtu: Reading in cells: %d\n", ncells);  
  int *cellTypes = new int[ncells];
  vtkHexahedron *hex = vtkHexahedron::New();
  hex->GetPointIds()->SetNumberOfIds(8);

  for(i = 0; i < ncells; i++) {
    for(j = 0; j < 8; j++) {
      read(fd, &cell[j], sizeof(int));
      hex->GetPointIds()->SetId(j,cell[j]);
    }
    cells->InsertNextCell(hex);
    cellTypes[i] = hex->GetCellType();
  }
  
  close(fd);
  
  dataSet -> SetPoints(newPts);
  dataSet -> GetPointData() -> SetScalars(newScalars);
  dataSet -> SetCells(cellTypes, cells);
  
  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  writer -> SetInput(dataSet);
  writer -> SetFileName(argv[2]);
  writer -> Write();

  writer -> Delete();
  newPts -> Delete();
  newScalars -> Delete();
  dataSet -> Delete();
  cells -> Delete();
  
  //  printf("Done.\n");
 
  return 0;

}
