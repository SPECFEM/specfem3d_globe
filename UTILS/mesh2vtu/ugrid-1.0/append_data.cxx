
#include <stdio.h>

#include <vtk/vtkXMLUnstructuredGridReader.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkAppendFilter.h>

int
main(int argc, char *argv[]) {
  int i;
  vtkXMLUnstructuredGridReader *reader;
  vtkXMLUnstructuredGridWriter *writer = vtkXMLUnstructuredGridWriter::New();
  vtkAppendFilter *filter = vtkAppendFilter::New();

  for(i = 2; i < argc; i++) {
    printf("Appending: %s\n", argv[i]);
    reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(argv[i]);  
    filter->AddInput((vtkDataSet*)reader->GetOutput());
    reader->Delete();
  }
  printf("Writing: %s\n", argv[1]);
  writer->SetInput(filter->GetOutput());
  writer->SetFileName(argv[1]);
  writer->Write();
  
  return 0;
}
