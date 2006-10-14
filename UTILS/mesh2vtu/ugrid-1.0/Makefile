CXX=g++
VTK=/home/nschwarz/software/vtk42
VTK=/opt/seismo-util
INCLUDES=-I$(VTK)/include
CPPFLAGS=$(INCLUDES) -Wno-deprecated
LIBS=-L$(VTK)/lib/vtk -lvtkCommon -lvtkIO -lvtkRendering -lvtkFiltering -lvtkGraphics
OBJS=ugrid2.o
SRC=ugrid2.cxx
PROG=ugrid2

$(PROG): $(OBJS)
	$(CXX) $(CPPFLAGS) -o $(PROG) $(OBJS) $(LIBS)

$(OBJS): $(SRC)
	$(CXX) $(CPPFLAGS) -c $(SRC) -o $(OBJS)

append_data: append_data.cxx
	$(CXX) $(CPPFLAGS) -o append_data append_data.cxx $(LIBS)

mesh2vtu: mesh2vtu.o
	$(CXX) $(CPPFLAGS) -o mesh2vtu mesh2vtu.o $(LIBS)

mesh2vtu.o: mesh2vtu.cxx
	$(CXX) $(CPPFLAGS) -c mesh2vtu.cxx -o mesh2vtu.o
