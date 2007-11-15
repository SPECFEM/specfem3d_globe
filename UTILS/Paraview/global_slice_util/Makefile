F90 = ifort
F90_FLAGS = -O3

LIB = 

SUB_MOD = sub_slice_number
SRC_ALL = global_slice_number normal_plane make_az_stations

MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = ..

SUB_MOD_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(SUB_MOD))
SRC_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(SRC_ALL))

F90_OBJ = $(SUB_MOD_OBJ) $(SRC_OBJ)

all : $(SRC_ALL)

$(SRC_ALL): % : %.f90 $(SUB_MOD_OBJ) constants.h
	$(F90) -o $(BIN_DIR)/x$@ $(F90_FLAGS) $*.f90 -module $(MOD_DIR) $(SUB_MOD_OBJ) $(LIB)

$(SUB_MOD_OBJ): $(OBJ_DIR)/%.o : %.f90
	$(F90) -c -o $(OBJ_DIR)/$*.o $*.f90 -module $(MOD_DIR)

clean: 
	\rm -f $(MOD_DIR)/*.mod $(OBJ_DIR)/*.o *.o *.mod *~ garc_station.txt  xsection_translate.txt xglobal_slice_number xmake_az_stations xnormal_plane
