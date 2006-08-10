# Process this file with Pyrex to produce PyxMeshfem.c


"""Python bindings for the SPECFEM3D Global Solver."""


# include 'config.h' in order to get the definitions of FC_FUNC and FC_FUNC_
cdef extern from "config.h":
    pass


# external Fortran functions

cdef extern void meshfem3D_f "FC_FUNC(meshfem3d, MESHFEM3D)" () except *
def meshfem3D(arg):
    """Run the SPECFEM3D Global Mesher."""
    import PyxParameters
    PyxParameters.component = arg
    meshfem3D_f()


cdef extern void read_params_and_create_movie_f "FC_FUNC_(read_params_and_create_movie, READ_PARAMS_AND_CREATE_MOVIE)" () except *
def read_params_and_create_movie(arg):
    """Create a movie."""
    import PyxParameters
    PyxParameters.component = arg
    read_params_and_create_movie_f()


# end of file
