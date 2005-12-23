# Process this file with Pyrex to produce Specfem3DGlobeCode.c


"""Python bindings for the SPECFEM3D Global Solver."""


# include 'config.h' in order to get the definitions of FC_FUNC and FC_FUNC_
cdef extern from "config.h":
    pass


cdef extern from "string.h":
    char *strncpy(char *, char *, int)


cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *, int)


# In the future, this could be passed through the Fortran layer as an
# opaque context argument.
component = None


# replacements for Fortran functions

cdef public void read_value_integer "FC_FUNC_(read_value_integer, READ_VALUE_INTEGER)" (int *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = component.readValue(attrName)

cdef public void read_value_double_precision "FC_FUNC_(read_value_double_precision, READ_VALUE_DOUBLE_PRECISION)" (double *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = component.readValue(attrName)

cdef public void read_value_logical "FC_FUNC_(read_value_logical, READ_VALUE_LOGICAL)" (int *value, char *name, int nameLen) except *:
    attrName = PyString_FromStringAndSize(name, nameLen)
    value[0] = component.readValue(attrName)

cdef public void read_value_string "FC_FUNC_(read_value_string, READ_VALUE_STRING)" (char *value, char *name, int valueLen, int nameLen) except *:
    cdef char *vp
    cdef int vl, i
    attrName = PyString_FromStringAndSize(name, nameLen)
    v = component.readValue(attrName)
    vl = len(v)
    if vl > valueLen:
        raise ValueError("%s value '%s' is too long (%d bytes) for destination Fortran buffer (%d bytes)" % (attrName, v, vl, valueLen))
    vp = v
    strncpy(value, vp, vl)
    for i from vl <= i < valueLen:
        value[i] = c' '
    return

cdef public void open_parameter_file "FC_FUNC_(open_parameter_file, OPEN_PARAMETER_FILE)" () except *:
    return

cdef public void close_parameter_file "FC_FUNC_(close_parameter_file, CLOSE_PARAMETER_FILE)" () except *:
    return


# external Fortran functions

cdef extern void meshfem3D_f "FC_FUNC(meshfem3d, MESHFEM3D)" () except *
def meshfem3D(arg):
    """Run the SPECFEM3D Global Mesher."""
    global component
    component = arg
    meshfem3D_f()

cdef extern void specfem3D_f "FC_FUNC(specfem3d, SPECFEM3D)" () except *
def specfem3D(arg):
    """Run the SPECFEM3D Global Solver."""
    global component
    component = arg
    specfem3D_f()

cdef extern void check_parameters_f "FC_FUNC_(check_parameters, CHECK_PARAMETERS)" () except *
def check_parameters(arg):
    global component
    component = arg
    check_parameters_f()


# end of file
