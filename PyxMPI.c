/* Generated by Pyrex 0.9.3 on Mon Jun 12 00:04:05 2006 */

#include "Python.h"
#include "structmember.h"
#ifndef PY_LONG_LONG
  #define PY_LONG_LONG LONG_LONG
#endif
#include "stdlib.h"
#include "mpi.h"


typedef struct {PyObject **p; char *s;} __Pyx_InternTabEntry; /*proto*/
typedef struct {PyObject **p; char *s; long n;} __Pyx_StringTabEntry; /*proto*/
static PyObject *__Pyx_UnpackItem(PyObject *, int); /*proto*/
static int __Pyx_EndUnpack(PyObject *, int); /*proto*/
static int __Pyx_PrintItem(PyObject *); /*proto*/
static int __Pyx_PrintNewline(void); /*proto*/
static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb); /*proto*/
static void __Pyx_ReRaise(void); /*proto*/
static PyObject *__Pyx_Import(PyObject *name, PyObject *from_list); /*proto*/
static PyObject *__Pyx_GetExcValue(void); /*proto*/
static int __Pyx_ArgTypeTest(PyObject *obj, PyTypeObject *type, int none_allowed, char *name); /*proto*/
static int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type); /*proto*/
static int __Pyx_GetStarArgs(PyObject **args, PyObject **kwds, char *kwd_list[], int nargs, PyObject **args2, PyObject **kwds2); /*proto*/
static void __Pyx_WriteUnraisable(char *name); /*proto*/
static void __Pyx_AddTraceback(char *funcname); /*proto*/
static PyTypeObject *__Pyx_ImportType(char *module_name, char *class_name, long size);  /*proto*/
static int __Pyx_SetVtable(PyObject *dict, void *vtable); /*proto*/
static int __Pyx_GetVtable(PyObject *dict, void *vtabptr); /*proto*/
static PyObject *__Pyx_CreateClass(PyObject *bases, PyObject *dict, PyObject *name, char *modname); /*proto*/
static int __Pyx_InternStrings(__Pyx_InternTabEntry *t); /*proto*/
static int __Pyx_InitStrings(__Pyx_StringTabEntry *t); /*proto*/
static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name); /*proto*/

static PyObject *__pyx_m;
static PyObject *__pyx_b;
static int __pyx_lineno;
static char *__pyx_filename;
staticforward char **__pyx_f;

/* Declarations from mpi */


/* Declarations from PyxMPI */

staticforward PyTypeObject __pyx_type_6PyxMPI_MPI_Comm;

struct __pyx_obj_6PyxMPI_MPI_Comm {
  PyObject_HEAD
  MPI_Comm comm;
};

static PyTypeObject *__pyx_ptype_6PyxMPI_MPI_Comm = 0;

/* Implementation of PyxMPI */


static PyObject *__pyx_n_mpi;
static PyObject *__pyx_n_MPI_COMM_WORLD;
static PyObject *__pyx_n_MPI_Error;
static PyObject *__pyx_n_MPI_Init;
static PyObject *__pyx_n_MPI_Finalize;
static PyObject *__pyx_n_MPI_Comm_rank;
static PyObject *__pyx_n_EnvironmentError;

static int __pyx_f_6PyxMPI_8MPI_Comm___init__(PyObject *__pyx_v_self, PyObject *__pyx_args, PyObject *__pyx_kwds); /*proto*/
static int __pyx_f_6PyxMPI_8MPI_Comm___init__(PyObject *__pyx_v_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  int __pyx_r;
  static char *__pyx_argnames[] = {0};
  if (!PyArg_ParseTupleAndKeywords(__pyx_args, __pyx_kwds, "", __pyx_argnames)) return -1;
  Py_INCREF(__pyx_v_self);

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":17 */
  ((struct __pyx_obj_6PyxMPI_MPI_Comm *)__pyx_v_self)->comm = MPI_COMM_WORLD;

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1:;
  __Pyx_AddTraceback("PyxMPI.MPI_Comm.__init__");
  __pyx_r = -1;
  __pyx_L0:;
  Py_DECREF(__pyx_v_self);
  return __pyx_r;
}

static PyObject *__pyx_n_len;
static PyObject *__pyx_n_append;

static PyObject *__pyx_f_6PyxMPI_MPI_Init(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); /*proto*/
static PyObject *__pyx_f_6PyxMPI_MPI_Init(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_argv = 0;
  int __pyx_v_error;
  int __pyx_v_cargc;
  int __pyx_v_i;
  char (*(*__pyx_v_cargv));
  char (*(*__pyx_v_mycargv));
  PyObject *__pyx_v_myargv;
  PyObject *__pyx_v_arg;
  PyObject *__pyx_r;
  PyObject *__pyx_1 = 0;
  PyObject *__pyx_2 = 0;
  PyObject *__pyx_3 = 0;
  int __pyx_4;
  char (*__pyx_5);
  static char *__pyx_argnames[] = {"argv",0};
  if (!PyArg_ParseTupleAndKeywords(__pyx_args, __pyx_kwds, "O", __pyx_argnames, &__pyx_v_argv)) return 0;
  Py_INCREF(__pyx_v_argv);
  __pyx_v_myargv = Py_None; Py_INCREF(__pyx_v_myargv);
  __pyx_v_arg = Py_None; Py_INCREF(__pyx_v_arg);

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":30 */
  __pyx_1 = PyList_New(0); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 30; goto __pyx_L1;}
  Py_DECREF(__pyx_v_myargv);
  __pyx_v_myargv = __pyx_1;
  __pyx_1 = 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":33 */
  __pyx_1 = __Pyx_GetName(__pyx_b, __pyx_n_len); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 33; goto __pyx_L1;}
  __pyx_2 = PyTuple_New(1); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 33; goto __pyx_L1;}
  Py_INCREF(__pyx_v_argv);
  PyTuple_SET_ITEM(__pyx_2, 0, __pyx_v_argv);
  __pyx_3 = PyObject_CallObject(__pyx_1, __pyx_2); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 33; goto __pyx_L1;}
  Py_DECREF(__pyx_1); __pyx_1 = 0;
  Py_DECREF(__pyx_2); __pyx_2 = 0;
  __pyx_4 = PyInt_AsLong(__pyx_3); if (PyErr_Occurred()) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 33; goto __pyx_L1;}
  Py_DECREF(__pyx_3); __pyx_3 = 0;
  __pyx_v_cargc = __pyx_4;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":34 */
  __pyx_v_cargv = ((char (*(*)))malloc(((__pyx_v_cargc + 1) * (sizeof(char (*))))));

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":35 */
  for (__pyx_v_i = 0; __pyx_v_i < __pyx_v_cargc; ++__pyx_v_i) {

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":36 */
    __pyx_1 = PyInt_FromLong(__pyx_v_i); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 36; goto __pyx_L1;}
    __pyx_2 = PyObject_GetItem(__pyx_v_argv, __pyx_1); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 36; goto __pyx_L1;}
    Py_DECREF(__pyx_1); __pyx_1 = 0;
    Py_DECREF(__pyx_v_arg);
    __pyx_v_arg = __pyx_2;
    __pyx_2 = 0;

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":37 */
    __pyx_3 = PyObject_GetAttr(__pyx_v_myargv, __pyx_n_append); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 37; goto __pyx_L1;}
    __pyx_1 = PyTuple_New(1); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 37; goto __pyx_L1;}
    Py_INCREF(__pyx_v_arg);
    PyTuple_SET_ITEM(__pyx_1, 0, __pyx_v_arg);
    __pyx_2 = PyObject_CallObject(__pyx_3, __pyx_1); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 37; goto __pyx_L1;}
    Py_DECREF(__pyx_3); __pyx_3 = 0;
    Py_DECREF(__pyx_1); __pyx_1 = 0;
    Py_DECREF(__pyx_2); __pyx_2 = 0;

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":38 */
    __pyx_5 = PyString_AsString(__pyx_v_arg); if (PyErr_Occurred()) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 38; goto __pyx_L1;}
    (__pyx_v_cargv[__pyx_v_i]) = __pyx_5;
    __pyx_L2:;
  }
  __pyx_L3:;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":39 */
  (__pyx_v_cargv[__pyx_v_cargc]) = 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":42 */
  __pyx_v_mycargv = __pyx_v_cargv;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":43 */
  __pyx_v_error = MPI_Init((&__pyx_v_cargc),(&__pyx_v_cargv));

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":44 */
  __pyx_4 = (__pyx_v_error != MPI_SUCCESS);
  if (__pyx_4) {

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":45 */
    free(__pyx_v_mycargv);

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":46 */
    __pyx_3 = __Pyx_GetName(__pyx_m, __pyx_n_MPI_Error); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 46; goto __pyx_L1;}
    __pyx_1 = PyInt_FromLong(__pyx_v_error); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 46; goto __pyx_L1;}
    __Pyx_Raise(__pyx_3, __pyx_1, 0);
    Py_DECREF(__pyx_3); __pyx_3 = 0;
    Py_DECREF(__pyx_1); __pyx_1 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 46; goto __pyx_L1;}
    goto __pyx_L4;
  }
  __pyx_L4:;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":49 */
  if (PySequence_DelSlice(__pyx_v_argv, 0, 0x7fffffff) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 49; goto __pyx_L1;}

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":50 */
  for (__pyx_v_i = 0; __pyx_v_i < __pyx_v_cargc; ++__pyx_v_i) {

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":51 */
    __pyx_2 = PyObject_GetAttr(__pyx_v_argv, __pyx_n_append); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 51; goto __pyx_L1;}
    __pyx_3 = PyString_FromString((__pyx_v_cargv[__pyx_v_i])); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 51; goto __pyx_L1;}
    __pyx_1 = PyTuple_New(1); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 51; goto __pyx_L1;}
    PyTuple_SET_ITEM(__pyx_1, 0, __pyx_3);
    __pyx_3 = 0;
    __pyx_3 = PyObject_CallObject(__pyx_2, __pyx_1); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 51; goto __pyx_L1;}
    Py_DECREF(__pyx_2); __pyx_2 = 0;
    Py_DECREF(__pyx_1); __pyx_1 = 0;
    Py_DECREF(__pyx_3); __pyx_3 = 0;
    __pyx_L5:;
  }
  __pyx_L6:;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":52 */
  free(__pyx_v_mycargv);

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":54 */
  __pyx_r = Py_None; Py_INCREF(__pyx_r);
  goto __pyx_L0;

  __pyx_r = Py_None; Py_INCREF(__pyx_r);
  goto __pyx_L0;
  __pyx_L1:;
  Py_XDECREF(__pyx_1);
  Py_XDECREF(__pyx_2);
  Py_XDECREF(__pyx_3);
  __Pyx_AddTraceback("PyxMPI.MPI_Init");
  __pyx_r = 0;
  __pyx_L0:;
  Py_DECREF(__pyx_v_myargv);
  Py_DECREF(__pyx_v_arg);
  Py_DECREF(__pyx_v_argv);
  return __pyx_r;
}

static PyObject *__pyx_f_6PyxMPI_MPI_Finalize(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); /*proto*/
static PyObject *__pyx_f_6PyxMPI_MPI_Finalize(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  int __pyx_v_error;
  PyObject *__pyx_r;
  int __pyx_1;
  PyObject *__pyx_2 = 0;
  PyObject *__pyx_3 = 0;
  static char *__pyx_argnames[] = {0};
  if (!PyArg_ParseTupleAndKeywords(__pyx_args, __pyx_kwds, "", __pyx_argnames)) return 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":59 */
  __pyx_v_error = MPI_Finalize();

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":60 */
  __pyx_1 = (__pyx_v_error != MPI_SUCCESS);
  if (__pyx_1) {

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":61 */
    __pyx_2 = __Pyx_GetName(__pyx_m, __pyx_n_MPI_Error); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 61; goto __pyx_L1;}
    __pyx_3 = PyInt_FromLong(__pyx_v_error); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 61; goto __pyx_L1;}
    __Pyx_Raise(__pyx_2, __pyx_3, 0);
    Py_DECREF(__pyx_2); __pyx_2 = 0;
    Py_DECREF(__pyx_3); __pyx_3 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 61; goto __pyx_L1;}
    goto __pyx_L2;
  }
  __pyx_L2:;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":62 */
  __pyx_r = Py_None; Py_INCREF(__pyx_r);
  goto __pyx_L0;

  __pyx_r = Py_None; Py_INCREF(__pyx_r);
  goto __pyx_L0;
  __pyx_L1:;
  Py_XDECREF(__pyx_2);
  Py_XDECREF(__pyx_3);
  __Pyx_AddTraceback("PyxMPI.MPI_Finalize");
  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}

static PyObject *__pyx_f_6PyxMPI_MPI_Comm_rank(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); /*proto*/
static PyObject *__pyx_f_6PyxMPI_MPI_Comm_rank(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_comm = 0;
  int __pyx_v_error;
  int __pyx_v_rank;
  struct __pyx_obj_6PyxMPI_MPI_Comm *__pyx_v_c_comm;
  PyObject *__pyx_r;
  int __pyx_1;
  PyObject *__pyx_2 = 0;
  PyObject *__pyx_3 = 0;
  static char *__pyx_argnames[] = {"comm",0};
  if (!PyArg_ParseTupleAndKeywords(__pyx_args, __pyx_kwds, "O", __pyx_argnames, &__pyx_v_comm)) return 0;
  Py_INCREF(__pyx_v_comm);
  ((PyObject*)__pyx_v_c_comm) = Py_None; Py_INCREF(((PyObject*)__pyx_v_c_comm));

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":69 */
  if (!__Pyx_TypeTest(__pyx_v_comm, __pyx_ptype_6PyxMPI_MPI_Comm)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 69; goto __pyx_L1;}
  Py_INCREF(__pyx_v_comm);
  Py_DECREF(((PyObject *)__pyx_v_c_comm));
  ((PyObject *)__pyx_v_c_comm) = __pyx_v_comm;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":70 */
  __pyx_v_error = MPI_Comm_rank(__pyx_v_c_comm->comm,(&__pyx_v_rank));

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":71 */
  __pyx_1 = (__pyx_v_error != MPI_SUCCESS);
  if (__pyx_1) {

    /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":72 */
    __pyx_2 = __Pyx_GetName(__pyx_m, __pyx_n_MPI_Error); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 72; goto __pyx_L1;}
    __pyx_3 = PyInt_FromLong(__pyx_v_error); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 72; goto __pyx_L1;}
    __Pyx_Raise(__pyx_2, __pyx_3, 0);
    Py_DECREF(__pyx_2); __pyx_2 = 0;
    Py_DECREF(__pyx_3); __pyx_3 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 72; goto __pyx_L1;}
    goto __pyx_L2;
  }
  __pyx_L2:;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":73 */
  __pyx_2 = PyInt_FromLong(__pyx_v_rank); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 73; goto __pyx_L1;}
  __pyx_r = __pyx_2;
  __pyx_2 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; Py_INCREF(__pyx_r);
  goto __pyx_L0;
  __pyx_L1:;
  Py_XDECREF(__pyx_2);
  Py_XDECREF(__pyx_3);
  __Pyx_AddTraceback("PyxMPI.MPI_Comm_rank");
  __pyx_r = 0;
  __pyx_L0:;
  Py_DECREF(__pyx_v_c_comm);
  Py_DECREF(__pyx_v_comm);
  return __pyx_r;
}

static __Pyx_InternTabEntry __pyx_intern_tab[] = {
  {&__pyx_n_EnvironmentError, "EnvironmentError"},
  {&__pyx_n_MPI_COMM_WORLD, "MPI_COMM_WORLD"},
  {&__pyx_n_MPI_Comm_rank, "MPI_Comm_rank"},
  {&__pyx_n_MPI_Error, "MPI_Error"},
  {&__pyx_n_MPI_Finalize, "MPI_Finalize"},
  {&__pyx_n_MPI_Init, "MPI_Init"},
  {&__pyx_n_append, "append"},
  {&__pyx_n_len, "len"},
  {&__pyx_n_mpi, "mpi"},
  {0, 0}
};

static PyObject *__pyx_tp_new_6PyxMPI_MPI_Comm(PyTypeObject *t, PyObject *a, PyObject *k) {
  PyObject *o = (*t->tp_alloc)(t, 0);
  struct __pyx_obj_6PyxMPI_MPI_Comm *p = (struct __pyx_obj_6PyxMPI_MPI_Comm *)o;
  return o;
}

static void __pyx_tp_dealloc_6PyxMPI_MPI_Comm(PyObject *o) {
  struct __pyx_obj_6PyxMPI_MPI_Comm *p = (struct __pyx_obj_6PyxMPI_MPI_Comm *)o;
  (*o->ob_type->tp_free)(o);
}

static int __pyx_tp_traverse_6PyxMPI_MPI_Comm(PyObject *o, visitproc v, void *a) {
  int e;
  struct __pyx_obj_6PyxMPI_MPI_Comm *p = (struct __pyx_obj_6PyxMPI_MPI_Comm *)o;
  return 0;
}

static int __pyx_tp_clear_6PyxMPI_MPI_Comm(PyObject *o) {
  struct __pyx_obj_6PyxMPI_MPI_Comm *p = (struct __pyx_obj_6PyxMPI_MPI_Comm *)o;
  return 0;
}

static struct PyMethodDef __pyx_methods_6PyxMPI_MPI_Comm[] = {
  {0, 0, 0, 0}
};

static PyNumberMethods __pyx_tp_as_number_MPI_Comm = {
  0, /*nb_add*/
  0, /*nb_subtract*/
  0, /*nb_multiply*/
  0, /*nb_divide*/
  0, /*nb_remainder*/
  0, /*nb_divmod*/
  0, /*nb_power*/
  0, /*nb_negative*/
  0, /*nb_positive*/
  0, /*nb_absolute*/
  0, /*nb_nonzero*/
  0, /*nb_invert*/
  0, /*nb_lshift*/
  0, /*nb_rshift*/
  0, /*nb_and*/
  0, /*nb_xor*/
  0, /*nb_or*/
  0, /*nb_coerce*/
  0, /*nb_int*/
  0, /*nb_long*/
  0, /*nb_float*/
  0, /*nb_oct*/
  0, /*nb_hex*/
  0, /*nb_inplace_add*/
  0, /*nb_inplace_subtract*/
  0, /*nb_inplace_multiply*/
  0, /*nb_inplace_divide*/
  0, /*nb_inplace_remainder*/
  0, /*nb_inplace_power*/
  0, /*nb_inplace_lshift*/
  0, /*nb_inplace_rshift*/
  0, /*nb_inplace_and*/
  0, /*nb_inplace_xor*/
  0, /*nb_inplace_or*/
  0, /*nb_floor_divide*/
  0, /*nb_true_divide*/
  0, /*nb_inplace_floor_divide*/
  0, /*nb_inplace_true_divide*/
};

static PySequenceMethods __pyx_tp_as_sequence_MPI_Comm = {
  0, /*sq_length*/
  0, /*sq_concat*/
  0, /*sq_repeat*/
  0, /*sq_item*/
  0, /*sq_slice*/
  0, /*sq_ass_item*/
  0, /*sq_ass_slice*/
  0, /*sq_contains*/
  0, /*sq_inplace_concat*/
  0, /*sq_inplace_repeat*/
};

static PyMappingMethods __pyx_tp_as_mapping_MPI_Comm = {
  0, /*mp_length*/
  0, /*mp_subscript*/
  0, /*mp_ass_subscript*/
};

static PyBufferProcs __pyx_tp_as_buffer_MPI_Comm = {
  0, /*bf_getreadbuffer*/
  0, /*bf_getwritebuffer*/
  0, /*bf_getsegcount*/
  0, /*bf_getcharbuffer*/
};

statichere PyTypeObject __pyx_type_6PyxMPI_MPI_Comm = {
  PyObject_HEAD_INIT(0)
  0, /*ob_size*/
  "PyxMPI.MPI_Comm", /*tp_name*/
  sizeof(struct __pyx_obj_6PyxMPI_MPI_Comm), /*tp_basicsize*/
  0, /*tp_itemsize*/
  __pyx_tp_dealloc_6PyxMPI_MPI_Comm, /*tp_dealloc*/
  0, /*tp_print*/
  0, /*tp_getattr*/
  0, /*tp_setattr*/
  0, /*tp_compare*/
  0, /*tp_repr*/
  &__pyx_tp_as_number_MPI_Comm, /*tp_as_number*/
  &__pyx_tp_as_sequence_MPI_Comm, /*tp_as_sequence*/
  &__pyx_tp_as_mapping_MPI_Comm, /*tp_as_mapping*/
  0, /*tp_hash*/
  0, /*tp_call*/
  0, /*tp_str*/
  0, /*tp_getattro*/
  0, /*tp_setattro*/
  &__pyx_tp_as_buffer_MPI_Comm, /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT|Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_BASETYPE, /*tp_flags*/
  0, /*tp_doc*/
  __pyx_tp_traverse_6PyxMPI_MPI_Comm, /*tp_traverse*/
  __pyx_tp_clear_6PyxMPI_MPI_Comm, /*tp_clear*/
  0, /*tp_richcompare*/
  0, /*tp_weaklistoffset*/
  0, /*tp_iter*/
  0, /*tp_iternext*/
  __pyx_methods_6PyxMPI_MPI_Comm, /*tp_methods*/
  0, /*tp_members*/
  0, /*tp_getset*/
  0, /*tp_base*/
  0, /*tp_dict*/
  0, /*tp_descr_get*/
  0, /*tp_descr_set*/
  0, /*tp_dictoffset*/
  __pyx_f_6PyxMPI_8MPI_Comm___init__, /*tp_init*/
  0, /*tp_alloc*/
  __pyx_tp_new_6PyxMPI_MPI_Comm, /*tp_new*/
  0, /*tp_free*/
  0, /*tp_is_gc*/
  0, /*tp_bases*/
  0, /*tp_mro*/
  0, /*tp_cache*/
  0, /*tp_subclasses*/
  0, /*tp_weaklist*/
};

static struct PyMethodDef __pyx_methods[] = {
  {"MPI_Init", (PyCFunction)__pyx_f_6PyxMPI_MPI_Init, METH_VARARGS|METH_KEYWORDS, 0},
  {"MPI_Finalize", (PyCFunction)__pyx_f_6PyxMPI_MPI_Finalize, METH_VARARGS|METH_KEYWORDS, 0},
  {"MPI_Comm_rank", (PyCFunction)__pyx_f_6PyxMPI_MPI_Comm_rank, METH_VARARGS|METH_KEYWORDS, 0},
  {0, 0, 0, 0}
};

DL_EXPORT(void) initPyxMPI(void); /*proto*/
DL_EXPORT(void) initPyxMPI(void) {
  PyObject *__pyx_1 = 0;
  PyObject *__pyx_2 = 0;
  PyObject *__pyx_3 = 0;
  __pyx_m = Py_InitModule4("PyxMPI", __pyx_methods, 0, 0, PYTHON_API_VERSION);
  if (!__pyx_m) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4; goto __pyx_L1;};
  __pyx_b = PyImport_AddModule("__builtin__");
  if (!__pyx_b) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4; goto __pyx_L1;};
  if (PyObject_SetAttrString(__pyx_m, "__builtins__", __pyx_b) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4; goto __pyx_L1;};
  if (__Pyx_InternStrings(__pyx_intern_tab) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4; goto __pyx_L1;};
  if (PyType_Ready(&__pyx_type_6PyxMPI_MPI_Comm) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 12; goto __pyx_L1;}
  if (PyObject_SetAttrString(__pyx_m, "MPI_Comm", (PyObject *)&__pyx_type_6PyxMPI_MPI_Comm) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 12; goto __pyx_L1;}
  __pyx_ptype_6PyxMPI_MPI_Comm = &__pyx_type_6PyxMPI_MPI_Comm;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":20 */
  __pyx_1 = PyTuple_New(0); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
  __pyx_2 = PyObject_CallObject(((PyObject*)__pyx_ptype_6PyxMPI_MPI_Comm), __pyx_1); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
  Py_DECREF(__pyx_1); __pyx_1 = 0;
  if (PyObject_SetAttr(__pyx_m, __pyx_n_MPI_COMM_WORLD, __pyx_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 20; goto __pyx_L1;}
  Py_DECREF(__pyx_2); __pyx_2 = 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":23 */
  __pyx_1 = PyDict_New(); if (!__pyx_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 23; goto __pyx_L1;}
  __pyx_2 = __Pyx_GetName(__pyx_b, __pyx_n_EnvironmentError); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 23; goto __pyx_L1;}
  __pyx_3 = PyTuple_New(1); if (!__pyx_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 23; goto __pyx_L1;}
  PyTuple_SET_ITEM(__pyx_3, 0, __pyx_2);
  __pyx_2 = 0;
  __pyx_2 = __Pyx_CreateClass(__pyx_3, __pyx_1, __pyx_n_MPI_Error, "PyxMPI"); if (!__pyx_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 23; goto __pyx_L1;}
  Py_DECREF(__pyx_3); __pyx_3 = 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":24 */
  if (PyObject_SetAttr(__pyx_m, __pyx_n_MPI_Error, __pyx_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 23; goto __pyx_L1;}
  Py_DECREF(__pyx_2); __pyx_2 = 0;
  Py_DECREF(__pyx_1); __pyx_1 = 0;

  /* "/ibrixfs1/home/lstrand/dv/SPECFEM3D_GLOBE/PyxMPI.pyx":65 */
  return;
  __pyx_L1:;
  Py_XDECREF(__pyx_1);
  Py_XDECREF(__pyx_2);
  Py_XDECREF(__pyx_3);
  __Pyx_AddTraceback("PyxMPI");
}

static char *__pyx_filenames[] = {
  "PyxMPI.pyx",
};
statichere char **__pyx_f = __pyx_filenames;

/* Runtime support code */

static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name) {
    PyObject *result;
    result = PyObject_GetAttr(dict, name);
    if (!result)
        PyErr_SetObject(PyExc_NameError, name);
    return result;
}

static PyObject *__Pyx_CreateClass(
    PyObject *bases, PyObject *dict, PyObject *name, char *modname)
{
    PyObject *py_modname;
    PyObject *result = 0;
    
    py_modname = PyString_FromString(modname);
    if (!py_modname)
        goto bad;
    if (PyDict_SetItemString(dict, "__module__", py_modname) < 0)
        goto bad;
    result = PyClass_New(bases, dict, name);
bad:
    Py_XDECREF(py_modname);
    return result;
}

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb) {
    Py_XINCREF(type);
    Py_XINCREF(value);
    Py_XINCREF(tb);
    /* First, check the traceback argument, replacing None with NULL. */
    if (tb == Py_None) {
        Py_DECREF(tb);
        tb = 0;
    }
    else if (tb != NULL && !PyTraceBack_Check(tb)) {
        PyErr_SetString(PyExc_TypeError,
            "raise: arg 3 must be a traceback or None");
        goto raise_error;
    }
    /* Next, replace a missing value with None */
    if (value == NULL) {
        value = Py_None;
        Py_INCREF(value);
    }
    /* Next, repeatedly, replace a tuple exception with its first item */
    while (PyTuple_Check(type) && PyTuple_Size(type) > 0) {
        PyObject *tmp = type;
        type = PyTuple_GET_ITEM(type, 0);
        Py_INCREF(type);
        Py_DECREF(tmp);
    }
    if (PyString_Check(type))
        ;
    else if (PyClass_Check(type))
        ; /*PyErr_NormalizeException(&type, &value, &tb);*/
    else if (PyInstance_Check(type)) {
        /* Raising an instance.  The value should be a dummy. */
        if (value != Py_None) {
            PyErr_SetString(PyExc_TypeError,
              "instance exception may not have a separate value");
            goto raise_error;
        }
        else {
            /* Normalize to raise <class>, <instance> */
            Py_DECREF(value);
            value = type;
            type = (PyObject*) ((PyInstanceObject*)type)->in_class;
            Py_INCREF(type);
        }
    }
    else {
        /* Not something you can raise.  You get an exception
           anyway, just not what you specified :-) */
        PyErr_Format(PyExc_TypeError,
                 "exceptions must be strings, classes, or "
                 "instances, not %s", type->ob_type->tp_name);
        goto raise_error;
    }
    PyErr_Restore(type, value, tb);
    return;
raise_error:
    Py_XDECREF(value);
    Py_XDECREF(type);
    Py_XDECREF(tb);
    return;
}

static int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type) {
    if (!type) {
        PyErr_Format(PyExc_SystemError, "Missing type object");
        return 0;
    }
    if (obj == Py_None || PyObject_TypeCheck(obj, type))
        return 1;
    PyErr_Format(PyExc_TypeError, "Cannot convert %s to %s",
        obj->ob_type->tp_name, type->tp_name);
    return 0;
}

static int __Pyx_InternStrings(__Pyx_InternTabEntry *t) {
    while (t->p) {
        *t->p = PyString_InternFromString(t->s);
        if (!*t->p)
            return -1;
        ++t;
    }
    return 0;
}

#include "compile.h"
#include "frameobject.h"
#include "traceback.h"

static void __Pyx_AddTraceback(char *funcname) {
    PyObject *py_srcfile = 0;
    PyObject *py_funcname = 0;
    PyObject *py_globals = 0;
    PyObject *empty_tuple = 0;
    PyObject *empty_string = 0;
    PyCodeObject *py_code = 0;
    PyFrameObject *py_frame = 0;
    
    py_srcfile = PyString_FromString(__pyx_filename);
    if (!py_srcfile) goto bad;
    py_funcname = PyString_FromString(funcname);
    if (!py_funcname) goto bad;
    py_globals = PyModule_GetDict(__pyx_m);
    if (!py_globals) goto bad;
    empty_tuple = PyTuple_New(0);
    if (!empty_tuple) goto bad;
    empty_string = PyString_FromString("");
    if (!empty_string) goto bad;
    py_code = PyCode_New(
        0,            /*int argcount,*/
        0,            /*int nlocals,*/
        0,            /*int stacksize,*/
        0,            /*int flags,*/
        empty_string, /*PyObject *code,*/
        empty_tuple,  /*PyObject *consts,*/
        empty_tuple,  /*PyObject *names,*/
        empty_tuple,  /*PyObject *varnames,*/
        empty_tuple,  /*PyObject *freevars,*/
        empty_tuple,  /*PyObject *cellvars,*/
        py_srcfile,   /*PyObject *filename,*/
        py_funcname,  /*PyObject *name,*/
        __pyx_lineno,   /*int firstlineno,*/
        empty_string  /*PyObject *lnotab*/
    );
    if (!py_code) goto bad;
    py_frame = PyFrame_New(
        PyThreadState_Get(), /*PyThreadState *tstate,*/
        py_code,             /*PyCodeObject *code,*/
        py_globals,          /*PyObject *globals,*/
        0                    /*PyObject *locals*/
    );
    if (!py_frame) goto bad;
    py_frame->f_lineno = __pyx_lineno;
    PyTraceBack_Here(py_frame);
bad:
    Py_XDECREF(py_srcfile);
    Py_XDECREF(py_funcname);
    Py_XDECREF(empty_tuple);
    Py_XDECREF(empty_string);
    Py_XDECREF(py_code);
    Py_XDECREF(py_frame);
}
