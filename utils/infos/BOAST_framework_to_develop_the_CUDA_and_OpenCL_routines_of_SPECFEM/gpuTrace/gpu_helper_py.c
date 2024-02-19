#include <Python.h>

#include "ocl_helper.h"
#include "cuda_helper.h"

#define PYTHON_MOD_NAME "parse_gpu_program"

#ifndef PYTHON_MOD_PATH
#define PYTHON_MOD_PATH /home/kevin/travail/sample/gpuTrace
#endif

#define STR(_str_) #_str_
#define xSTR(_x_, _str_) _x_ ## _str_
#define LSTR_(_str_) xSTR(L, #_str_)
#define LSTR(_str_) LSTR_(_str_)

static PyObject *ocl_pf_parse, *ocl_pf_prep_progr;
static PyObject *cuda_pf_get_lookup_table;

static void init_helper(void);
void cuda_init_helper(void) { init_helper(); }
void ocl_init_helper(void) { init_helper(); }

static void init_helper(void) {
  static int initialized = 0;

  if (initialized) {
    return;
  } else {
    initialized = 1;
  }

  PyObject *pDict, *pModule;

  Py_Initialize();

  PyObject *sys_path;
  PyObject *path;

  sys_path = PySys_GetObject("path");
  path = PyUnicode_FromWideChar(LSTR(PYTHON_MOD_PATH), -1);
  PyList_Append(sys_path, path);

  pModule = PyImport_ImportModule(PYTHON_MOD_NAME);
  PyErr_Print();

  pDict = PyModule_GetDict(pModule);

  ocl_pf_parse = PyDict_GetItemString(pDict, "ocl_parse");
  ocl_pf_prep_progr = PyDict_GetItemString(pDict, "ocl_prepare_program");

  cuda_pf_get_lookup_table = PyDict_GetItemString(pDict, "cuda_get_raw_lookup_table");

  Py_DECREF(pDict);
  Py_DECREF(pModule);
}

void ocl_handle_program(void *program,
                    unsigned int count,
                    const char **strings,
                    const size_t *lengths) {
  PyObject *p_progr_uid = PyUnicode_FromFormat("%p", program);
  PyObject *p_program_lines = PyTuple_New(count);
  PyObject *p_params;
  int i;

  for (i = 0; i < count; i++) {
    PyObject *line;

    if (!lengths || !lengths[i]) {
      line = PyUnicode_FromString(strings[i]);
    } else {
      line = PyUnicode_FromStringAndSize(strings[i], lengths[i]);
    }
    PyTuple_SetItem(p_program_lines, i, line);
  }

  p_params = PyTuple_Pack(2, p_progr_uid, p_program_lines);

  PyObject_Call(ocl_pf_prep_progr, p_params, NULL);
  PyErr_Print();


  Py_DECREF(p_params);
  Py_DECREF(p_program_lines);
  Py_DECREF(p_progr_uid);
}

char **ocl_handle_create_kernel(void *program, void *kernel, const char *name) {
  PyObject *p_progr_uid = PyUnicode_FromFormat("%p", program);
  PyObject *p_kern_uid = PyUnicode_FromFormat("%p", kernel);
  PyObject *p_kern_name = PyUnicode_FromString(name);
  PyObject *p_result, *p_params = PyTuple_Pack(2, p_progr_uid, p_kern_name);
  Py_ssize_t py_result_size;
  int result_size;
  char **param_types_names;
  int i;

  p_result = PyObject_Call(ocl_pf_parse, p_params, NULL);
  PyErr_Print();

  if (!p_result) {
    printf("ERROR, no result ...\n");
    perror("help...");
    exit(-1);
  }
  py_result_size = PyList_Size(p_result);
  result_size = Py_SAFE_DOWNCAST(py_result_size, Py_ssize_t, int);

  param_types_names = malloc(sizeof(char *)*(result_size + 1));

  for (i = 0; i < result_size; i++) {
    PyObject *p_ascii_str = PyUnicode_AsASCIIString(PyList_GetItem(p_result, i));

    param_types_names[i] = PyBytes_AsString(p_ascii_str);
  }
  param_types_names[result_size] = NULL;

  Py_DECREF(p_result);
  Py_DECREF(p_params);
  Py_DECREF(p_progr_uid);
  Py_DECREF(p_kern_uid);
  Py_DECREF(p_kern_name);

  return param_types_names;
}

struct kernel_lookup_s *cuda_get_lookup_table(void) {
  PyObject *p_params = PyTuple_New(0);
  PyObject *p_result = PyObject_Call(cuda_pf_get_lookup_table, p_params, NULL);
  Py_ssize_t py_result_size;
  PyObject *p_result_items;
  struct kernel_lookup_s *lookup_table;

  int i, j;
  int nb_kernels;

  PyErr_Print();
  if (!p_result) {
    printf("ERROR, no result ...\n");
    perror("help...");
    exit(-1);
  }

  py_result_size = PyDict_Size(p_result);
  nb_kernels = Py_SAFE_DOWNCAST(py_result_size, Py_ssize_t, int);

  lookup_table = malloc(sizeof(struct kernel_lookup_s) * (nb_kernels + 1));

  p_result_items = PyDict_Items(p_result);
  lookup_table[nb_kernels].address = NULL;
  for (i = 0; i < nb_kernels; i++) {
    PyObject *py_addr_info = PyList_GetItem(p_result_items, i);
    PyObject *py_addr = PyTuple_GetItem(py_addr_info, 0);
    PyObject *py_info_tpl = PyTuple_GetItem(py_addr_info, 1);
    PyObject *py_name = PyUnicode_AsASCIIString(PyTuple_GetItem(py_info_tpl, 0));
    PyObject *py_params_lst = PyTuple_GetItem(py_info_tpl, 1);
    Py_ssize_t py_nb_params = PyList_Size(py_params_lst);

    lookup_table[i].address = PyLong_AsVoidPtr(py_addr);
    lookup_table[i].name = PyBytes_AsString(py_name);
    lookup_table[i].nb_params = Py_SAFE_DOWNCAST(py_nb_params, Py_ssize_t, int);

    lookup_table[i].params = malloc(sizeof(struct param_info_s) * lookup_table[i].nb_params);
    for (j = 0; j < lookup_table[i].nb_params; j++) {
      PyObject *py_type_name = PyList_GetItem(py_params_lst, j);
      PyObject *py_param_type = PyUnicode_AsASCIIString(PyTuple_GetItem(py_type_name, 0));
      PyObject *py_param_name = PyUnicode_AsASCIIString(PyTuple_GetItem(py_type_name, 1));

      lookup_table[i].params[j].name = PyBytes_AsString(py_param_name);
      lookup_table[i].params[j].type = PyBytes_AsString(py_param_type);
    }
  }

  Py_DECREF(p_result);
  Py_DECREF(p_result_items);
  Py_DECREF(p_params);

  return lookup_table;
}
