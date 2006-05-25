
#include <Python.h>
#include <stdio.h>
#include "config.h"

extern void initPyxParameters(void);
#ifdef USE_MPI
extern void initPyxMeshfem(void);
extern void initPyxSpecfem(void);
extern void initPyxMPI(void);
#endif

static int g_status;
int g_argc;
char **g_argv;

struct _inittab inittab[] = {
    { "PyxParameters", initPyxParameters },
#ifdef USE_MPI
    { "PyxMeshfem", initPyxMeshfem },
    { "PyxSpecfem", initPyxSpecfem },
    { "PyxMPI", initPyxMPI },
#endif
    { 0, 0 }
};


#define FC_PY_MAIN FC_FUNC_(fc_py_main, FC_PY_MAIN)
void FC_PY_MAIN()
{
    /* start Python */
    g_status = Py_Main(g_argc, g_argv);
}


int main(int argc, char **argv)
{
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }
    
    g_argc = argc;
    g_argv = argv;
#define main 42
#if FC_MAIN == main
    /* start Python */
    FC_PY_MAIN();
#else
    /* call the Fortran trampoline (which, in turn, starts Python) */
    FC_MAIN();
#endif
    
    return g_status;
}


/* end of file */
