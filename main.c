
#include <Python.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "config.h"


extern void initPyxParameters(void);
extern void initPyxMeshfem(void);
#ifdef BUILDING_SOLVER
extern void initPyxSpecfem(void);
#endif

static int g_status;

#define COMMAND \
"import sys; " \
"path = sys.argv[1]; " \
"requires = sys.argv[2]; " \
"entry = sys.argv[3]; " \
"path = path.split(':'); " \
"path.extend(sys.path); " \
"sys.path = path; " \
"from merlin import loadObject; " \
"entry = loadObject(entry); " \
"entry(sys.argv[3:], kwds={'requires': requires})"

/* include the implementation of _mpi */
#include "mpi/_mpi.c"

struct _inittab inittab[] = {
    { "PyxParameters", initPyxParameters },
    { "PyxMeshfem", initPyxMeshfem },
    { "_mpi", init_mpi },
#ifdef BUILDING_SOLVER
    { "PyxSpecfem", initPyxSpecfem },
#endif
    { 0, 0 }
};


#define FC_PY_MAIN FC_FUNC_(fc_py_main, FC_PY_MAIN)
void FC_PY_MAIN()
{
    /* run the Python command */
    g_status = PyRun_SimpleString(COMMAND) != 0;
}


int main(int argc, char **argv)
{
#ifdef USE_MPI
    /* initialize MPI */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "%s: MPI_Init failed! Exiting ...", argv[0]);
        return 1;
    }
#endif

    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }

    if (argc < 3 || strcmp(argv[1], "--pyre-start") != 0) {
        return Py_Main(argc, argv);
    }

    /* make sure 'sys.executable' is set to the path of this program  */
    Py_SetProgramName(argv[0]);

    /* initialize Python */
    Py_Initialize();

    /* initialize sys.argv */
    PySys_SetArgv(argc - 1, argv + 1);

#define main 42
#if FC_MAIN == main
    /* start Python */
    FC_PY_MAIN();
#else
    /* call the Fortran trampoline (which, in turn, starts Python) */
    FC_MAIN();
#endif

    /* shut down Python */
    Py_Finalize();

#ifdef USE_MPI
    /* shut down MPI */
    MPI_Finalize();
#endif

    return g_status;
}


/* end of file */
