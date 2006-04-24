
#include <Python.h>
#include <stdio.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "config.h"

extern void initPyxParameters(void);
#ifdef USE_MPI
extern void initPyxMeshfem(void);
extern void initPyxSpecfem(void);
extern void initPyxMPI(void);
#endif

static int status;

struct _inittab inittab[] = {
    { "PyxParameters", initPyxParameters },
#ifdef USE_MPI
    { "PyxMeshfem", initPyxMeshfem },
    { "PyxSpecfem", initPyxSpecfem },
    { "PyxMPI", initPyxMPI },
#endif
    { 0, 0 }
};


#define FC_RUN_PYTHON_SCRIPT FC_FUNC_(run_python_script, RUN_PYTHON_SCRIPT)
void FC_RUN_PYTHON_SCRIPT()
{
    /* run the Python script */
#ifndef SCRIPT
#define SCRIPT Specfem
#endif
#define STR(s) #s
#define COMMAND(s) "from Specfem3DGlobe."STR(s)" import "STR(s)"; app = "STR(s)"(); app.run()"
    status = PyRun_SimpleString(COMMAND(SCRIPT)) != 0;
}


int main(int argc, char **argv)
{
#ifdef USE_MPI
    /*
      We are a worker on a compute node -- i.e., we are a process
      started by 'mpirun'.  Call MPI_Init() now to clean/set-up
      'argv'.
    */
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
	fprintf(stderr, "%s: MPI_Init failed! Exiting ...\n", argv[0]);
	return status;
    }
#else
    /*
      We are either the launcher, or the scheduler started directly by
      the user on the login node.  Don't call MPI_Init(), as MPICH-GM
      will die with SIGPIPE ("<MPICH-GM> Error: Need to obtain the job
      magic number in GMPI_MAGIC !").
    */
#endif
    
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }

    /* initialize Python */
    Py_Initialize();
    
    /* initialize sys.argv */
    PySys_SetArgv(argc, argv);
    
#define main 42
#if FC_MAIN == main
    /* run the Python script */
    FC_RUN_PYTHON_SCRIPT();
#else
    /* call the Fortran trampoline (which runs the Python script) */
    FC_MAIN();
#endif
    
    /* shut down Python */
    Py_Finalize();
    
#ifdef USE_MPI
    /* shut down MPI */
    MPI_Finalize();
#endif

    return status;
}


/* end of file */
