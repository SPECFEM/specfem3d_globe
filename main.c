
#include <Python.h>
#include <mpi.h>
#include <stdio.h>
#include "config.h"


#define STR(s) #s
#define RUN_SCRIPT(s) "from Specfem3DGlobe."STR(s)" import "STR(s)"; app = "STR(s)"(); app.run()"
#ifndef SCRIPT
#define SCRIPT Specfem
#endif

extern void initSpecfem3DGlobeCode(void);

struct _inittab inittab[] = {
    { "Specfem3DGlobeCode", initSpecfem3DGlobeCode },
    { 0, 0 }
};


int main(int argc, char **argv)
{
    int status;
    FILE *fp;
    
    /* initialize MPI */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "%s: MPI_Init failed! Exiting ...", argv[0]);
        return 1;
    }
    
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }
    
    /* initialize Python */
    Py_Initialize();
    
    /* initialize sys.argv */
    PySys_SetArgv(argc, argv);
    
    /* run the Python command */
    status = PyRun_SimpleString(RUN_SCRIPT(SCRIPT)) != 0;
    
    /* shut down Python and MPI */
    Py_Finalize();
    MPI_Finalize();
    
    return status;
}


/* end of file */
