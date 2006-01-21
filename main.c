
#include <Python.h>
#include <stdio.h>
#include "config.h"


#define STR(s) #s
#define RUN_SCRIPT(s) "from Specfem3DGlobe."STR(s)" import "STR(s)"; app = "STR(s)"(); app.run()"
#ifndef SCRIPT
#define SCRIPT Specfem
#endif

extern void initSpecfem3DGlobeCode(void);
extern void initPyxMPI(void);

static int status;

struct _inittab inittab[] = {
    { "Specfem3DGlobeCode", initSpecfem3DGlobeCode },
    { "PyxMPI", initPyxMPI },
    { 0, 0 }
};


int main(int argc, char **argv)
{
    /* add our extension module */
    if (PyImport_ExtendInittab(inittab) == -1) {
        fprintf(stderr, "%s: PyImport_ExtendInittab failed! Exiting ...", argv[0]);
        return 1;
    }
    
    /* initialize Python */
    Py_Initialize();
    
    /* initialize sys.argv */
    PySys_SetArgv(argc, argv);
    
    /* call the Fortran trampoline */
    FC_MAIN();
    
    /* shut down Python */
    Py_Finalize();
    
    return status;
}


void FC_FUNC_(run_python_script, RUN_PYTHON_SCRIPT)()
{
    /* run the Python command */
    status = PyRun_SimpleString(RUN_SCRIPT(SCRIPT)) != 0;
}


/* end of file */
