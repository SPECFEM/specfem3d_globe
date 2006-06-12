
#include <Python.h>

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

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


void handler()
{
    pid_t pid;
    char *argv[10];
    char pidstr[42];
    int status;
    struct sigaction sa;
    int tfd;
    char gdbScriptName[] = "gdbXXXXXX";
    static char gdbScript[] = "bt\nkill\n";
    
#if 0
    sa.sa_handler = SIG_DFL;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO;
    sigaction(SIGABRT, &sa, NULL);
#endif
    
    tfd = mkstemp(gdbScriptName);
    if (tfd == -1) {
        abort();
    }
    write(tfd, gdbScript, sizeof(gdbScript));
    close(tfd);
    
    pid = getpid();
    sprintf(pidstr, "%d", pid);
    
    argv[0] = "/usr/bin/gdb";
    argv[1] = "-batch";
    argv[2] = "-x";
    argv[3] = gdbScriptName;
    argv[4] = g_argv[0];
    argv[5] = pidstr;
    argv[6] = 0;
    
    int cpid = fork();
    if (cpid) {
        /* parent */
        wait(&status);
        exit(0);
    } else {
        /* child */
        execv(argv[0], argv);
    }
}


void trapSignals()
{
    struct sigaction sa;
    
    sa.sa_sigaction = handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO | SA_ONESHOT;
    sigaction(SIGSEGV, &sa, NULL);
    /*sigaction(SIGABRT, &sa, NULL);*/
    
    return;
}


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
    
    trapSignals();
    
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
