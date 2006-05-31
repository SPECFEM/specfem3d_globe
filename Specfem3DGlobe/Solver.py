#!/usr/bin/env python


from cig.addyndum.components import Component
from pyre.units.time import minute


class Solver(Component):

    
    # name by which the user refers to this component
    componentName = "solver"
    
    
    #
    #--- parameters
    #
    
    import pyre.inventory as pyre
    import cig.addyndum.inventory as addyndum
    from CMTSolution import cmtValidator

    cmtSolution                   = addyndum.inputFile("cmt-solution",
                                                       default="DATA/CMTSOLUTION")
    stations                      = addyndum.inputFile("stations",
                                                       default="DATA/STATIONS")

    outputFile                    = addyndum.outputFile("output-file",
                                                        default="${output-dir}/output_solver.txt")

    ABSORBING_CONDITIONS          = pyre.bool("absorbing-conditions")
    MOVIE_SURFACE                 = pyre.bool("movie-surface")
    MOVIE_VOLUME                  = pyre.bool("movie-volume")
    RECEIVERS_CAN_BE_BURIED       = pyre.bool("receivers-can-be-buried")
    PRINT_SOURCE_TIME_FUNCTION    = pyre.bool("print-source-time-function")
    dry                           = pyre.bool("dry")

    HDUR_MOVIE                    = pyre.float("hdur-movie")
    record_length                 = pyre.dimensional("record-length", default=0.0*minute)

    NTSTEP_BETWEEN_FRAMES         = pyre.int("ntstep-between-frames")
    NTSTEP_BETWEEN_OUTPUT_INFO    = pyre.int("ntstep-between-output-info")
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = pyre.int("ntstep-between-output-seismos")
    NUMBER_OF_RUNS                = pyre.int("number-of-runs")
    NUMBER_OF_THIS_RUN            = pyre.int("number-of-this-run")

    SAVE_FORWARD                  = pyre.bool("save-forward")
    simulation_type               = pyre.str("simulation-type", validator=pyre.choice(['forward', 'adjoint', 'both']), default='forward')
    
    
    #
    #--- configuration
    #
    
    def _configure(self):
        Component._configure(self)
        
        # convert to minutes
        self.RECORD_LENGTH_IN_MINUTES = self.record_length / minute
        
        # convert to the ID numbers understood by the Fortran code
        st = { 'forward':1, 'adjoint':2, 'both':3 }
        self.SIMULATION_TYPE = st[self.simulation_type]
        
        return


    #
    #--- final initialization
    #
    
    def _init(self):
        Component._init(self)

        from os.path import abspath
        self.CMTSOLUTION = abspath(self.cmtSolution.name)
        self.STATIONS = abspath(self.stations.name)


    #
    #--- building
    #
    
    def build(self, script, srcdir):
        """Build the solver."""

        import os, os.path, sys
        from os.path import abspath

        outputDir = abspath(script.outputDir)
        pyspecfem3D = abspath(script.interpreter)

        wd = os.getcwd()
        print "cd", srcdir
        os.chdir(srcdir)
        
        # create the include file for the solver
        self.createheader(script)
        
        # now finally build the solver
        argv = ['make', 'OUTPUT_DIR=' + outputDir, pyspecfem3D]
        print ' '.join(argv)
        status = os.spawnvp(os.P_WAIT, argv[0], argv)
        if status != 0:
            sys.exit("%s: %s: exit %d" % (sys.argv[0], argv[0], status))

        print "cd", wd
        os.chdir(wd)

        return
    
    
    def createheader(self, script):
        """Create the include file for the solver."""

        import os, shutil, sys
        
        # This path is hardwired into the Fortran source.
        oldHeader = 'OUTPUT_FILES/values_from_mesher.h'
        
        # First generate the header into a temporary file.
        from tempfile import mktemp
        newHeader  = mktemp()
        self.HEADER_FILE = newHeader
        from PyxParameters import create_header_file
        create_header_file(script) # call into Fortran
        
        # Did the header file change?
        argv = ['diff', oldHeader, newHeader]
        print ' '.join(argv)
        status = os.spawnvp(os.P_WAIT, argv[0], argv)
        if status == 0:
            # Nope!  Nothing to do here.
            os.remove(newHeader)
            return
        if status != 1:
            # diff countered a problem
            os.remove(newHeader)
            sys.exit("%s: %s: exit %d" % (sys.argv[0], argv[0], status))

        # Replace the old header with the new one.
        print "mv", newHeader, oldHeader
        shutil.move(newHeader, oldHeader)

        return
    
    
    #
    #--- execution
    #
    
    def execute(self, script):
        """Execute the solver."""
        from PyxSpecfem import specfem3D
        if self.dry:
            print >> outputFile, "execute", specfem3D
        else:
            specfem3D(script) # call into Fortran
        return


# end of file
