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
    from CMTSolution import cmtValidator

    cmtSolution                   = pyre.inputFile("cmt-solution", default="DATA/CMTSOLUTION", validator=cmtValidator)
    stations                      = pyre.inputFile("stations", default="DATA/STATIONS")

    outputFile                    = pyre.outputFile("output-file", default="output_solver.txt")

    ABSORBING_CONDITIONS          = pyre.bool("absorbing-conditions")
    MOVIE_SURFACE                 = pyre.bool("movie-surface")
    MOVIE_VOLUME                  = pyre.bool("movie-volume")
    RECEIVERS_CAN_BE_BURIED       = pyre.bool("receivers-can-be-buried")
    PRINT_SOURCE_TIME_FUNCTION    = pyre.bool("print-source-time-function")

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
    #--- building
    #
    
    def build(self, script):
        """Build the solver."""

        import os
        
        # create the include file for the solver
        self.createheader(script)
        
        # now finally build the solver
        pyspecfem3D = script.computeExe
        argv = ['make', 'OUTPUT_DIR=' + script.outputDir, pyspecfem3D]
        print ' '.join(argv)
        status = os.spawnvp(os.P_WAIT, argv[0], argv)
        if status != 0:
            sys.exit("%s: %s: exit %d" % (sys.argv[0], args[0], status))
        
        return
    
    
    def createheader(self, script):
        """Create the include file for the solver."""
        self.CMTSOLUTION = self.cmtSolution.name
        self.STATIONS = self.stations.name
        from PyxParameters import create_header_file
        create_header_file(script) # call into Fortran
    
    
    #
    #--- execution
    #
    
    def execute(self, script):
        """Execute the solver."""
        self.CMTSOLUTION = self.cmtSolution.name
        self.STATIONS = self.stations.name
        from PyxSpecfem import specfem3D
        #specfem3D(script) # call into Fortran
        print "execute", specfem3D


# end of file
