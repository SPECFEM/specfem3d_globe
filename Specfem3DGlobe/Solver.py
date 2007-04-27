#!/usr/bin/env python


from pyre.components import Component
from pyre.units.time import minute


class Solver(Component):

    
    # name by which the user refers to this component
    name = "solver"
    
    
    #
    #--- parameters
    #
    
    import pyre.inventory as pyre
    from CMTSolution import cmtValidator

    cmtSolution                   = pyre.inputFile("cmt-solution",
                                                   default="DATA/CMTSOLUTION",
                                                   validator=cmtValidator)
    stations                      = pyre.inputFile("stations",
                                                   default="DATA/STATIONS")

    ABSORBING_CONDITIONS          = pyre.bool("absorbing-conditions")
    MOVIE_SURFACE                 = pyre.bool("movie-surface")
    MOVIE_VOLUME                  = pyre.bool("movie-volume")
    RECEIVERS_CAN_BE_BURIED       = pyre.bool("receivers-can-be-buried")
    PRINT_SOURCE_TIME_FUNCTION    = pyre.bool("print-source-time-function")
    dry                           = pyre.bool("dry")

    HDUR_MOVIE                    = pyre.float("hdur-movie")
    record_length                 = pyre.dimensional("record-length", default=0.0*minute)

    NTSTEP_BETWEEN_FRAMES         = pyre.int("ntstep-between-frames", default=100, validator=pyre.greaterEqual(1))
    NTSTEP_BETWEEN_OUTPUT_INFO    = pyre.int("ntstep-between-output-info", default=200, validator=pyre.greaterEqual(1))
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = pyre.int("ntstep-between-output-seismos", default=5000000, validator=pyre.greaterEqual(1))
    NTSTEP_BETWEEN_READ_ADJSRC    = pyre.int("ntstep-between-read-adjsrc", default=1000, validator=pyre.greaterEqual(1))
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

        from os.path import abspath, join
        self.CMTSOLUTION = abspath(self.cmtSolution.name)
        self.STATIONS = abspath(self.stations.name)


    def setOutputDirectories(self, LOCAL_PATH, OUTPUT_FILES):
        from os.path import abspath, join
        
        self.LOCAL_PATH = LOCAL_PATH
        self.OUTPUT_FILES = OUTPUT_FILES

        # always written by the mesher
        self.HEADER_FILE = abspath(join(OUTPUT_FILES, 'values_from_mesher.h'))


    #
    #--- building
    #
    
    def build(self, script, srcdir):
        """Build the solver."""

        import os, os.path, sys
        from os.path import abspath, join

        outputDir = abspath(script.outputDir)
        pyspecfem3D = abspath(script.mpiExecutable)

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
        from os.path import exists
        from PyxParameters import create_header_file
        
        # This path is hardwired into the Fortran source.
        oldHeader = 'OUTPUT_FILES/values_from_mesher.h'

        # If the header doesn't exist, simply create it.
        if not exists(oldHeader):
            self.HEADER_FILE = oldHeader
            create_header_file(script) # call into Fortran
            return
        
        # First generate the header into a temporary file.
        from tempfile import mktemp
        newHeader  = mktemp()
        self.HEADER_FILE = newHeader
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
            print "execute", specfem3D
        else:
            specfem3D(script) # call into Fortran
        return


    #
    #--- clean-up
    #
    
    def processOutputFiles(self):
        """process output files"""
        
        if self.dry:
            return
        
        import os, tarfile
        from os.path import basename, join
        from cig.seismo.sac import asc2sac

        seismogramArchive = join(self.OUTPUT_FILES, "seismograms.tar.gz")
        archiveOut = open(seismogramArchive, "w")
        skipList = ['mpipyspecfem3D', basename(archiveOut.name)]

        # Archive output files, deleting seismograms after archiving.

        filesIn = []
        seismogramsIn = []
        for name in os.listdir(self.OUTPUT_FILES):
            if name in skipList:
                continue
            pathname = join(self.OUTPUT_FILES, name)
            if name.endswith(".sem") or name.endswith(".semd"):
                seismogramsIn.append((pathname, name))
                name_sac = name + ".sac"
                pathname_sac = join(self.OUTPUT_FILES, name_sac)
                asc2sac(pathname, pathname_sac)
                seismogramsIn.append((pathname_sac, name_sac))
            else:
                filesIn.append((pathname, name))
        if len(filesIn) == 0:
            self._warning.log("No output files!")
            archiveOut.close()
            os.remove(archiveOut.name)
            return

        tgzOut = tarfile.open(archiveOut.name, "w:gz", archiveOut)
        
        # Archive seismograms.

        for name, arcname in seismogramsIn:
            tgzOut.add(name, arcname)
        
        # Archive other output files.
        
        for name, arcname in filesIn:
            tgzOut.add(name, arcname)
        
        tgzOut.close()
        archiveOut.close()

        # Delete seismograms.
        
        for name, arcname in seismogramsIn:
            os.remove(name)

        return


# end of file
