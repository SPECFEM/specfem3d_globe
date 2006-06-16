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
    headerFile                    = addyndum.outputFile("header-file",
                                                        default="${output-dir}/values_from_mesher.h")
    seismogramArchive             = addyndum.outputFile("seismogram-archive",
                                                        default="${output-dir}/seismograms.tar.gz")
    scratchSeismogramArchive      = addyndum.scratchFile("scratch-seismogram-archive",
                                                         default="${scratch-dir}/seismograms-${rank}.tar.gz")
    

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
        self.HEADER_FILE = abspath(self.headerFile.name) # always written by the mesher

        # filled-in by _init() in Specfem.py
        self.LOCAL_PATH = None
        self.OUTPUT_FILES = None


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
            print >> self.outputFile, "execute", specfem3D
        else:
            specfem3D(script) # call into Fortran
        return


    #
    #--- clean-up
    #
    
    def finiForComputeNode(self, context):
        """collect seismograms (part 1/2)"""

        if self.dry:
            return
        
        import os, shutil, tarfile
        from os.path import basename, join
        from glob import glob

        archive = self.scratchSeismogramArchive
        
        pattern = join(self.LOCAL_PATH, "*.sem*")
        files = []
        for name in glob(pattern):
            files.append(name)
        if len(files) == 0:
            # no seismograms on this node
            archive.close()
            os.remove(archive.name)
            return
        
        # A compressed tar file is between 10-15% of the size of the
        # raw data files.  By taring and compressing it now -- in
        # parallel, on local filesystems -- we sharply reduce the
        # amount of data we have to shovel over the network.
        
        tgz = tarfile.open(archive.name, "w:gz", archive)
        for name in files:
            arcname = basename(name)
            tgz.add(name, arcname)
        tgz.close()
        archive.close()

        # Copy the archive to the shared filesystem.

        outputDir = context.application.outputDir
        src = archive.name
        dst = join(outputDir, basename(src))
        shutil.copyfile(src, dst)

        return


    def finiForLauncherNode(self, context):
        """collect seismograms (part 2/2)"""
        
        if self.dry:
            return
        
        import os, tarfile
        from os.path import join
        from glob import glob

        archiveOut = self.seismogramArchive

        # Search for the intermediate seismogram archives delivered
        # from the compute nodes.
        
        pattern = join(self.OUTPUT_FILES, "seismograms-*.tar.gz")
        archivesIn = []
        for name in glob(pattern):
            archivesIn.append(name)
        if len(archivesIn) == 0:
            self._warning.log("No seismograms!")
            archiveOut.close()
            os.remove(archiveOut.name)
            return

        # Rearchive the seismograms in one, big archive.
        
        tgzOut = tarfile.open(archiveOut.name, "w:gz", archiveOut)
        for archiveIn in archivesIn:
            tgzIn = tarfile.open(archiveIn, "r:gz")
            for member in tgzIn.getmembers():
                seismogram = tgzIn.extractfile(member)
                tgzOut.addfile(member, seismogram)
            tgzIn.close()
        tgzOut.close()
        archiveOut.close()

        # Delete the intermediate seismogram archives.
        
        for archiveIn in archivesIn:
            os.remove(archiveIn)

        return


# end of file
