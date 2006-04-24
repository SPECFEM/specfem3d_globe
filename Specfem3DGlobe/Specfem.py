#!/usr/bin/env python


from Script import Script
import os, sys


class Specfem(Script):
    
    #
    #--- inventory
    #
    
    class Inventory(Script.Inventory):

        from pyre.inventory import facility, str
        from Mesher import Mesher
        from Model import ModelFacility
        from Properties import OutputDir
        from Solver import Solver

        LOCAL_PATH                    = str("local-path")
        OUTPUT_FILES                  = OutputDir("output-dir", default="OUTPUT_FILES")
        
        mesher                        = facility("mesher", factory=Mesher, args=["mesher"])
        model                         = ModelFacility("model")
        solver                        = facility("solver", factory=Solver, args=["solver"])

    #
    #--- methods performed on the login node
    #

    def mesh(self):
        """Schedule a launch of the mesher."""
        self.setupOutputDir()
        self.schedule('launchMesher') # bsub
        
    def solve(self):
        """Build, then schedule a launch of the solver."""
        self.setupOutputDir()
        self.build()
        self.schedule('launchSolver') # bsub
        
    def go(self):
        """Build, then schedule a launch of the mesher and the solver."""
        self.setupOutputDir()
        self.build()
        self.schedule('launchMesherSolver') # bsub
    
    def setupOutputDir(self):
        
        from os import makedirs
        from os.path import isdir, join
        outputDir = self.inventory.OUTPUT_FILES
        
        # the launcher is simply a copy of the current executable
        xspecfem3D = sys.argv[0]
        specfem3Dlauncher = self.launcherExe
        self.shellCommand('cp', xspecfem3D, specfem3Dlauncher)

        # copy the Python scripts
        import Specfem3DGlobe
        srcScriptDir = Specfem3DGlobe.__path__[0]
        destScriptDir = join(outputDir, "Specfem3DGlobe")
        if not isdir(destScriptDir):
            makedirs(destScriptDir)
        for srcScript in os.listdir(srcScriptDir):
            destScript = join(destScriptDir, srcScript)
            srcScript = join(srcScriptDir, srcScript)
            if not isdir(srcScript):
                self.shellCommand('cp', srcScript, destScript)
        return

    def build(self):
        """Build the solver."""

        outputDir = self.inventory.OUTPUT_FILES
        
        # create the include file for the solver
        self.createheader()

        # now finally build the solver
        specfem3D = self.computeExe
        self.shellCommand('make', 'SIMULATION_DIR=' + outputDir, specfem3D)
        return

    def createheader(self):
        """Create the include file for the solver."""
        from PyxParameters import create_header_file
        create_header_file(self)

    def hello(self):
        self.setupOutputDir()
        self.build()
        self.schedule('launchHello') # bsub
        
    #
    #--- methods performed when the batch job runs
    #
    
    def launchMesher(self):
        """Launch the mesher."""
        self.checkOutputDir("output_mesher.txt")
        self.launch('excecuteMesher') # mpirun
    
    def launchSolver(self):
        """Launch the solver."""
        self.checkOutputDir("output_solver.txt")
        self.launch('executeSolver') # mpirun
    
    def launchMesherSolver(self):
        """Launch the mesher and the solver."""
        self.checkOutputDir("output_mesher.txt", "output_solver.txt")
        self.launch('executeMesherSolver') # mpirun
    
    def launchHello(self):
        self.launch('executeHello') # mpirun
    
    def checkOutputDir(self, *outputFilenames):
        from os import remove
        from os.path import join
        outputDir = self.inventory.OUTPUT_FILES
        for outputFilename in outputFilenames:
            if outputFilename is None:
                continue
            temp = join(outputDir, outputFilename)
            try:
                f = open(temp, 'w')
            except IOError:
                self.raisePropertyValueError('output-dir')
            f.close()
            remove(temp)
        return
    
    #
    #--- methods performed in parallel on the compute nodes
    #
    
    def excecuteMesher(self):
        """Execute the mesher.  Only performed on compute nodes."""
        from PyxMeshfem import meshfem3D
        meshfem3D(self)

    def executeSolver(self):
        """Execute the solver.  Only performed on compute nodes."""
        from PyxSpecfem import specfem3D
        specfem3D(self)

    def executeMesherSolver(self):
        """Execute the mesher and the solver.  Only performed on compute nodes."""
        from PyxMeshfem import meshfem3D
        from PyxSpecfem import specfem3D
        meshfem3D(self)
        specfem3D(self)

    def executeHello(self):
        import socket
        hostname = socket.gethostname()
        rank = 0 #MPI_Comm_rank()
        size = 0 #MPI_Comm_size()
        print "[%03d/%03d] Hello from '%s'!" % (rank, size, hostname)
        print "[%03d/%03d] Good-bye from '%s'!" % (rank, size, hostname)

    #
    #--- high-level Pyre Component initialization
    #

    def _init(self):
        
        super(Specfem, self)._init()

        # declare the exectuables used for launching and computing
        from os.path import join
        outputDir = self.inventory.OUTPUT_FILES
        self.launcherExe = join(outputDir, "specfem3Dlauncher")
        self.computeExe = join(outputDir, "specfem3D")

        # get the model identifier
        self.MODEL = self.inventory.model.className
        
        # compute the total number of processors needed
        nodes = self.inventory.mesher.nproc()
        self.inventory.launcher.inventory.nodes = nodes
        self.inventory.launcher.nodes = nodes
        self.inventory.scheduler.inventory.nodes = nodes
        self.inventory.scheduler.nodes = nodes
        
        # validate absorbing conditions
        if self.inventory.solver.inventory.ABSORBING_CONDITIONS:
            NCHUNKS = self.inventory.mesher.inventory.NCHUNKS
            if NCHUNKS == 6:
                raise ValueError("cannot have absorbing conditions in the full Earth")
            elif NCHUNKS == 3:
                raise ValueError("absorbing conditions not supported for three chunks yet")
        return

    #
    #--- low-level object initialization
    #

    def __init__(self):
        super(Specfem, self).__init__("Specfem3DGlobe")


# main
if __name__ == "__main__":
    script = Specfem()
    script.run()


# end of file
