#!/usr/bin/env python


from cig.addyndum.applications import ParallelScript


class Specfem(ParallelScript):
    
    
    # name by which the user refers to this application
    componentName = "Specfem3DGlobe"
    
    
    #
    #--- inventory
    #
    
    import pyre.inventory as pyre
    import cig.addyndum.inventory as addyndum
    
    from Mesher import Mesher
    from Model import model
    from Solver import Solver

    outputDir                     = addyndum.outputDir("output-dir", default="OUTPUT_FILES")
        
    model                         = model("model", default="isotropic_prem")
    mesher                        = addyndum.facility("mesher", factory=Mesher)
    solver                        = addyndum.facility("solver", factory=Solver)
    
    command                       = pyre.str("command", validator=pyre.choice(['mesh', 'solve', 'go']), default="go")

    # default to LSF for Caltech cluster
    scheduler                     = addyndum.scheduler("scheduler", default="lsf")

    
    #
    #--- configuration
    #
    
    def _configure(self):
        ParallelScript._configure(self)
        
        # declare the interpreter to be used on the compute nodes
        from os.path import join
        self.interpreter = join(self.outputDir, "pyspecfem3D") # includes solver

        # validate absorbing conditions
        if self.solver.ABSORBING_CONDITIONS:
            NCHUNKS = self.mesher.NCHUNKS
            if NCHUNKS == 6:
                raise ValueError("cannot have absorbing conditions in the full Earth")
            elif NCHUNKS == 3:
                raise ValueError("absorbing conditions not supported for three chunks yet")
        
        return


    #
    #--- final initialization
    #
    
    def _init(self):
        ParallelScript._init(self)
        self.OUTPUT_FILES = self.outputDir
        self.LOCAL_PATH = self.scratchDir
        self.solver.LOCAL_PATH = self.LOCAL_PATH
        self.solver.OUTPUT_FILES = self.OUTPUT_FILES


    #
    #--- executed on the login node in response to the user's command
    #
    
    def onLoginNode(self, context):
        
        import sys
        import shutil
        import journal
        
        pyspecfem3D = sys.executable # the current executable
        
        # build the solver
        self.solver.build(self, context.kwds['srcdir'])
        
        # compute the total number of processors needed
        self.nodes = self.mesher.nproc()

        self.schedule(context) # bsub
        
        return


    #
    #--- executed when the batch job is scheduled
    #

    def onLauncherNode(self, context):
        self.launch(context) # mpirun
        # tar up results
        return

    
    #
    #--- executed in parallel on the compute nodes
    #
    
    def onComputeNodes(self, context):

        # execute mesher and/or solver
        
        if self.command == "mesh" or self.command == "go":
            self.mesher.execute(self)
        
        if self.command == "solve" or self.command == "go":
            self.solver.execute(self)

        return

    
    #
    #--- support for reading Par_file
    #
    

    def collectUserInput(self, registry):
        # read Par_files given on the command line
        from ParFileCodec import ParFileCodec
        from os.path import isfile, splitext
        argv = self.argv
        self.argv = []
        codec = ParFileCodec(self.name)
        for arg in argv:
            base, ext = splitext(arg)
            if not ext and isfile(arg):
                shelf = codec.open(base)
                paramRegistry = shelf['inventory'].getFacility(self.name)
                if paramRegistry:
                    self.updateConfiguration(paramRegistry)
            else:
                self.argv.append(arg)
        super(Specfem, self).collectUserInput(registry)
        return



def main(*args, **kwds):
    script = Specfem()
    script.run(*args, **kwds)


# end of file
