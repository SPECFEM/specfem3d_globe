#!/usr/bin/env python


from cig.addyndum.applications import Script, ParallelScript


class SpecfemScript(Script):
    
    
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
    
    #
    #--- configuration
    #
    
    def _configure(self):
        super(SpecfemScript, self)._configure()
        
        # validate absorbing conditions
        if self.solver.ABSORBING_CONDITIONS:
            NCHUNKS = self.mesher.NCHUNKS
            if NCHUNKS == 6:
                raise ValueError("cannot have absorbing conditions in the full Earth")
            elif NCHUNKS == 3:
                raise ValueError("absorbing conditions not supported for three chunks yet")
        
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
        super(SpecfemScript, self).collectUserInput(registry)
        return


class ParallelSpecfemScript(SpecfemScript, ParallelScript):

    #
    #--- inventory
    #
    
    import cig.addyndum.inventory as addyndum
    
    # default to LSF for Caltech cluster
    scheduler = addyndum.scheduler("scheduler", default="lsf")


    #
    #--- final initialization
    #
    
    def _init(self):
        super(ParallelSpecfemScript, self)._init()
        self.OUTPUT_FILES = self.outputDir
        self.LOCAL_PATH = self.scratchDir
        self.solver.LOCAL_PATH = self.LOCAL_PATH
        self.solver.OUTPUT_FILES = self.OUTPUT_FILES



class Specfem(ParallelSpecfemScript):
    
    #
    #--- inventory
    #
    
    import pyre.inventory as pyre
    
    command = pyre.str("command", validator=pyre.choice(['mesh', 'solve', 'go']), default="go")


    #
    #--- configuration
    #
    
    def _configure(self):
        super(Specfem, self)._configure()
        
        # declare the interpreter to be used on the compute nodes
        from os.path import join
        self.interpreter = join(self.outputDir, "pyspecfem3D") # includes solver


    #
    #--- executed on the login node in response to the user's command
    #
    
    def onLoginNode(self, context):
        
        import sys
        import shutil
        import journal
        
        # build the solver
        import __main__
        from os.path import dirname
        srcdir = dirname(__main__.__file__)
        self.solver.build(self, srcdir)
        
        # compute the total number of processors needed
        self.nodes = self.mesher.nproc()

        self.schedule(context) # bsub
        
        return


    #
    #--- executed when the batch job is scheduled
    #

    def onLauncherNode(self, context):
        self.launch(context) # mpirun
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


class MovieScript(ParallelSpecfemScript):

    import pyre.inventory as pyre
    
    formatName = pyre.str("format", default='OpenDX', validator=pyre.choice(['OpenDX', 'AVS-multi', 'AVS-single', 'GMT']))
    beginning  = pyre.int("beginning", default=1, validator=pyre.greaterEqual(1))
    end        = pyre.int("end", default=1, validator=pyre.greaterEqual(1))

    
    def _configure(self):
        super(MovieScript, self)._configure()

        if not self.solver.MOVIE_SURFACE:
            raise ValueError("movie frames were not saved by the solver")

        # count number of movie frames
        nframes = 0
        for it in xrange(self.beginning, self.end + 1):
            if it % self.solver.NTSTEP_BETWEEN_FRAMES == 0:
                nframes = nframes + 1
        if nframes == 0:
            raise ValueError('null number of frames')
        
        formatCode = { 'OpenDX':1, 'AVS-multi':2, 'AVS-single':3, 'GMT':4 }
        self.format = formatCode[self.formatName]


    def onLoginNode(self, context):
        self.nodes = 1
        self.schedule(context) # bsub

        
    def onLauncherNode(self, context):
        from PyxMeshfem import read_params_and_create_movie
        read_params_and_create_movie(self)


# entry points

def main(*args, **kwds):
    script = Specfem()
    script.run(*args, **kwds)


def create_movie_AVS_DX(*args, **kwds):
    script = MovieScript()
    script.run(*args, **kwds)


# end of file
