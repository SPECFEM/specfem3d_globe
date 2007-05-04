#!/usr/bin/env python


from pyre.applications import Script
from mpi import Application as MPIApplication


class SpecfemScript(Script):
    
    
    # name by which the user refers to this application
    name = "Specfem3DGlobe"
    
    
    #
    #--- inventory
    #
    
    import pyre.inventory as pyre
    
    from Mesher import Mesher
    from Solver import Solver

    outputDir                     = pyre.str("output-dir", default="OUTPUT_FILES")
    scratchDir                    = pyre.str("scratch-dir", default="/scratch")
        
    model                         = pyre.facility("model", default="isotropic_prem")
    mesher                        = pyre.facility("mesher", factory=Mesher)
    solver                        = pyre.facility("solver", factory=Solver)

    
    #
    #--- validation
    #
    
    def _validate(self, context):
        super(SpecfemScript, self)._validate(context)
        
        # validate absorbing conditions
        if self.solver.ABSORBING_CONDITIONS:
            NCHUNKS = self.mesher.NCHUNKS
            if NCHUNKS == 6:
                context.error(ValueError("cannot have absorbing conditions in the full Earth"))
            elif NCHUNKS == 3:
                context.error(ValueError("absorbing conditions not supported for three chunks yet"))

        return


    #
    #--- initialization
    #
    
    def _init(self):

        from os import makedirs
        from os.path import abspath, isdir
        from pyre.util import expandMacros

        macros = { 'job.name': self.job.task }
        self.outputDir = expandMacros(self.outputDir, macros)
        self.mpiExecutable = expandMacros(self.mpiExecutable, macros)
        
        if not isdir(self.outputDir):
            makedirs(self.outputDir)
        
        return


    #
    #--- support for reading Par_file
    #

    def collectUserInput(self, registry, context):
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
        super(SpecfemScript, self).collectUserInput(registry, context)
        return


    #
    #--- .odb files
    #

    def _getPrivateDepositoryLocations(self):
        from os.path import dirname, isdir, join
        models = join(dirname(__file__), 'models')
        assert isdir(models)
        return [models]


class ParallelSpecfemScript(SpecfemScript, MPIApplication):

    #
    #--- final initialization
    #
    
    def _init(self):
        super(ParallelSpecfemScript, self)._init()
        self.OUTPUT_FILES = self.outputDir
        self.LOCAL_PATH = self.scratchDir
        self.solver.setOutputDirectories(
            LOCAL_PATH = self.LOCAL_PATH,
            OUTPUT_FILES = self.OUTPUT_FILES
            )


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
        self.mpiExecutable = join(self.outputDir, "mpipyspecfem3D") # includes solver

        # compute the total number of processors needed
        self.nodes = self.mesher.nproc()
        
        return


    #
    #--- executed on the login node in response to the user's command
    #
    
    def onLoginNode(self, *args, **kwds):
        
        import sys
        import shutil
        import journal
        
        # build the solver
        import __main__
        from os.path import dirname
        srcdir = dirname(__main__.__file__)
        self.solver.build(self, srcdir)
        
        # schedule the job (bsub)
        self.scheduleJob(*args, **kwds)
        
        return


    #
    #--- executed when the batch job is scheduled
    #

    def onLauncherNode(self, *args, **kwds):
        self.launchParallelRun(*args, **kwds) # mpirun
        self.solver.processOutputFiles()
        return

    
    #
    #--- executed in parallel on the compute nodes
    #
    
    def onComputeNodes(self, *args, **kwds):
        import os

        # There is a race condition here if the scratch directory is
        # on a shared filesystem.  For now, incorporating ${rank} is
        # necessary.

        makedirs(self.LOCAL_PATH)

        # execute mesher and/or solver
        
        if self.command == "mesh" or self.command == "go":
            self.mesher.execute(self)
        
        if self.command == "solve" or self.command == "go":
            self.solver.execute(self)

        # Once again, there is a race condition here if the scratch
        # directory is on a shared filesystem.  We could avoid the
        # OSError exceptions on the individual files, at least, if we
        # knew which nodes were reponsible for which files...

        from os.path import join

        try:
            files = os.listdir(self.LOCAL_PATH)
            for file in files:
                try:
                    os.remove(join(self.LOCAL_PATH, file))
                except OSError:
                    pass
        except OSError:
            pass

        # Don't try to remove the intermediate-level directories that
        # may have been created by makedirs() -- they might be in use
        # by another job anyway.
        try:
            os.rmdir(self.LOCAL_PATH)
        except OSError:
            pass

        return


class MovieScript(ParallelSpecfemScript):

    import pyre.inventory as pyre
    
    formatName = pyre.str("format", default='OpenDX', validator=pyre.choice(['OpenDX', 'AVS-multi', 'AVS-single', 'GMT']))
    beginning  = pyre.int("beginning", default=1, validator=pyre.greaterEqual(1))
    end        = pyre.int("end", default=1, validator=pyre.greaterEqual(1))


    def _validate(self, context):
        super(MovieScript, self)._validate(context)

        if not self.solver.MOVIE_SURFACE:
            context.error(ValueError("movie frames were not saved by the solver"))

        # count number of movie frames
        nframes = 0
        for it in xrange(self.beginning, self.end + 1):
            if it % self.solver.NTSTEP_BETWEEN_FRAMES == 0:
                nframes = nframes + 1
        if nframes == 0:
            context.error(ValueError('null number of frames'),
                          items=[self.metainventory.beginning,
                                 self.metainventory.end])
        
        return

    
    def _configure(self):
        super(MovieScript, self)._configure()

        formatCode = { 'OpenDX':1, 'AVS-multi':2, 'AVS-single':3, 'GMT':4 }
        self.format = formatCode[self.formatName]

        self.nodes = 1

        return


    def onLauncherNode(self, *args, **kwds):
        from PyxMeshfem import read_params_and_create_movie
        read_params_and_create_movie(self)



# this shall be moved to the framework

def makedirs(name, mode=0777):
    """Like os.makedirs(), but multi-{thread,process} safe."""
    import os.path as path
    from os import curdir, mkdir
    head, tail = path.split(name)
    if not tail:
        head, tail = path.split(head)
    if head and tail and not path.exists(head):
        makedirs(head, mode)
        if tail == curdir:           # xxx/newdir/. exists if xxx/newdir exists
            return
    try:
        mkdir(name, mode)
    except OSError, error:
        # 17 == "File exists"
        if hasattr(error, 'errno') and error.errno == 17:
            pass
        else:
            raise
    return



# entry points

def main(*args, **kwds):
    script = Specfem()
    script.run(*args, **kwds)


def create_movie_AVS_DX(*args, **kwds):
    script = MovieScript()
    script.run(*args, **kwds)


# end of file
