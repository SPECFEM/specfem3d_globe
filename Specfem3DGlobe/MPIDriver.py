#!/usr/bin/env python


from Driver import Driver
from mpi.Application import Application as MPIApplication


class MPIDriver(Driver, MPIApplication):

    """Base class for scripts which run SPECFEM3D Fortran code in
    parallel under MPI."""

    class Inventory(Driver.Inventory, MPIApplication.Inventory): pass

    def __init__(self, outputFilename):
        super(MPIDriver, self).__init__()
        self.outputFilename = outputFilename
    
    def _init(self):
        super(MPIDriver, self)._init()
        # make sure the output directory is writable
        if self.inventory.mode == "server": # NYI: MPI_Comm_rank() == 0???
            self.checkOutputDir()
        return
    
    def checkOutputDir(self):
        from os import remove
        from os.path import join
        outputDir = self.inventory.OUTPUT_FILES
        temp = join(outputDir, self.outputFilename)
        try:
            f = open(temp, 'w')
        except IOError:
            self.raisePropertyValueError('output-dir')
        f.close()
        remove(temp)



# end of file
