#!/usr/bin/env python


from cig.addyndum.components import Component
from pyre.units.angle import deg


class Mesher(Component):

    
    # name by which the user refers to this component
    componentName = "mesher"

    
    #
    #--- parameters
    #
    
    import pyre.inventory as pyre
    import cig.addyndum.inventory as addyndum

    outputFile                    = addyndum.outputFile("output-file",
                                                        default="${output-dir}/output_mesher.txt")
    
    SAVE_MESH_FILES               = pyre.bool("save-files")
    dry                           = pyre.bool("dry")

    angular_width_eta             = pyre.dimensional("angular-width-eta", default=90.0*deg)
    angular_width_xi              = pyre.dimensional("angular-width-xi", default=90.0*deg)
    center_latitude               = pyre.dimensional("center-latitude", default=0.0*deg)
    center_longitude              = pyre.dimensional("center-longitude", default=0.0*deg)
    gamma_rotation_azimuth        = pyre.dimensional("gamma-rotation-azimuth", default=0.0*deg)

    NCHUNKS                       = pyre.int("nchunks", validator=pyre.choice([1,2,3,6]), default=6)
    NEX_ETA                       = pyre.int("nex-eta", default=64)
    NEX_XI                        = pyre.int("nex-xi", default=64)
    NPROC_ETA                     = pyre.int("nproc-eta", validator=pyre.greaterEqual(1), default=1)
    NPROC_XI                      = pyre.int("nproc-xi", validator=pyre.greaterEqual(1), default=1)
    
    
    #
    #--- configuration
    #
    
    def _configure(self):
        Component._configure(self)

        # convert to degrees
        self.ANGULAR_WIDTH_ETA_IN_DEGREES = self.angular_width_eta / deg
        self.ANGULAR_WIDTH_XI_IN_DEGREES  = self.angular_width_xi / deg
        self.CENTER_LATITUDE_IN_DEGREES   = self.center_latitude / deg
        self.CENTER_LONGITUDE_IN_DEGREES  = self.center_longitude / deg
        self.GAMMA_ROTATION_AZIMUTH       = self.gamma_rotation_azimuth / deg

        # this MUST be 90 degrees for two chunks or more to match geometrically
        if (self.NCHUNKS > 1 and
            self.ANGULAR_WIDTH_XI_IN_DEGREES != 90.0):
            raise ValueError("'angular-width-xi' must be 90 degrees for more than one chunk")
        
        # this can be any value in the case of two chunks
        if (self.NCHUNKS > 2 and
            self.ANGULAR_WIDTH_ETA_IN_DEGREES != 90.0):
            raise ValueError("'angular-width-eta' must be 90 degrees for more than two chunks")
        
        # check that topology is correct if more than two chunks
        if (self.NCHUNKS > 2 and
            self.NPROC_XI != self.NPROC_ETA):
            raise ValueError("'nproc-xi' and 'nproc-eta' must be equal for more than two chunks")

        # check that size can be coarsened in depth twice (block size multiple of 8)
        if ((self.NEX_XI / 8) % self.NPROC_XI) != 0:
            raise ValueError("'nex-xi' must be a multiple of 8*nproc-xi")
        if ((self.NEX_ETA / 8) % self.NPROC_ETA) != 0:
            raise ValueError("'nex-eta' must be a multiple of 8*nproc-eta")
        if self.NEX_XI % 8 != 0:
            raise ValueError("'nex-xi' must be a multiple of 8")
        if self.NEX_ETA % 8 != 0:
            raise ValueError("'nex-eta' must be a multiple of 8")

        # check that sphere can be cut into slices without getting negative Jacobian
        if self.NEX_XI < 48:
            raise ValueError("'nex-xi' must be greater than 48 to cut the sphere into slices with positive Jacobian")
        if self.NEX_ETA < 48:
            raise ValueError("'nex-eta' must be greater than 48 to cut the sphere into slices with positive Jacobian")

        # number of elements in each slice (i.e. per processor)
        NEX_PER_PROC_XI = self.NEX_XI / self.NPROC_XI
        NEX_PER_PROC_ETA = self.NEX_ETA / self.NPROC_ETA
        
        # one doubling layer in outer core (block size multiple of 16)
        if NEX_PER_PROC_XI % 16 != 0:
            raise ValueError("'nex-xi / nproc-xi' (i.e., elements per processor) must be a multiple of 16 for outer core doubling")
        if NEX_PER_PROC_ETA % 16 != 0:
            raise ValueError("'nex-eta / nproc-eta' (i.e., elements per processor) must be a multiple of 16 for outer core doubling")

        # check that number of elements per processor is the same in both directions
        if self.NCHUNKS > 2 and NEX_PER_PROC_XI != NEX_PER_PROC_ETA:
            raise ValueError("must have the same number of elements per processor in both directions for more than two chunks")
        
        return

    
    def nproc(self):
        """Return the total number of processors needed."""
        return self.NCHUNKS * self.NPROC_XI * self.NPROC_ETA

    
    #
    #--- execution
    #
    
    def execute(self, script):
        """Execute the mesher."""
        from PyxMeshfem import meshfem3D
        if self.dry:
            print >> self.outputFile, "execute", meshfem3D
        else:
            meshfem3D(script) # call into Fortran
        return


# end of file
