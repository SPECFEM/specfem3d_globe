#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.units.angle import deg
import Specfem3DGlobeCode


class Mesher(Component):
    
    class Inventory(Component.Inventory):

        from pyre.inventory import bool, choice, dimensional, greaterEqual, int
    
        SAVE_MESH_FILES               = bool("save-files")
        
        angular_width_eta             = dimensional("angular-width-eta", default=90.0*deg)
        angular_width_xi              = dimensional("angular-width-xi", default=90.0*deg)
        center_latitude               = dimensional("center-latitude", default=0.0*deg)
        center_longitude              = dimensional("center-longitude", default=0.0*deg)
        gamma_rotation_azimuth        = dimensional("gamma-rotation-azimuth", default=0.0*deg)
        
        NCHUNKS                       = int("nchunks", validator=choice([1,2,3,6]), default=6)
        NEX_ETA                       = int("nex-eta", default=64)
        NEX_XI                        = int("nex-xi", default=64)
        NPROC_ETA                     = int("nproc-eta", validator=greaterEqual(1), default=1)
        NPROC_XI                      = int("nproc-xi", validator=greaterEqual(1), default=1)
        
    def __init__(self, name):
        Component.__init__(self, name, "mesher")

    def _init(self):
        Component._init(self)

        # convert to degrees
        self.ANGULAR_WIDTH_ETA_IN_DEGREES = self.inventory.angular_width_eta / deg
        self.ANGULAR_WIDTH_XI_IN_DEGREES  = self.inventory.angular_width_xi / deg
        self.CENTER_LATITUDE_IN_DEGREES   = self.inventory.center_latitude / deg
        self.CENTER_LONGITUDE_IN_DEGREES  = self.inventory.center_longitude / deg
        self.GAMMA_ROTATION_AZIMUTH       = self.inventory.gamma_rotation_azimuth / deg

        # copy to locals for convenience
        NCHUNKS   = self.inventory.NCHUNKS
        NEX_ETA   = self.inventory.NEX_ETA
        NEX_XI    = self.inventory.NEX_XI
        NPROC_ETA = self.inventory.NPROC_ETA
        NPROC_XI  = self.inventory.NPROC_XI

        # this MUST be 90 degrees for two chunks or more to match geometrically
        if (NCHUNKS > 1 and
            self.ANGULAR_WIDTH_XI_IN_DEGREES != 90.0):
            raise ValueError("'angular-width-xi' must be 90 degrees for more than one chunk")
        
        # this can be any value in the case of two chunks
        if (NCHUNKS > 2 and
            self.ANGULAR_WIDTH_ETA_IN_DEGREES != 90.0):
            raise ValueError("'angular-width-eta' must be 90 degrees for more than two chunks")
        
        # check that topology is correct if more than two chunks
        if (NCHUNKS > 2 and
            NPROC_XI != NPROC_ETA):
            raise ValueError("'nproc-xi' and 'nproc-eta' must be equal for more than two chunks")

        # check that size can be coarsened in depth twice (block size multiple of 8)
        if ((NEX_XI / 8) % NPROC_XI) != 0:
            raise ValueError("'nex-xi' must be a multiple of 8*nproc-xi")
        if ((NEX_ETA / 8) % NPROC_ETA) != 0:
            raise ValueError("'nex-eta' must be a multiple of 8*nproc-eta")
        if NEX_XI % 8 != 0:
            raise ValueError("'nex-xi' must be a multiple of 8")
        if NEX_ETA % 8 != 0:
            raise ValueError("'nex-eta' must be a multiple of 8")

        # check that sphere can be cut into slices without getting negative Jacobian
        if NEX_XI < 48:
            raise ValueError("'nex-xi' must be greater than 48 to cut the sphere into slices with positive Jacobian")
        if NEX_ETA < 48:
            raise ValueError("'nex-eta' must be greater than 48 to cut the sphere into slices with positive Jacobian")

        # number of elements in each slice (i.e. per processor)
        NEX_PER_PROC_XI = NEX_XI / NPROC_XI
        NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA
        
        # one doubling layer in outer core (block size multiple of 16)
        if NEX_PER_PROC_XI % 16 != 0:
            raise ValueError("'nex-xi / nproc-xi' (i.e., elements per processor) must be a multiple of 16 for outer core doubling")
        if NEX_PER_PROC_ETA % 16 != 0:
            raise ValueError("'nex-eta / nproc-eta' (i.e., elements per processor) must be a multiple of 16 for outer core doubling")

        # check that number of elements per processor is the same in both directions
        if NCHUNKS > 2 and NEX_PER_PROC_XI != NEX_PER_PROC_ETA:
            raise ValueError("must have the same number of elements per processor in both directions for more than two chunks")
        
        return

    def nproc(self):
        """Return the total number of processors needed."""
        return (self.inventory.NCHUNKS *
                self.inventory.NPROC_XI *
                self.inventory.NPROC_ETA)


# end of file
