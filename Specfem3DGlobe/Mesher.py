#!/usr/bin/env python


from pyre.components.Component import Component
from pyre.units.angle import deg
import Specfem3DGlobeCode


class Mesher(Component):
    
    class Inventory(Component.Inventory):

        from pyre.inventory import bool, choice, dimensional, int
    
        SAVE_MESH_FILES               = bool("save-files")
        
        angular_width_eta             = dimensional("angular-width-eta", default=0.0*deg)
        angular_width_xi              = dimensional("angular-width-xi", default=0.0*deg)
        center_latitude               = dimensional("center-latitude", default=0.0*deg)
        center_longitude              = dimensional("center-longitude", default=0.0*deg)
        gamma_rotation_azimuth        = dimensional("gamma-rotation-azimuth", default=0.0*deg)
        
        NCHUNKS                       = int("nchunks", validator=choice([1,2,3,6]))
        NEX_ETA                       = int("nex-eta")
        NEX_XI                        = int("nex-xi")
        NPROC_ETA                     = int("nproc-eta")
        NPROC_XI                      = int("nproc-xi")
        
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

        # an example of validation
        if (self.inventory.NCHUNKS >= 3 and
            self.inventory.NPROC_XI != self.inventory.NPROC_ETA):
            raise ValueError("nproc-xi and nproc-eta must be equal for more than two chunks")

    def nproc(self):
        """Return the total number of processors needed."""
        return (self.inventory.NCHUNKS *
                self.inventory.NPROC_XI *
                self.inventory.NPROC_ETA)


# end of file
