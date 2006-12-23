#!/usr/bin/env python


#
# base class for all models
#


from pyre.components import Component


class Model(Component):
    
    # parameters common to all models
    import pyre.inventory as pyre
    from TopoBathy import TopoBathy
    ATTENUATION        = pyre.bool("attenuation")
    ELLIPTICITY        = pyre.bool("ellipticity")
    GRAVITY            = pyre.bool("gravity")
    OCEANS             = pyre.bool("oceans")
    ROTATION           = pyre.bool("rotation")
    TOPOGRAPHY         = pyre.bool("topography")
    topoBathy          = pyre.facility("topo-bathy", factory=TopoBathy)
    
    # configuration
    def _xconfigure(self):
        Component._configure(self)
        if not (self.TOPOGRAPHY or self.OCEANS):
            # We don't need topo-bathy data.
            self.jettison(Model.topoBathy)
        return


# end of file
