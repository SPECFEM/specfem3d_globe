#!/usr/bin/env python


from cig.addyndum.components import Component


class TopoBathy(Component):

    
    # name by which the user refers to this component
    componentName = "topo-bathy"

    
    # parameters
    import pyre.inventory as pyre
    import cig.addyndum.inventory as addyndum
    
    NX_BATHY              = pyre.int("nx", default=5400)
    NY_BATHY              = pyre.int("ny", default=2700)
    RESOLUTION_TOPO_FILE  = pyre.int("resolution", default=4); RESOLUTION_TOPO_FILE.meta['tip'] = "resolution of topography file in minutes"
    dataFile              = addyndum.inputFile("data-file",
                                               default="DATA/topo_bathy/topo_bathy_etopo4_smoothed_window_7.dat")
    
    # final initialization
    def _init(self):
        Component._init(self)
        self.PATHNAME_TOPO_FILE = self.dataFile.name


# end of file
