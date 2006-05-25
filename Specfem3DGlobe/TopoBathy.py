#!/usr/bin/env python


from cig.addyndum.components import Component


class TopoBathy(Component):

    # name by which the user refers to this component
    componentName = "topo-bathy"

    # parameters
    import pyre.inventory as pyre
    NX_BATHY              = pyre.int("nx", default=5400)
    NY_BATHY              = pyre.int("ny", default=2700)
    RESOLUTION_TOPO_FILE  = pyre.int("resolution", default=4); RESOLUTION_TOPO_FILE.meta['tip'] = "resolution of topography file in minutes"
    dataFile              = pyre.inputFile("data-file",
                                            default="DATA/topo_bathy/topo_bathy_etopo4_smoothed_window_7.dat")

    # configuration
    def _configure(self):
        Component._configure(self)
        #self.PATHNAME_TOPO_FILE = self.dataFile.name


# end of file
