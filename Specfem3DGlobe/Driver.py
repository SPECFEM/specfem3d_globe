#!/usr/bin/env python


from Script import Script


class Driver(Script):

    """Base class for scripts which call into SPECFEM3D Fortran code."""
    
    def __init__(self):
        super(Driver, self).__init__()
        import Specfem3DGlobeCode as code
        self.code = code
    
    def readValue(self, name):
        """Callback from Fortran into Python."""
        l = name.split('.')
        o = self
        for n in l:
            try:
                o = getattr(o, n)
            except AttributeError:
                o = getattr(o.inventory, n)
        return o


# end of file
