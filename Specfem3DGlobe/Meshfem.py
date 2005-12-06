#!/usr/bin/env python


from Script import Script
import Specfem3DGlobeCode


class Meshfem(Script):
    
    def __init__(self):
        Script.__init__(self, "meshfem")

    def main(self, *args, **kwds):
        Specfem3DGlobeCode.meshfem3D(self)


# end of file
