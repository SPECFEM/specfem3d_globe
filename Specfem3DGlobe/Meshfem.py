#!/usr/bin/env python


from MPIDriver import MPIDriver


class Meshfem(MPIDriver):
    
    def __init__(self):
        super(Meshfem, self).__init__("output_mesher.txt")

    def main(self, *args, **kwds):
        #self.code.meshfem3D(self)
        self.code.xxxxfem3D(self)


# end of file
