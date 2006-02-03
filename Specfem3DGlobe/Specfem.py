#!/usr/bin/env python


from MPIDriver import MPIDriver


class Specfem(MPIDriver):
    
    def __init__(self):
        super(Specfem, self).__init__("output_solver.txt")

    def main(self, *args, **kwds):
        self.code.specfem3D(self)


# end of file
