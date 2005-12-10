#!/usr/bin/env python


from pyre.applications.Script import Script as PyreScript
import sys


def myexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    import pdb
    pdb.post_mortem(traceback)


class Script(PyreScript):

    class Inventory(PyreScript.Inventory):

        from pyre.inventory import bool, facility, str
        from Mesher import Mesher
        from Model import ModelFacility
        from Solver import Solver
        
        ABSORBING_CONDITIONS          = bool("absorbing-conditions")

        LOCAL_PATH                    = str("local-path")
        
        mesher                        = facility("mesher", factory=Mesher, args=["mesher"])
        model                         = ModelFacility("model")
        solver                        = facility("solver", factory=Solver, args=["solver"])
    
    def _init(self):
        PyreScript._init(self)
        self.MODEL = self.inventory.model.className

    def run(self, *args, **kwds):
        for i in xrange(0, len(sys.argv)):
            if sys.argv[i] == '-pdb':
                if sys.stdin.isatty():
                    sys.excepthook = myexcepthook
                del sys.argv[i]
                break
        PyreScript.run(self, *args, **kwds)

    def readValue(self, name):
        l = name.split('.')
        o = self
        for n in l:
            try:
                o = getattr(o, n)
            except AttributeError:
                o = getattr(o.inventory, n)
        return o


# end of file
