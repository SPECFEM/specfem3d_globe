#!/usr/bin/env python


from CodecConfig import CodecConfig
from mpi.Application import Application
import sys


# patch Pyre to our liking
import PyrePatches


def myexcepthook(type, value, traceback):
    sys.__excepthook__(type, value, traceback)
    import pdb
    pdb.post_mortem(traceback)


class Script(Application):

    class Inventory(Application.Inventory):

        from pyre.inventory import bool, facility, str
        from Launchers import LauncherFacility
        from Mesher import Mesher
        from Model import ModelFacility
        from Solver import Solver
        
        ABSORBING_CONDITIONS          = bool("absorbing-conditions")

        LOCAL_PATH                    = str("local-path")
        
        launcher                      = LauncherFacility("launcher")
        mesher                        = facility("mesher", factory=Mesher, args=["mesher"])
        model                         = ModelFacility("model")
        solver                        = facility("solver", factory=Solver, args=["solver"])
    
    def __init__(self):
        Application.__init__(self, "Specfem3DGlobe")

    def _init(self):
        Application._init(self)
        self.MODEL = self.inventory.model.className
        # compute the total number of processors needed
        nodes = self.inventory.mesher.nproc()
        self.inventory.launcher.inventory.nodes = nodes
        self.inventory.launcher.nodes = nodes
        
    def run(self, *args, **kwds):
        for i in xrange(0, len(sys.argv)):
            if sys.argv[i] == '-pdb':
                if sys.stdin.isatty():
                    sys.excepthook = myexcepthook
                del sys.argv[i]
                break
        try:
            Application.run(self, *args, **kwds)
        except PyrePatches.PropertyValueError, error:
            self.reportPropertyValueError(error)

    def reportPropertyValueError(self, error):
        import journal
        import linecache
        j = journal.journal()
        j.device.renderer = PyrePatches.PropertyValueError.Renderer()
        property = error.args[0]
        locator = error.args[1]
        entry = journal.diagnostics.Entry.Entry()
        meta = entry.meta
        meta["name"] = property.name
        meta["error"] = error
        if locator:
            meta["filename"] = locator.source
            try:
                line = locator.line
            except AttributeError:
                meta["src"] = ""
                meta["line"] = "???"
            else:
                src = linecache.getline(locator.source, locator.line)
                src = src.rstrip()
                meta["src"] = src
                meta["line"] = locator.line
        else:
            meta["filename"] = "???"
            meta["line"] = "???"
            meta["src"] = ""
        j.record(entry)
        
    def initializeCurator(self, curator, registry):
        cfg = CodecConfig()
        curator.registerCodecs(cfg)
        return Application.initializeCurator(self, curator, registry)
        
    def collectUserInput(self, registry):
        curator = self.getCurator()
        configRegistry = curator.getTraits(self.name, extraDepositories=[], encoding='cfg')
        registry.update(configRegistry)
        
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
