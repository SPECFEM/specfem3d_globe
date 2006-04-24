#!/usr/bin/env python


from CodecConfig import CodecConfig
import journal
from pyre.applications.Script import Script as PyreScript
import os, sys


# patch Pyre to our liking
import PyrePatches


class Script(PyreScript):

    class Inventory(PyreScript.Inventory):
        
        from pyre.inventory import str
        from Launchers import LauncherFacility
        from Schedulers import SchedulerFacility
        
        launcher                      = LauncherFacility("launcher", default="mpich")
        scheduler                     = SchedulerFacility("scheduler", default="lsf")
        
        command                       = str("command", default="go")
        
    class ShellError(RuntimeError): pass

    #
    #--- top-level entry points for all cases (login, batch, compute)
    #
    
    def run(self, *args, **kwds):
        # This method is run *before* the Pyre framework has processed
        # our command-line arguments and input files.
        try:
            super(Script, self).run(*args, **kwds)
        except PyrePatches.PropertyValueError, error:
            self.reportPropertyValueError(error)
        return

    def main(self, *args, **kwds):
        # This method is run *after* the Pyre framework has processed
        # our command-line arguments and input files.  All the
        # properties and components in our inventory are ready.
        command = self.inventory.command
        try:
            method = getattr(self, command)
        except AttributeError:
            print "%s: unknown command: %s" % (sys.argv[0], command)
            sys.exit(1)
        else:
            try:
                method()
            except Script.ShellError, error:
                print "%s: %s" % (sys.argv[0], error)
                sys.exit(1)
        return

    #
    #--- methods performed on the login node
    #
    
    def schedule(self, method):
        scheduler = self.inventory.scheduler
        scheduler.inventory.executable = self.launcherExe
        scheduler.schedule(self, method)

    #
    #--- methods performed when the batch job runs
    #
    
    def launch(self, method):
        launcher = self.inventory.launcher
        launcher.inventory.executable = self.computeExe
        launcher.launch(self, method)

    #
    #--- performing
    #

    def perform(self, method):
        method = getattr(self, method)
        method()

    #
    #--- shell
    #

    def shellCommand(self, *args):
        print ' '.join(args)
        status = os.spawnvp(os.P_WAIT, args[0], args)
        if status != 0:
            raise Script.ShellError("%s: exit %d" % (args[0], status))
        return
    
    #
    #--- error handling
    #

    def raisePropertyValueError(self, traitName):
        desc = self.inventory.getTraitDescriptor(traitName)
        trait = self.inventory.getTrait(traitName)
        raise PyrePatches.PropertyValueError, (trait, desc.locator)
    
    def reportPropertyValueError(self, error):
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

    #
    #--- support for '.cfg' files and specifying parameter files on the command line
    #
    
    def initializeCurator(self, curator, registry):
        cfg = CodecConfig()
        curator.registerCodecs(cfg)
        return super(Script, self).initializeCurator(curator, registry)
        
    def collectUserInput(self, registry):
        # read INI-style .cfg files
        curator = self.getCurator()
        configRegistry = curator.getTraits(self.name, extraDepositories=[], encoding='cfg')
        self.updateConfiguration(configRegistry)
        # read parameter files given on the command line
        from os.path import isfile, splitext
        for arg in self.argv:
            if isfile(arg):
                base, ext = splitext(arg)
                encoding = ext[1:] # NYI: not quite
                codec = self.getCurator().codecs.get(encoding)
                if codec:
                    shelf = codec.open(base)
                    paramRegistry = shelf['inventory'].getFacility(self.name)
                    if paramRegistry:
                        self.updateConfiguration(paramRegistry)
                else:
                    self._error.log("unknown encoding: %s" % ext)
            else:
                self._error.log("cannot open '%s'" % arg)
        return

    #
    #--- support for reading Python variables from Fortran
    #
    
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
    
    #
    #--- high-level Pyre Component initialization
    #

    def _init(self):
        super(Script, self)._init()
        self.launcherExe = None
        self.computeExe = None
        
    #
    #--- low-level object initialization
    #

    def __init__(self, name):
        super(Script, self).__init__(name)
        self._error = journal.error(self.name)


# end of file
