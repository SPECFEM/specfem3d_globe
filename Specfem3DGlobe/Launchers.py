#!/usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MPI Launchers:
#     (Replacement) Launcher for MPICH
#     Launcher for LAM/MPI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from mpi.Launcher import Launcher
from pyre.inventory.Facility import Facility
import os, sys


defaultEnvironment = "[EXPORT_ROOT,LD_LIBRARY_PATH,PYTHONPATH]"


class LauncherNone(Launcher):
    class Inventory(Launcher.Inventory):
        from pyre.inventory import bool, str
        dry = bool("dry", default=False)
    def __init__(self):
        Launcher.__init__(self, "none")
    def launch(self, script, method):
        if self.inventory.dry:
            print '%s.%s()' % (script, method)
            return
        script.perform(method)


class LauncherMPI(Launcher):


    class Inventory(Launcher.Inventory):

        from pyre.inventory import bool, list, str

        dry = bool("dry", default=False)
        dry.meta['tip'] = "prints the command line and exits"
        
        debug = bool("debug", default=False)

        Launcher.Inventory.nodes.meta['tip'] = """number of machine nodes"""
        Launcher.Inventory.nodelist.meta['tip'] = """a comma-separated list of machine names in square brackets (e.g., [101-103,105,107])"""
        nodegen = str("nodegen")
        nodegen.meta['tip'] = """a printf-style format string, used in conjunction with 'nodelist' to generate the list of machine names (e.g., "n%03d")"""
        
        extra = str("extra")
        extra.meta['tip'] = "extra arguments to pass to mpirun"
        
        command = str("command", default="mpirun")
        executable = str("executable", default=sys.argv[0])
        arguments = list("arguments", default=sys.argv[1:])


    def launch(self, script, method):
        args = self._buildArgumentList(method)
        if not args:
            return self.inventory.dry
        
        command = " ".join(args)
        self._info.log("executing: {%s}" % command)
        
        if self.inventory.dry:
            print command
            return True
        
        os.system(command)
        return True

            
    def _buildArgumentList(self, method):
        if not self.nodes:
            self.nodes = len(self.nodelist)

        if self.nodes < 1:
            self.nodes = 1
        self.nodes = 2

        # build the command
        args = []
        args.append(self.inventory.command)
        self._appendMpiRunArgs(args)

        args.append(self.inventory.executable)
        args += self.inventory.arguments
        #args.append("--mode=worker")
        args.append("--command=" + method)

        return args

    
    def _appendMpiRunArgs(self, args):
        args.append(self.inventory.extra)
        args.append("-np %d" % self.nodes)
        
        # use only the specific nodes specified explicitly
        if self.nodelist:
            self._appendNodeListArgs(args)


class LauncherMPICH(LauncherMPI):


    class Inventory(LauncherMPI.Inventory):
 
        from pyre.inventory import str

        machinefile = str("machinefile", default="mpirun.nodes")
        machinefile.meta['tip'] = """filename of machine file"""


    def __init__(self):
        LauncherMPI.__init__(self, "mpich")


    def _appendNodeListArgs(self, args):
        machinefile = self.inventory.machinefile
        nodegen = self.inventory.nodegen
        file = open(machinefile, "w")
        for node in self.nodelist:
            file.write((nodegen + '\n') % node)
        file.close()
        args.append("-machinefile %s" % machinefile)


class LauncherLAMMPI(LauncherMPI):


    class Inventory(LauncherMPI.Inventory):

        from pyre.inventory import list

        environment = list("environment", default=defaultEnvironment)
        environment.meta['tip'] = """a comma-separated list of environment variables to export to the batch job"""


    def __init__(self):
        LauncherMPI.__init__(self, "lam-mpi")


    def _appendMpiRunArgs(self, args):
        args.append("-x %s" % ','.join(self.inventory.environment))
        super(LauncherLAMMPI, self)._appendMpiRunArgs(args)


    def _appendNodeListArgs(self, args):
        nodegen = self.inventory.nodegen
        args.append("n" + ",".join([(nodegen) % node for node in self.nodelist]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LauncherFacility
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


builtInLauncherClasses = {
    'none':    LauncherNone,
    'lam-mpi': LauncherLAMMPI,
    'mpich':   LauncherMPICH,
    }

class LauncherFacility(Facility):

    def __init__(self, name, default=None):
        if default is None:
            default = "none"
        Facility.__init__(self, name, default=default)
        return
    
    def _retrieveComponent(self, instance, componentName):
        cls = builtInLauncherClasses.get(componentName, None)
        if cls is None:
            return Facility._retrieveComponent(self, instance, componentName)
        launcher = cls()
        launcher.aliases.append(self.name)
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')
        return launcher, locator

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# main
if __name__ == "__main__":
    from pyre.applications.Script import Script
    class TestApp(Script):
        class Inventory(Script.Inventory):
            launcher = LauncherFacility("launcher", default="mpich")
        def main(self, *args, **kwds):
            launcher = self.inventory.launcher
            if launcher:
                print ' '.join(launcher._buildArgumentList())
            return
        def __init__(self):
            Script.__init__(self, "TestApp")
    app = TestApp()
    app.run()


# end of file 
