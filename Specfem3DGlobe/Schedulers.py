#!/usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are Schedulers for batch schedulers found on the TeraGrid.
# These should be incorporated into Pythia eventually.

# With something like StringTemplate by Terence Parr and Marq Kole,
# the batch scripts could be generated entirely from an
# inventory-data-driven template.

#     http://www.stringtemplate.org/doc/python-doc.html

# This code uses a hybrid approach, mixing Python logic with primitive
# templates powered by pyre.util.expandMacros().
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from pyre.components.Component import Component
from pyre.inventory.Facility import Facility
from pyre.util import expandMacros
import os, sys


def hms(t):
    return (int(t / 3600), int((t % 3600) / 60), int(t % 60))

defaultEnvironment = "[EXPORT_ROOT,LD_LIBRARY_PATH,PYTHONPATH]"


class Scheduler(Component):
    class Inventory(Component.Inventory):
        from pyre.inventory import bool, int, list, str
        dry = bool("dry", default=False)
        executable = str("executable", default=sys.argv[0])
        arguments = list("arguments", default=sys.argv[1:])
        nodes = int("nodes", default=0)
    def __init__(self, name):
        Component.__init__(self, name, "scheduler")
    def schedule(self, script, method):
        raise NotImplementedError


class SchedulerNone(Scheduler):
    def __init__(self):
        Scheduler.__init__(self, "none")
    def schedule(self, script, method):
        if self.inventory.dry:
            print '%s.%s()' % (script, method)
            return
        script.perform(method)


# for debugging
class KernelScheduler(Scheduler):
    def __init__(self):
        Scheduler.__init__(self, "kernel")
    def schedule(self, script, method):
        args = [self.inventory.executable] + self.inventory.arguments
        args.append("--command=" + method)
        if self.inventory.dry:
            print ' '.join(args)
            return
        status = os.spawnvp(os.P_WAIT, args[0], args)
        if status != 0:
            print "%s: %s: exit %d" % (sys.argv[0], args[0], status)
            os.exit(1)
        return


class BatchScheduler(Scheduler):


    class Inventory(Scheduler.Inventory):

        from pyre.inventory import bool, dimensional, list, str
        from pyre.units.time import minute

        task = str("task")

        # Ignore 'nodegen' so that the examples will work without modification.
        nodegen = str("nodegen")
        nodegen.meta['tip'] = """(ignored)"""

        walltime = dimensional("walltime", default=0*minute)
        mail = bool("mail", default=False)
        queue = str("queue")

        directory = str("directory", default="${cwd}")
        script = str("script", default="${directory}/script")
        stdout = str("stdout", default="${directory}/stdout")
        stderr = str("stderr", default="${directory}/stderr")
        environment = list("environment", default=defaultEnvironment)

        cwd = str("cwd", default=os.getcwd())


    def schedule(self, script, method):

        self.inventory.executable = os.path.abspath(self.inventory.executable)

        if self.inventory.dry:
            print self._buildScript(method)
            print "# submit with:"
            print "#", self._buildBatchCommand()
            return True
        
        # write the script
        scriptFile = open(expandMacros("${script}", self.inv), "w")
        scriptFile.write(self._buildScript(method))
        scriptFile.close()

        # build the batch command
        command = self._buildBatchCommand()
        self._info.log("executing: {%s}" % command)

        os.system(command)
        return True


    def __init__(self, name):
        Scheduler.__init__(self, name)
        
        # Used to recursively expand ${macro) in format strings using my inventory.
        class InventoryAdapter(object):
            def __init__(self, launcher):
                self.launcher = launcher
            def __getitem__(self, key):
                return expandMacros(str(self.launcher.inventory.getTraitValue(key)), self)
        self.inv = InventoryAdapter(self)
        
        return

    
    def _configure(self):
        self.nodes = self.inventory.nodes


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} ${script}", self.inv)


    def _buildScript(self, method):
        script = [
            "#!/bin/sh",
            ]
        self._buildScriptDirectives(script)
        args = [self.inventory.executable] + self.inventory.arguments
        args.append("--command=" + method)
        args = ' '.join(args)
        script += [
            expandMacros('''\

cd ${directory}
''' + args, self.inv)
            ]
        script = "\n".join(script) + "\n"
        return script


# Note: mpi.LauncherPBS in Pythia-0.8 does not work!

class SchedulerPBS(BatchScheduler):


    class Inventory(BatchScheduler.Inventory):
        
        from pyre.inventory import str
        
        command = str("command", default="mpirun -np ${nodes} -machinefile $PBS_NODEFILE") # Sub-launcher?
        batch_command = str("batch-command", default="qsub")


    def __init__(self):
        BatchScheduler.__init__(self, "pbs")


    def _buildScriptDirectives(self, script):
        
        queue = self.inventory.queue
        if queue:
            script.append("#PBS -q %s" % queue)

        task = self.inventory.task
        if task:
            script.append("#PBS -N %s" % task)

        if self.inventory.stdout:
            script.append(expandMacros("#PBS -o ${stdout}", self.inv))
        if self.inventory.stderr:
            script.append(expandMacros("#PBS -e ${stderr}", self.inv))

        resourceList = self._buildResourceList()

        script += [
            "#PBS -V", # export qsub command environment to the batch job
            "#PBS -l %s" % resourceList,
            ]

        return script



    def _buildResourceList(self):

        resourceList = [
            "nodes=%d" % self.nodes,
            ]

        walltime = self.inventory.walltime.value
        if walltime:
            resourceList.append("walltime=%d:%02d:%02d" % hms(walltime))

        resourceList = ",".join(resourceList)

        return resourceList


class SchedulerLSF(BatchScheduler):


    class Inventory(BatchScheduler.Inventory):
        
        from pyre.inventory import list, str
        
        command = str("command", default="mpijob mpirun")
        batch_command = str("batch-command", default="bsub")
        bsub_options = list("bsub-options")


    def __init__(self):
        BatchScheduler.__init__(self, "lsf")


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} < ${script}", self.inv)


    def _buildScriptDirectives(self, script):

        # LSF scripts must have a job name; otherwise strange "/bin/sh: Event not found"
        # errors occur (tested on TACC's Lonestar system).
        task = self.inventory.task
        if not task:
            task = "jobname"
        script.append("#BSUB -J %s" % task)
        
        queue = self.inventory.queue
        if queue:
            script.append("#BSUB -q %s" % queue)

        walltime = self.inventory.walltime.value
        if walltime:
            script.append("#BSUB -W %d:%02d" % hms(walltime)[0:2])
        
        if self.inventory.stdout:
            script.append(expandMacros("#BSUB -o ${stdout}", self.inv))
        if self.inventory.stderr:
            script.append(expandMacros("#BSUB -e ${stderr}", self.inv))
            
        script += [
            "#BSUB -n %d" % self.nodes,
            ]

        script += ["#BSUB " + option for option in self.inventory.bsub_options]

        return script


class SchedulerGlobus(BatchScheduler):


    class Inventory(BatchScheduler.Inventory):

        from pyre.inventory import str
        
        batch_command = str("batch-command", default="globusrun")
        resource = str("resource", default="localhost")


    def _buildBatchCommand(self):
        return expandMacros("${batch-command} -b -r ${resource} -f ${script}", self.inv)


    def __init__(self):
        BatchScheduler.__init__(self, "globus")


    def _buildScript(self, method):
        script = [
            expandMacros('''\
&   (jobType=mpi)
    (executable="${executable}")
    (count=${nodes})
    (directory="${directory}")
    (stdout="${stdout}")
    (stderr="${stderr}")''', self.inv),
            ]
        
        script.append('    (environment = %s)' % self._buildEnvironment())

        # add the arguments
        args = self.inventory.arguments
        args.append("--command=" + method)
        command = '    (arguments= ' + ' '.join([('"%s"' % arg) for arg in args]) + ')'
        script.append(command)

        script = '\n'.join(script) + '\n'

        return script


    def _buildEnvironment(self):
        from os import environ
        #vars = environ.keys()
        vars = self.inventory.environment
        env = [('(%s "%s")' % (var, environ.get(var,""))) for var in vars]
        env = ' '.join(env)
        return env


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SchedulerFacility
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


builtInSchedulerClasses = {
    'none':    SchedulerNone,
    'kernel':  KernelScheduler,
    'globus':  SchedulerGlobus,
    'lsf':     SchedulerLSF,
    'pbs':     SchedulerPBS,
    }

class SchedulerFacility(Facility):

    def __init__(self, name, default=None):
        if default is None:
            default = "none"
        Facility.__init__(self, name, default=default)
        return
    
    def _retrieveComponent(self, instance, componentName):
        cls = builtInSchedulerClasses.get(componentName, None)
        if cls is None:
            return Facility._retrieveComponent(self, instance, componentName)
        scheduler = cls()
        scheduler.aliases.append(self.name)
        import pyre.parsing.locators
        locator = pyre.parsing.locators.simple('built-in')
        return scheduler, locator


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# main
if __name__ == "__main__":
    from pyre.applications.Script import Script
    class TestApp(Script):
        class Inventory(Script.Inventory):
            scheduler = SchedulerFacility("scheduler", default="none")
        def main(self, *args, **kwds):
            scheduler = self.inventory.scheduler
            if scheduler:
                scheduler.inventory.dry = True
                scheduler.schedule(self, 'launch')
            return
        def __init__(self):
            Script.__init__(self, "TestApp")
    
    app = TestApp()
    app.run()


# end of file 
