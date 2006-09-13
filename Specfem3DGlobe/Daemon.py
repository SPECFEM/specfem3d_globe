#!/usr/bin/env python


multipartBoundary   = '----------eArThQuAkE$'


# 0         1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME

class Job(object):
    
    def __init__(self, jobId, stat, jobName, simulationNumber):
        self.jobId = jobId
        self.stat = stat
        self.jobName = jobName
        self.simulationNumber = simulationNumber

    def __repr__(self): return repr(self.__dict__)


def bjobs():
    from os import popen
    stream = popen('bjobs -a', 'r')
    jobs = {}
    for line in stream:
        jobId = line[0:7].strip()
        if jobId == 'JOBID' or jobId == '':
            continue
        jobId = int(jobId)
        stat = line[16:21].strip()
        jobName = line[57:67].strip()
        if jobName.startswith('web-'):
            simulationNumber = int(jobName.split('-')[1])
            jobs[simulationNumber] = Job(jobId, stat, jobName, simulationNumber)
    stream.close()
    return jobs


class Sim(object):
    def __init__(self, id, job):
        self.id = id
        self.job = job


# /!\ These SimStatus* class names are stored in the database on the web server!

class SimStatus(object):

    def __init__(self, daemon, sim):
        self.daemon = daemon
        self.portal = daemon.portal
        self.sim = sim
    
    def display(cls):
        """Return the string displayed to the web user."""
        raise NotImplementedError()
    display = classmethod(display)

    def poll(self):
        stat = ""
        if self.sim.job:
            stat = self.sim.job.stat
        self.daemon._info.log("simulation %d %s %s" % (self.sim.id, self.display(), stat))
        self._poll()

    def _poll(self): pass
    
    def postStatusChange(self, newStatus):
        import urllib
        id = self.sim.id
        body = urllib.urlencode({'status': newStatus.__name__})
        headers = {"Content-Type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain"}
        self.post(id, body, headers)
        return

    def postStatusChangeAndUploadOutput(self, newStatus, output):
        # Based upon a recipe by Wade Leftwich.
        id = self.sim.id
        fields = {'status': newStatus.__name__}
        files = [('output', 'output.txt', 'application/x-gtar', 'gzip', output)]
        body = self.encodeMultipartFormData(fields, files)
        headers = {"Content-Type": "multipart/form-data; boundary=%s" % multipartBoundary,
                   "Content-Length": str(len(body)),
                   "Accept": "text/plain"}
        self.post(id, body, headers)
        return

    def encodeMultipartFormData(self, fields, files):
        import shutil
        from StringIO import StringIO
        stream = StringIO()
        def line(s=''): stream.write(s + '\r\n')
        for key, value in fields.iteritems():
            line('--' + multipartBoundary)
            line('Content-Disposition: form-data; name="%s"' % key)
            line()
            line(value)
        for (key, filename, contentType, contentEncoding, content) in files:
            line('--' + multipartBoundary)
            line('Content-Disposition: form-data; name="%s"; filename="%s"' % (key, filename))
            line('Content-Type: %s' % contentType)
            line('Content-Encoding: %s' % contentEncoding)
            line()
            shutil.copyfileobj(content, stream)
            line()
        line('--' + multipartBoundary + '--')
        line()
        return stream.getvalue()

    def post(self, id, body, headers):
        self.daemon._info.log("POST %d" % id)
        import httplib
        conn = httplib.HTTPConnection(self.portal.host)
        conn.request("POST", (self.portal.simStatusUrl % id), body, headers)
        response = conn.getresponse()
        data = response.read()
        conn.close()
        self.daemon._info.log("response %s" % data)
        return


class SimStatusNew(SimStatus):
    def display(cls): return 'new'
    display = classmethod(display)

    def _poll(self):
        import os
        import urllib2
        import shutil
        from os.path import isdir, join
        
        id = self.sim.id
        portal = self.portal

        simDir = join(self.daemon.simulationRoot, str(id))
        if isdir(simDir):
            pass
        else:
            os.makedirs(simDir)

        def localCopy(inputFile):
            copy = join(simDir, inputFile)
            infile = urllib2.urlopen(portal.inputFileUrl % (id, inputFile))
            outfile = open(copy, 'wb')
            shutil.copyfileobj(infile,  outfile)
            outfile.close()
            infile.close()
            return copy
            
        simulation = Simulation()
        simulation.daemonArgv = [
            localCopy('parameters.pml'),
            '--output-dir=' + simDir,
            '--solver.cmt-solution=' + localCopy('events.txt'),
            '--solver.stations=' + localCopy('stations.txt'),
            '--solver.seismogram-archive=' + join(simDir, 'output.tar.gz'),
            '--job.name=web-%d' % id,
            ]
        
        # schedule
        try:
            simulation.run()
        except Exception, e:
            self.daemon._error.log("simulation %d: error: %s" % (self.sim.id, e))
            self.postStatusChange(SimStatusDone)
            return
        
        self.postStatusChange(SimStatusPending)
        
        return


class SimStatusPending(SimStatus):
    def display(cls): return 'pending'
    display = classmethod(display)

    def _poll(self):
        if not self.sim.job:
            pass
        elif self.sim.job.stat == 'PEND':
            return
        elif self.sim.job.stat == 'RUN':
            self.postStatusChange(SimStatusRunning)
            return
        # The job passed through the 'RUN' state without our noticing.
        status = SimStatusRunning(self.daemon, self.sim)
        status.poll()
        return


class SimStatusRunning(SimStatus):
    def display(cls): return 'running'
    display = classmethod(display)

    def _poll(self):
        from os.path import join, getsize
        if not self.sim.job:
            pass
        elif self.sim.job.stat == 'RUN':
            return
        elif self.sim.job.stat != 'DONE':
            # error
            self.postStatusChange(SimStatusDone)
            return
        id = self.sim.id
        simDir = join(self.daemon.simulationRoot, str(id))
        outputName = join(simDir, 'output.tar.gz')
        size = 0
        try:
            size = getsize(outputName)
        except Exception, e:
            self.daemon._error.log("simulation %d: error: %s" % (self.sim.id, e))
        if size > 0:
            output = open(outputName, 'rb')
            self.postStatusChangeAndUploadOutput(SimStatusDone, output)
            output.close()
        else:
            # error
            self.postStatusChange(SimStatusDone)
        return


class SimStatusDone(SimStatus):
    def display(cls): return 'done'
    display = classmethod(display)

    def poll(self): pass


simStatusList = [
    SimStatusNew,
    SimStatusPending,
    SimStatusRunning,
    SimStatusDone,
    ]

# /!\ This is imported by the Django app!
STATUS_CHOICES = tuple([(cls.__name__, cls.display()) for cls in simStatusList])


from cig.addyndum.applications import Script
from cig.addyndum.components import Component
from Specfem3DGlobe.Specfem import Specfem


class Simulation(Specfem):

    def __init__(self):
        super(Simulation, self).__init__()
        self.daemonArgv = None
    
    def processCommandline(self, registry, parser=None):
        import sys
        if parser is None:
            parser = self.createCommandlineParser()
        if self.daemonArgv is None:
            argv = sys.argv
        else:
            argv = self.daemonArgv
        help, unprocessedArguments = parser.parse(registry, argv=argv)
        return help, unprocessedArguments


class WebPortal(Component):
    
    componentName = "Specfem3DGlobe-WebPortal"

    import pyre.inventory as pyre
    host     = pyre.str("host", default="localhost:8000")
    urlRoot  = pyre.str("url-root", default="/specfem3dglobe/")

    def _configure(self):
        self.urlPrefix           = 'http://%s%s' % (self.host, self.urlRoot)
        self.simulationsUrl      = self.urlPrefix + 'simulations/list.py'
        self.inputFileUrl        = self.urlPrefix + 'simulations/%d/%s'
        self.simStatusUrl        = self.urlRoot + 'simulations/%d/status/'


class Daemon(Script):

    componentName = "Specfem3DGlobe-Daemon"

    import cig.addyndum.inventory as addyndum
    portal = addyndum.facility("portal", factory=WebPortal)
    simulationRoot = addyndum.outputDir("simulation-root", default="output")

    def main(self, *args, **kwds):
        self._info.activate()
        self._info.log("~~~~~~~~~~ daemon started ~~~~~~~~~~")
        import urllib2
        infile = urllib2.urlopen(self.portal.simulationsUrl)
        simulations = infile.read()
        infile.close()
        simulations = eval(simulations)
        jobs = bjobs()
        for sim in simulations:
            cls = globals()[sim['status']]
            id = int(sim['id'])
            job = jobs.get(id)
            sim = Sim(id, job)
            status = cls(self, sim)
            status.poll()
        return


def main(*args, **kwds):
    daemon = Daemon()
    daemon.run(*args, **kwds)


# end of file
