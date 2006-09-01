#!/usr/bin/env python


multipartBoundary   = '----------eArThQuAkE$'


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

    def poll(self): pass
    
    def postStatusChange(self, newStatus):
        import urllib
        id = self.sim['id']
        body = urllib.urlencode({'status': newStatus.__name__})
        headers = {"Content-Type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain"}
        self.post(id, body, headers)
        return

    def postStatusChangeAndUploadOutput(self, newStatus, output):
        # Based upon a recipe by Wade Leftwich.
        id = self.sim['id']
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
        import httplib
        conn = httplib.HTTPConnection(self.portal.host)
        conn.request("POST", (self.portal.simStatusUrl % id), body, headers)
        response = conn.getresponse()
        data = response.read()
        conn.close()
        return


class SimStatusNew(SimStatus):
    def display(cls): return 'new'
    display = classmethod(display)

    def poll(self):
        import os
        import urllib2
        import shutil
        from os.path import isdir, join
        
        id = self.sim['id']
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
            
        simulation = Simulation([
            localCopy('parameters.pml'),
            '--output-dir=' + simDir,
            '--solver.cmt-solution=' + localCopy('events.txt'),
            '--solver.stations=' + localCopy('stations.txt'),
            '--solver.seismogram-archive=' + join(simDir, 'output.tar.gz'),
            '--scheduler.dry',
            '--launcher.dry',
            ])
        
        # schedule
        simulation.run()
        
        self.postStatusChange(SimStatusPending)
        
        return


class SimStatusPending(SimStatus):
    def display(cls): return 'pending'
    display = classmethod(display)

    def poll(self):
        id = self.sim['id']
        print "check job status (waiting)", id
        self.postStatusChange(SimStatusRunning)
        return


class SimStatusRunning(SimStatus):
    def display(cls): return 'running'
    display = classmethod(display)

    def poll(self):
        from os.path import join
        id = self.sim['id']
        print "check job status (running)", id
        simDir = join(self.daemon.simulationRoot, str(id))
        outputName = join(simDir, 'output.tar.gz')
        output = open(outputName, 'rb')
        self.postStatusChangeAndUploadOutput(SimStatusDone, output)
        output.close()
        return


class SimStatusDone(SimStatus):
    def display(cls): return 'done'
    display = classmethod(display)


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

    def __init__(self, daemonArgv):
        super(Simulation, self).__init__()
        self.daemonArgv = daemonArgv
    
    def processCommandline(self, registry, parser=None):
        if parser is None:
            parser = self.createCommandlineParser()
        help, unprocessedArguments = parser.parse(registry, argv=self.daemonArgv)
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
        self.simStatusUrl        = self.urlRoot + 'simulations/%s/status/'


class Daemon(Script):

    componentName = "Specfem3DGlobe-Daemon"

    import cig.addyndum.inventory as addyndum
    portal = addyndum.facility("portal", factory=WebPortal)
    simulationRoot = addyndum.outputDir("simulation-root", default="output")

    def main(self, *args, **kwds):
        import urllib2
        infile = urllib2.urlopen(self.portal.simulationsUrl)
        simulations = infile.read()
        infile.close()
        simulations = eval(simulations)
        for sim in simulations:
            cls = globals()[sim['status']]
            status = cls(self, sim)
            status.poll()
        return


def main(*args, **kwds):
    daemon = Daemon()
    daemon.run(*args, **kwds)


# end of file
