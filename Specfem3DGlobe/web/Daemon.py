#!/usr/bin/env python


# 
# source address
#
host                = 'localhost:8000'
urlRoot             = '/specfem3dglobe/'
urlPrefix           = 'http://%s%s' % (host, urlRoot)
simulationsUrl      = urlPrefix + 'simulations/list.py'
inputFileUrl        = urlPrefix + 'simulations/%d/%s'
inputFiles          = ['parameters.pml', 'events.txt', 'stations.txt']
simStatusUrl        = urlRoot + 'simulations/%s/status/'

#
# destination 
#
sSaveDir            = './output'


# /!\ These SimStatus* class names are stored in the database on the web server!

class SimStatus(object):
    
    def display(cls):
        """Return the string displayed to the web user."""
        raise NotImplementedError()
    display = classmethod(display)

    def poll(self, sim): pass
    
    def postStatusChange(self, sim, newStatus):
        import httplib
        import urllib
        id = sim['id']
        params = urllib.urlencode({'status': newStatus.__name__})
        headers = {"Content-type": "application/x-www-form-urlencoded",
                   "Accept": "text/plain"}
        conn = httplib.HTTPConnection(host)
        conn.request("POST", (simStatusUrl % id), params, headers)
        response = conn.getresponse()
        data = response.read()
        conn.close()
        return


class SimStatusNew(SimStatus):
    def display(cls): return 'new'
    display = classmethod(display)

    def poll(self, sim):
        # get input files and save
        import urllib2
        import shutil
        id = sim['id']
        for inputFile in inputFiles:
            infile = urllib2.urlopen(inputFileUrl % (id, inputFile))
            outfile = open('%s/%d-%s' % (sSaveDir, id, inputFile), 'w')
            shutil.copyfileobj(infile,  outfile)
            outfile.close()
            infile.close()
        print "run simulation", id
        self.postStatusChange(sim, SimStatusPending)
        return


class SimStatusPending(SimStatus):
    def display(cls): return 'pending'
    display = classmethod(display)

    def poll(self, sim):
        id = sim['id']
        print "check job status (waiting)", id
        self.postStatusChange(sim, SimStatusRunning)
        return


class SimStatusRunning(SimStatus):
    def display(cls): return 'running'
    display = classmethod(display)

    def poll(self, sim):
        id = sim['id']
        print "check job status (running)", id
        self.postStatusChange(sim, SimStatusDone)
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


def main():
    import urllib2
    infile = urllib2.urlopen(simulationsUrl)
    simulations = infile.read()
    infile.close()
    simulations = eval(simulations)
    for sim in simulations:
        cls = globals()[sim['status']]
        status = cls()
        status.poll(sim)
    return

if __name__ == '__main__':
    main()


# end of file
