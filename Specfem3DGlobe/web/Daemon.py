#!/usr/bin/env python

import urllib2
import shutil

# 
# source address
#
urlPrefix           = 'http://localhost:8000/specfem3dglobe/'
simulationsUrl      = urlPrefix + 'simulations/list.py'
inputFileUrl        = urlPrefix + 'simulations/%s/%s'
inputFiles          = ['parameters.pml', 'events.txt', 'stations.txt']

#
# destination 
#
sSaveDir            = './output'


infile = urllib2.urlopen(simulationsUrl)
simulations = infile.read()
infile.close()
simulations = eval(simulations)


def newSimulation(sim):
    id = sim['id']
    status = sim['status']
    for inputFile in inputFiles:
        infile = urllib2.urlopen(inputFileUrl % (id, inputFile))
        outfile = open('%s/%d-%s' % (sSaveDir, id, inputFile), 'w')
        shutil.copyfileobj(infile,  outfile)
        outfile.close()
        infile.close()
    print "run simulation", id
    return

#
# for each sim, get input files and save
#
for sim in simulations:
    switch = {
        'new': lambda sim: newSimulation(sim),
        }
    switch[sim['status']](sim)


# Here is an example session that shows how to "POST" requests:

# >>> import httplib, urllib
# >>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
# >>> headers = {"Content-type": "application/x-www-form-urlencoded",
# ...            "Accept": "text/plain"}
# >>> conn = httplib.HTTPConnection("musi-cal.mojam.com:80")
# >>> conn.request("POST", "/cgi-bin/query", params, headers)
# >>> response = conn.getresponse()
# >>> print response.status, response.reason
# 200 OK
# >>> data = response.read()
# >>> conn.close()


# end of file
