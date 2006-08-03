#!/usr/bin/env python

import urllib2
import pysqlite2.dbapi2 as sqlite

# 
# database related
#
sDatabaseFile       = '../../data/data.db'
sTableName          = 'Specfem3DGlobe_simulation'
nReady              = 2
nPending            = 3

# 
# source address
#
sUrlPrefix          = 'http://localhost:8000/specfem3dglobe/simulations/'
sUrlPostfix         = '.pml'
lsUpdateID          = []

#
# destination 
#
sSaveDir            = './'
sFileExt            = '.xml'
sMode               = 'w'

#
# read 'ready' records from database, and save them in list
# for now, let's assume that we have 1 as ready sim id.
#
conn = sqlite.connect(sDatabaseFile)
cur = conn.cursor()
cur.execute('select id from %s where status = %d' % (sTableName,nReady)) 

#
# for each sim id, get xml file and save
#
for id in cur:
    try:
        infile = urllib2.urlopen('%s%d%s' % (sUrlPrefix,id[0],sUrlPostfix))
    except urllib2.HTTPError:
        # print error to log
        pass
    else:
        outfile = open('%s%d%s' % (sSaveDir,id[0],sFileExt),sMode)
        outfile.write(infile.read())
        infile.close()
        outfile.close()
        lsUpdateID.append(id[0])

for id in lsUpdateID:
    cur.execute('update %s set status = %d where id = %d' % (sTableName,nPending,id))

conn.commit()
cur.close()
conn.close()
