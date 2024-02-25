#!/usr/bin/python

# requesting data for events and convert them to sac format
# written by Qinya Liu, Nov 2009, for global shakemovie website

import os, sys, glob, getopt, datetime

usage='\nUsage: dhi2mseed.py -m CMT -d duration (min) [-j DHI2mseed.jar_file -p iris.txt -s station_file -o out_datadir]\n'

try:
  opts,arg=getopt.getopt(sys.argv[1:],'j:m:d:p:s:o:')
except getopt.GetoptError:
  sys.exit(usage)
if (len(opts) == 0):
  sys.exit(usage)

# dir where _GSN.dataless or II/IU.dataless is placed
dataless_dir='/data2/Datalib/stations/'
jar_file='/data2/Datalib/iris_request/DHI2mseed.jar'
iris_file='/data2/Datalib/iris_request/iris.txt'
sta_file='/data2/Datalib/stations/STATIONS_LH_II_IU'
out_dir='data'

for o,a in opts:
  if o == '-j':
    jar_file=a
  if o == '-m':
    cmt_file=a
  if o == '-d':
    dur=float(a)
  if o == '-p':
    iris_file=a
  if o == '-s':
    sta_file=a
  if o == '-o':
    out_dir=a

if (not os.path.isfile(jar_file)):
  sys.exit(jar_file+' does not exist')
if (not os.path.isfile(cmt_file)):
  sys.exit(cmt_file+' does not exist')
if (not os.path.isfile(iris_file)):
  sys.exit(iris_file+' does not exist')
if (not os.path.isfile(sta_file)):
  sys.exit(sta_file+' does not exist')
if (not os.path.isdir(out_dir)):
  os.mkdir(out_dir)

# write station request file
f=open(out_dir+'/scnl.list','w')
lines = file(sta_file, 'r').readlines()
for i in range(0,len(lines)):
  sp = lines[i].split()
#  print "%s LH* %s *" % (sp[0], sp[1])
  f.write("%s LH* %s *\n" % (sp[0], sp[1]))
f.close()

# calculate request start time
[year,mon,day,hr,min,sec]=os.popen("awk 'NR == 1 {print $0}' "+cmt_file).readline().split()[1:7]

starttime=(datetime.datetime(int(year), int(mon), int(day), int(hr), int(min), int(float(sec)))+datetime.timedelta(minutes=-3)).isoformat(' ')

command='java -jar '+jar_file+' -starttime "'+starttime+'" -duration '+str(dur)+' -p '+iris_file+' -f '+out_dir+'/scnl.list -o '+out_dir

print '\n**** requesting data from IRIS ***'
print command+'\n'
if (os.system(command+' |grep writingi 1>&2') != 0):
  sys.exit('Error requesting data')

print '\n**** convert miniSEED files to sac files ****'
if (not os.path.isdir(dataless_dir)):
  sys.exit(dataless_dir+' does not exist')

# use different dataless seed file for files from different networks
net_input={};
for file in glob.glob(out_dir+'/*.mseed'):
  bfile=os.path.basename(file)
  net=bfile.split('.')[2]
  if (net in net_input.keys()):
    net_input[net]=net_input[net]+bfile+'\n\n\nd\n\n\n\n\n\n\n\n\n3\n\n\n\n\n\n'
  else:
    net_input[net]=bfile+'\n\n\nd\n\n\n\n\n\n\n\n\n3\n\n\n\n\n\n'

for net,ni in net_input.iteritems():
  dataless_file=dataless_dir+'/'+net+'.dataless'
#  dataless_file=dataless_dir+'/_GSN.dataless'
  if (not os.path.isfile(dataless_file)):
    print dataless_file+' does not exist'; exit()
  print 'Using dataless file ' + dataless_file+' for network '+ net

  os.environ['ALT_RESPONSE_FILE']=dataless_file

#  print 'cd '+out_dir+'; rdseed <<EOF\n'+ni+'Quit\nEOF\n' # debug
  if (os.system('cd '+out_dir+'; rdseed <<EOF 2> /dev/null | grep Writing 1>&2 \n'+ni+'Quit\nEOF\n') != 0):
    sys.exit('Error converting miniseed to sac for '+net)

# move original mseed files to a sub-directory
os.system('rm -f scnl.list; cd '+out_dir+'; mkdir -p mseed_files; mv *.mseed mseed_files/')
