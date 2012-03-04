#!/usr/bin/python
# this script finds the corresponding name info for STATIONS_LH from downloaded
# IRIS station info file
# - Qinya Liu
import os
t=" "
lines=file("lh2.txt",'r').readlines()
for i in range(0,len(lines)):
  sp = lines[i].split()
  str = "egrep '^%s .*%s ' stations_info.txt " % (sp[0],sp[1])
  error=os.system(str+">/dev/null")
  if (error != 0):
    print "error ====" + t.join(sp)
  else:
    sp2=os.popen(str + " | awk 'NR == 1 {print $0}'").read().split()
    print "%-8s %-5s %9.4f  %9.4f  %-40s" % (sp2[1],sp2[0],float(sp2[2]),float(sp2[3]),t.join(sp2[4:]))

