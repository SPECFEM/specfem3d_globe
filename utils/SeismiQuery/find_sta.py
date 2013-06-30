#!/usr/bin/python
# this script converts the output downloaded from IRIS SeismicQuery to
# a proper STATIONS file  - Qinya Liu
import os
lines=file("lh2.txt",'r').readlines()
for i in range(0,len(lines)):
  sp = lines[i].split()
  str = "egrep '%s .*%s' lh.txt " % (sp[0],sp[1])
  sp2 = os.popen(str).read().split()
  print "%-10s %-5s %10.4f  %10.4f  %8.1f  %6.1f" % (sp2[1],sp2[0],float(sp2[7]),float(sp2[8]),float(sp2[9]),float(sp2[10]))

