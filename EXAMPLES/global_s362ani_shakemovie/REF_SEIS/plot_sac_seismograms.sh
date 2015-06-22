#!/bin/bash

sac <<EOF
r ./IU.GRFO.MXZ.sem.sac ./IU.YAK.MXZ.sem.sac ./IU.COLA.MXZ.sem.sac ./IU.RAIO.MXZ.sem.sac ./IU.RAO.MXZ.sem.sac
p1
saveimg seismograms.pdf
quit
EOF

echo ""
echo "plotted: seismograms.pdf"
echo

