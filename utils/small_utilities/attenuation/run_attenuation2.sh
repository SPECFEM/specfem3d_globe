#!/bin/sh

t1="1000.0" # Max Period
t2="20.0"   # Min Period
n="3"       # Number of Standard Linear Solids
Q_IC=84.6   # Inner Core
Q_LM=312.0  # Lower Mantle
Q_UM=143.0  # Upper Mantle
Q_LVZ=80.0  # Low Velocity Zone
Q_S=600.0   # Surface

LD_LIBRARY_PATH=/opt/intel/lib
./attenuation_test <<EOF
$t1
$t2
$n
$Q_IC
$Q_LM
$Q_UM
$Q_LVZ
$Q_S
EOF

