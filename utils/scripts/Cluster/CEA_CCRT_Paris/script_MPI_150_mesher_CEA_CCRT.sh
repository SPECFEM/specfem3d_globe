#!/bin/bash

########### -n : number of tasks that will be used in parallel mode (default=1)
########### -N : number of nodes to allocate for parallel usage (default is chosen by the underlying system)
########### -c : number of cores per parallel task to allocate (default=1)
#MSUB -n 150      # Reservation des coeurs de processeurs
#MSUB -N 19      # nombre total de noeuds a utiliser (mettre n/8 car on prend 4 coeurs de CPUs par noeud)
#MSUB -c  1     # 1 coeur par tache MPI

#MSUB -r MESHFEM150          # Nom du job
#################################MSUB -T 3600                # Limite de temps elapsed du job en secondes
#MSUB -T 1799                # Limite de temps elapsed du job en secondes
#BSUB -x              # mode exclusif
#MSUB -o output_night_mesher_150.o          # Sortie standard
#MSUB -e output_night_mesher_150.e          # Sortie d'erreur
##############MSUB -@ dimitri.komatitsch@univ-pau.fr:end    # envoie un mail a l'adresse indiquee en fin de job
##############MSUB -j oe                 # join output and error
#MSUB -k n                  # do not keep older output and error from previous jobs with the same name

set -x
cd $BRIDGE_MSUB_PWD

echo $PBS_NODEFILE

##########################

mpirun -np 150 ./xmeshfem3D

