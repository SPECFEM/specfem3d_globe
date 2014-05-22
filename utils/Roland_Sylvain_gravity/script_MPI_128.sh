#!/bin/bash

########### -n : number of tasks that will be used in parallel mode (default=1)
########### -N : number of nodes to allocate for parallel usage (default is chosen by the underlying system)
########### -c : number of cores per parallel task to allocate (default=1)
#MSUB -n 600     # Reservation des coeurs de processeurs
#MSUB -N 38     # nombre total de noeuds a utiliser (mettre n/16 car il y a 16 coeurs par noeud)
#MSUB -c  1     # 1 coeur par tache MPI

#################### #MSUB -E '-R select[rank!=1004&rank!=1005]'

##################MSUB -A gen7165
#MSUB -A gen6351
#MSUB -r SPECFEM128          # Nom du job                
#MSUB -T 86400                # Limite de temps elapsed du job en secondes
#BSUB -x              # mode exclusif
#MSUB -o output_night_%I.o          # Sortie standard et %I est le job_ID
#MSUB -e output_night_%I.e          # Sortie d'erreur et %I est le job_ID
##############MSUB -j oe                 # join output and error
#MSUB -k n                  # do not keep older output and error from previous jobs with the same name

set -x

# ${BRIDGE_MSUB_PWD} est une variable d'environnement representant le repertoire de soumission
cd ${BRIDGE_MSUB_PWD}

#export CUDA_PROFILE=1
#export CUDA_PROFILE_LOG=dimitri_CUDA_profile.log
#export CUDA_PROFILE_CONFIG=dimitri_CUDA_profile_config.txt

##########################

ccc_mprun ./bin/xmeshfem3D

