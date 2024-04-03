# Nom arbitraire du travail LoadLeveler
# @ job_name = Maxwell_3D_FDTD_PML_OpenMP

# Type de travail
# @ job_type = serial

# Choix de l'interpreteur de commande
# @ shell = /bin/ksh

# Fichier de sortie standard du travail
# @ output = $(job_name).$(jobid)

# Fichier de sortie d'erreur du travail
# @ error =  $(job_name).$(jobid)

# Temps CPU max. en seconde (pour 1 heure elapsed, compter OMP_NUM_THREADS heures)
# @ cpu_limit = 43200

# Memoire max. en data utilisee
# @ data_limit = 11.5Gb

# Memoire max. en stack utilisee
# @ stack_limit = 1.2Gb,1.2Gb

# Nombre de processeurs a affecter aux threads
# OpenMP (ici 4, voir plus bas la variable OMP_NUM_THREADS).
# @ resources = ConsumableCpus(4)

# @ queue

# Pour avoir l'echo des commandes
set -x

# Repertoire temporaire de travail
################cd $TMPDIR

# La variable LOADL_STEP_INITDIR est automatiquement positionnee par
# LoadLeveler au repertoire dans lequel on tape la commande llsubmit
################cp $LOADL_STEP_INITDIR/source.f .

# Repertoire temporaire de travail
cd $HOME/code_3D/with_PML_OpenMP_4tasks

# Compilation et edition de liens d'un programme OpenMP au format libre
rm -f xonde3D
xlf_r -qsmp=omp -O4 -qfree=f90 -qsuffix=f=f90 -o xonde3D onde3d_mathieu_maxwell_PML_12oct2005.f90

# La memoire STACK max. (defaut 4Mo) utilisee (ici 64 Mo) par
# les variables privees de chaque thread
export XLSMPOPTS=stack=65536000

# Variables d'environnement indiquant le nombre de threads OpenMP
# (indiquer une valeur identique a celle positionnee plus haut
# dans la directive threads_per_task)
export OMP_NUM_THREADS=4

# Execution du programme OpenMP
./xonde3D

