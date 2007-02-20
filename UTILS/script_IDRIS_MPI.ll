
# Nombre de processus MPI demand�s
# @ total_tasks = 8

# Temps CPU max. en seconde par processus
# @ cpu_limit = 400

# M�moire max. utilis�e par processus
# @ data_limit = 700mb
###################### @ data_limit   = 9.2Gb

###################### @ stack_limit  = 6.2Gb,6.2Gb

##### nom du repertoire ou est stocke le code dans le home
REPERTOIRE_CODE = SPECFEM3D_BASIN_Carcione

#----------------------------------------------------

# Nom arbitraire du travail LoadLeveler
# @ job_name = run_SPECFEM3D_Carcione_Copper_MPI

# Type de travail
# @ job_type = parallel

# Fichier de sortie standard du travail
# @ output = $(job_name).$(jobid)

# Fichier de sortie d'erreur du travail
# @ error =  $(job_name).$(jobid)

# @ queue

#----------------------------------------------------

# Pour avoir l'�cho des commandes
set -x

# R�pertoire temporaire de travail
cd $TMPDIR

################# La variable LOADL_STEP_INITDIR est automatiquement positionn�e par
################# LoadLeveler au r�pertoire dans lequel on tape la commande llsubmit
################cp $LOADL_STEP_INITDIR/a.out .

cp -r -p $HOME/$REPERTOIRE_CODE .
cd $REPERTOIRE_CODE

# vider les repertoires
rm -r -f OUTPUT_FILES DATABASES_MPI_DIMITRI
mkdir OUTPUT_FILES DATABASES_MPI_DIMITRI

# compiler et executer le mailleur en MPI
make clean
make meshfem3D
./xmeshfem3D

# compiler et executer le solver en MPI
make clean
make specfem3D
./xspecfem3D

