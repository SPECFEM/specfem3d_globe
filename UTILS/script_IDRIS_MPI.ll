
# Choix du shell
# @ shell = /bin/ksh

# Nombre de processus MPI demandés
# @ total_tasks = 24

# Temps CPU max. par processus MPI en hh:mm:ss
# @ cpu_limit = 00:55:00

# Mémoire max. utilisée par processus
# @ data_limit = 1100mb
###################### @ data_limit   = 9.2Gb

###################### @ stack_limit  = 6.2Gb,6.2Gb

# Nom arbitraire du travail LoadLeveler
# @ job_name = run_SPECFEM3D_acoustic_MPI

#----------------------------------------------------

# Type de travail
# @ job_type = parallel

# Fichier de sortie standard du travail
# @ output = $(job_name).$(jobid)

# Fichier de sortie d'erreur du travail
# @ error =  $(job_name).$(jobid)

# @ queue

#----------------------------------------------------

# Pour avoir l'écho des commandes
set -x

##### nom du repertoire ou est stocke le code dans le home
##### et nom du sous-repertoire ou est stockee la base de donnees
export repertoire_code=SPHERE_ACOUST_3D
export repertoire_database=DATABASES_MPI_DIMITRI

# vider les repertoires dans le home
rm -r -f $HOME/$repertoire_code/OUTPUT_FILES $HOME/$repertoire_code/$repertoire_database

# copier les codes source depuis le home vers le repertoire temporaire
rm -r -f $TMPDIR/$repertoire_code
cp -r -p $HOME/$repertoire_code $TMPDIR

# creer les nouveaux repertoires temporaires
mkdir $TMPDIR/$repertoire_code/OUTPUT_FILES $TMPDIR/$repertoire_code/$repertoire_database

cd $TMPDIR/$repertoire_code

# compiler et executer le mailleur en MPI
make clean
make meshfem3D
./xmeshfem3D

# compiler et executer le solver en MPI
make clean
make specfem3D
./xspecfem3D

# supprimer la base de donnees creee
rm -r -f $TMPDIR/$repertoire_code/$repertoire_database

# recuperer le job ID
export myjobid=$( echo $LOADL_STEP_ID | cut -d'.' -f4 )

# deplacer tous les resultats dans le workdir en ajoutant le job ID au nom
# sortir d'abord dans home pour pouvoir detruire le repertoire courant de tmpdir
cd $HOME
mv $TMPDIR/$repertoire_code $WORKDIR/${repertoire_code}_${myjobid}

