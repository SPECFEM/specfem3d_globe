
# Choix du shell
# @ shell = /bin/ksh

# Nombre de processus MPI demandes
# @ total_tasks = 24

# Temps CPU maximum par processus MPI en hh:mm:ss
# @ cpu_limit = 04:59:59

# Memoire maximum utilisee par processus dans data et dans stack
# @ data_limit = 2.0Gb

# @ stack_limit = 1.0Gb,1.0Gb

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

# Pour avoir l'echo des commandes
set -x

##### nom du repertoire depuis ou est lance le code
##### et nom du sous-repertoire ou sera stockee la base de donnees
##### La variable LOADL_STEP_INITDIR est automatiquement positionnee par
##### LoadLeveler au repertoire dans lequel on tape la commande llsubmit
export repertoire_code=$( basename $LOADL_STEP_INITDIR )
export repertoire_database=DATABASES_MPI_DIMITRI

# vider les sous-repertoires dans le repertoire de depart
rm -r -f $LOADL_STEP_INITDIR/OUTPUT_FILES $LOADL_STEP_INITDIR/$repertoire_database

# copier les codes source depuis le repertoire de depart vers le repertoire temporaire
rm -r -f $TMPDIR/$repertoire_code
cp -r -p $LOADL_STEP_INITDIR $TMPDIR

# creer les nouveaux sous-repertoires dans le repertoire temporaire
mkdir $TMPDIR/$repertoire_code/OUTPUT_FILES $TMPDIR/$repertoire_code/$repertoire_database

# aller dans le repertoire temporaire
cd $TMPDIR/$repertoire_code

# compiler le mailleur et l'executer en MPI
make clean
make meshfem3D
./bin/xmeshfem3D

# compiler le solver et l'executer en MPI
make clean
make specfem3D
./bin/xspecfem3D

# deplacer les sismogrammes dans le repertoire de travail
mv $TMPDIR/$repertoire_code/$repertoire_database/*.semd $TMPDIR/$repertoire_code

# supprimer la base de donnees creee car elle est de tres grande taille
rm -r -f $TMPDIR/$repertoire_code/$repertoire_database

# recuperer le job ID
export myjobid=$( echo $LOADL_STEP_ID | cut -d'.' -f4 )

# deplacer tous les resultats dans le workdir en ajoutant le job ID au nom
# sortir d'abord dans home pour pouvoir detruire le repertoire courant de tmpdir
cd $HOME
mv $TMPDIR/$repertoire_code $WORKDIR/${repertoire_code}_${myjobid}

