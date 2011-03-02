
# Choix du shell
# @ shell = /bin/ksh

# Temps CPU max. par processus MPI en hh:mm:ss
# @ cpu_limit = 00:59:59

# Mémoire max. utilisée par processus dans data et dans stack
# @ data_limit = 14.5Gb

# @ stack_limit = 1.0Gb,1.0Gb

# Nom arbitraire du travail LoadLeveler
# @ job_name = run_CPML_serial

#----------------------------------------------------

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

# vider les repertoires dans le home
rm -r -f $LOADL_STEP_INITDIR/OUTPUT_FILES $LOADL_STEP_INITDIR/$repertoire_database

# copier les codes source depuis le home vers le repertoire temporaire
rm -r -f $TMPDIR/$repertoire_code
cp -r -p $LOADL_STEP_INITDIR $TMPDIR

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
./bin/xspecfem3D

# deplacer les sismogrammes dans le repertoire de travail
mv $TMPDIR/$repertoire_code/$repertoire_database/*.semd $TMPDIR/$repertoire_code

# supprimer la base de donnees creee
rm -r -f $TMPDIR/$repertoire_code/$repertoire_database

# recuperer le job ID
export myjobid=$( echo $LOADL_STEP_ID | cut -d'.' -f4 )

# deplacer tous les resultats dans le workdir en ajoutant le job ID au nom
# sortir d'abord dans home pour pouvoir detruire le repertoire courant de tmpdir
cd $HOME
mv $TMPDIR/$repertoire_code $WORKDIR/${repertoire_code}_${myjobid}

