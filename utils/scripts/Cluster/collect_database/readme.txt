-------------------------------
collect files on a cluster
-------------------------------

the scripts here help to copy files from remote compute nodes of a cluster to a local (master) machine.

- setup:
the scripts assume that when running the mesher/solver on the cluster, output files will be stored on local
scratch filesystems on each compute node:
/scratch/$ENV{USER}/DATABASES_MPI/
or
/scratch/$ENV{USER}/DATABASES_MPI$jobid/

all the hostnames of the necessary nodes need to be listed in a `machinefile` or `lsf_machine_file`.
The `machinefile` lists hostnames of the compute nodes, the `lsf_machine_file` is created by the LSF scheduler.
The scripts using `lsf_machine_file` require a format like:
  ip-addressN number-of-processes ip-addressN-1 number-of-processes ..
  127.0.0.5 1 127.0.0.4 1 127.0.0.3 1 127.0.0.2 1 127.0.0.1 1
assuming that process 0 runs on node listed last and process N runs on first node listed.

- requirements:
the scripts need to have `scp` for the file transfer to the local directory ./

- author:
scripts written by Qinya Liu, May 2007, Caltech



