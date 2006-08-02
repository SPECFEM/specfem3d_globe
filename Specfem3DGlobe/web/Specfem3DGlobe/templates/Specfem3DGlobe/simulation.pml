<?xml version="1.0"?>
<!--
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
! {LicenseText}
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-->

<!DOCTYPE inventory>

<inventory>

    <component name="Specfem3DGlobe">
        <facility name="model">isotropic_prem</facility>
        <property name="nodes">24</property>
        <property name="scratch-dir">/scratch/${user}/${job.name}-${job.id}-${rank}</property>

        <component name="lsf">
            <property name="bsub-options">[-a mpich_gm]</property>
            <property name="wait">true</property>
        </component>


        <component name="solver">
            <property name="number-of-runs">1</property>
            <property name="record-length">2.0*minute</property>
            <property name="save-forward">False</property>
            <property name="ntstep-between-output-seismos">5000000</property>
            <property name="stations">STATIONS</property>
            <property name="absorbing-conditions">False</property>
            <property name="output-file">OUTPUT_FILES/output_solver.txt</property>
            <property name="movie-volume">False</property>
            <property name="header-file">OUTPUT_FILES/values_from_mesher.h</property>
            <property name="number-of-this-run">1</property>
            <property name="receivers-can-be-buried">False</property>
            <property name="print-source-time-function">False</property>
            <property name="cmt-solution">CMTSOLUTION</property>
            <property name="seismogram-archive">OUTPUT_FILES/seismograms.tar.gz</property>

            <!-- here is an example -->
            <property name="ntstep-between-frames">{{ simulation.ntstep_between_frames }}</property>

            <property name="ntstep-between-output-info">200</property>
            <property name="hdur-movie">0.0</property>
            <property name="simulation-type">forward</property>
            <property name="movie-surface">False</property>
        </component>


        <component name="launcher">
            <property name="command">mpirun.lsf --gm-no-shmem --gm-copy-env --gm-recv hybrid --gm-kill 10</property>
        </component>


        <component name="job">
            <property name="script">OUTPUT_FILES/batch-script.py</property>
            <property name="stdout">OUTPUT_FILES/stdout.txt</property>
            <property name="queue">normal</property>
            <property name="stderr">OUTPUT_FILES/stderr.txt</property>
            <property name="interpreter">OUTPUT_FILES/pyspecfem3D</property>
            <property name="mpiScript">OUTPUT_FILES/mpi-script.py</property>
            <property name="walltime">30*minute</property>
        </component>


        <component name="mesher">
            <property name="nex-eta">64</property>
            <property name="nchunks">6</property>
            <property name="nex-xi">64</property>
            <property name="output-file">OUTPUT_FILES/output_mesher.txt</property>
            <property name="center-longitude">10.0*deg</property>
            <property name="center-latitude">40.0*deg</property>
            <property name="angular-width-xi">90.0*deg</property>
            <property name="nproc-eta">2</property>
            <property name="gamma-rotation-azimuth">20.0*deg</property>
            <property name="nproc-xi">2</property>
            <property name="angular-width-eta">90.0*deg</property>
            <property name="save-files">False</property>
        </component>


        <component name="model">
            <property name="topography">False</property>
            <property name="oceans">False</property>
            <property name="gravity">False</property>
            <property name="attenuation">False</property>
            <property name="ellipticity">False</property>
            <property name="rotation">False</property>
        </component>

    </component>

</inventory>

<!-- version-->
<!-- $Id$-->

<!-- Generated automatically by Renderer on Tue Aug  1 12:41:27 2006-->

<!-- End of file -->
