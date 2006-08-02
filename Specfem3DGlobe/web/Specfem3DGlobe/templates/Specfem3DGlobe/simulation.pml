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
        <facility name="model">{{ simulation.model.get_type_display }}</facility>
        <property name="nodes">24</property>
        <property name="scratch-dir">/scratch/${user}/${job.name}-${job.id}-${rank}</property>

        <component name="lsf">
            <property name="bsub-options">[-a mpich_gm]</property>
            <property name="wait">true</property>
        </component>


        <component name="solver">
            <property name="number-of-runs">1</property>
            <property name="record-length">{{ simulation.record_length }}*minute</property>
            <property name="save-forward">{{ simulation.save_forward }}</property>
            <property name="ntstep-between-output-seismos">{{ simulation.ntstep_between_output_seismos }}</property>
            <property name="stations">STATIONS</property>
            <property name="absorbing-conditions">{{ simulation.absorbing_conditions }}</property>
            <property name="output-file">OUTPUT_FILES/output_solver.txt</property>
            <property name="movie-volume">{{ simulation.movie_volume }}</property>
            <property name="header-file">OUTPUT_FILES/values_from_mesher.h</property>
            <property name="number-of-this-run">1</property>
            <property name="receivers-can-be-buried">{{ simulation.receivers_can_be_buried }}</property>
            <property name="print-source-time-function">{{ simulation.print_source_time_function }}</property>
            <property name="cmt-solution">CMTSOLUTION</property>
            <property name="seismogram-archive">OUTPUT_FILES/seismograms.tar.gz</property>

            <!-- here is an example -->
            <property name="ntstep-between-frames">{{ simulation.ntstep_between_frames }}</property>

            <property name="ntstep-between-output-info">{{ simulation.ntstep_between_output_info }}</property>
            <property name="hdur-movie">{{ simulation.hdur_movie }}</property>
            <property name="simulation-type">{{ simulation.model.get_type_display }}</property>
            <property name="movie-surface">{{ simulation.movie_surface }}</property>
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
            <property name="nex-eta">{{ simulation.mesh.nex_eta }}</property>
            <property name="nchunks">{{ simulation.mesh.nchunks }}</property>
            <property name="nex-xi">{{ simulation.mesh.nex_xi }}</property>
            <property name="output-file">OUTPUT_FILES/output_mesher.txt</property>
            <property name="center-longitude">{{ simulation.mesh.center_longitude }}*deg</property>
            <property name="center-latitude">{{ simulation.mesh.center_latitude }}*deg</property>
            <property name="angular-width-xi">{{ simulation.mesh.angular_width_xi }}*deg</property>
            <property name="nproc-eta">{{ simulation.mesh.nproc_eta }}</property>
            <property name="gamma-rotation-azimuth">{{ simulation.mesh.gamma_rotation_azimuth }}*deg</property>
            <property name="nproc-xi">{{ simulation.mesh.nproc_xi }}</property>
            <property name="angular-width-eta">{{ simulation.mesh.angular_width_eta }}*deg</property>
            <property name="save-files">{{ simulation.mesh.save_files }}</property>
        </component>


        <component name="model">
            <property name="topography">{{ simulation.model.topography }}</property>
            <property name="oceans">{{ simulation.model.oceans }}</property>
            <property name="gravity">{{ simulation.model.gravity }}</property>
            <property name="attenuation">{{ simulation.model.attenuation }}</property>
            <property name="ellipticity">{{ simulation.model.ellipticity }}</property>
            <property name="rotation">{{ simulation.model.rotation }}</property>
        </component>

    </component>

</inventory>

<!-- version-->
<!-- $Id$-->

<!-- Generated automatically by Renderer on Tue Aug  1 12:41:27 2006-->

<!-- End of file -->
