<?xml version="1.0"?>

<!DOCTYPE inventory>

<inventory>

    <component name="Specfem3DGlobe">
        <facility name="model">{{ simulation.model.get_type_id }}</facility>

        <component name="solver">
            <property name="simulation-type">{{ simulation.get_simulation_type_id }}</property>
            <property name="save-forward">{{ simulation.save_forward }}</property>
            <property name="record-length">{{ simulation.record_length }}*minute</property>
            <property name="absorbing-conditions">{{ simulation.absorbing_conditions }}</property>
            <property name="movie-surface">{{ simulation.movie_surface }}</property>
            <property name="movie-volume">{{ simulation.movie_volume }}</property>
            <property name="ntstep-between-frames">{{ simulation.ntstep_between_frames }}</property>
            <property name="hdur-movie">{{ simulation.hdur_movie }}</property>
            <property name="number-of-runs">1</property>
            <property name="number-of-this-run">1</property>
            <property name="ntstep-between-output-info">{{ simulation.ntstep_between_output_info }}</property>
            <property name="ntstep-between-output-seismos">{{ simulation.ntstep_between_output_seismos }}</property>
            <property name="receivers-can-be-buried">{{ simulation.receivers_can_be_buried }}</property>
            <property name="print-source-time-function">{{ simulation.print_source_time_function }}</property>
        </component>


        <component name="mesher">
            <property name="nchunks">{{ simulation.mesh.nchunks }}</property>
            <property name="angular-width-xi">{{ simulation.mesh.angular_width_xi }}*deg</property>
            <property name="angular-width-eta">{{ simulation.mesh.angular_width_eta }}*deg</property>
            <property name="center-latitude">{{ simulation.mesh.center_latitude }}*deg</property>
            <property name="center-longitude">{{ simulation.mesh.center_longitude }}*deg</property>
            <property name="gamma-rotation-azimuth">{{ simulation.mesh.gamma_rotation_azimuth }}*deg</property>
            <property name="nex-xi">{{ simulation.mesh.nex_xi }}</property>
            <property name="nex-eta">{{ simulation.mesh.nex_eta }}</property>
            <property name="nproc-xi">{{ simulation.mesh.nproc_xi }}</property>
            <property name="nproc-eta">{{ simulation.mesh.nproc_eta }}</property>
            <property name="save-files">{{ simulation.mesh.save_files }}</property>
        </component>


        <component name="model">
            <property name="oceans">{{ simulation.model.oceans }}</property>
            <property name="ellipticity">{{ simulation.model.ellipticity }}</property>
            <property name="topography">{{ simulation.model.topography }}</property>
            <property name="gravity">{{ simulation.model.gravity }}</property>
            <property name="rotation">{{ simulation.model.rotation }}</property>
            <property name="attenuation">{{ simulation.model.attenuation }}</property>
        </component>

    </component>

</inventory>

<!-- Generated automatically by the Specfem Web Portal on {% now "F jS, Y \a\t H:i" %} -->

<!-- end of file -->
