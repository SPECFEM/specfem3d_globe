#!/bin/csh

# suppress files to leave regular version only

echo " "
echo "Are you sure you want to permanently remove the inlined version (y/n)?"
echo " "
set answer = $<

if ($answer == "y" || $answer == "Y" || $answer == "yes" || $answer == "YES" || $answer == "Yes") then

rm -f compute_forces_crust_mantle_inlined_5.f90 compute_forces_inner_core_inlined_5.f90 compute_forces_outer_core_inlined_5.f90
rm -f specfem3D_inlined.f90

mv -f compute_forces_crust_mantle_regular.f90 compute_forces_crust_mantle.f90
mv -f compute_forces_inner_core_regular.f90 compute_forces_inner_core.f90
mv -f compute_forces_outer_core_regular.f90 compute_forces_outer_core.f90
mv -f specfem3D_regular.f90 specfem3D.f90

rm -f copy_files_inlined_5.csh copy_files_regular.csh remove_classic.csh script_job_ES.csh 

rm -f -r DATA/util/inline_code DATA/util/scaling_nodes_ES DATA/util/events_ES_8.5sec

# suppress comment about inlined version in README file
cat README_SPECFEM3D_V3.3 | grep -v " NOTE " > _____tutu_____
mv -f _____tutu_____ README_SPECFEM3D_V3.3

  echo done permanently removing the inlined version

else

  echo okay, nothing has been removed

endif

