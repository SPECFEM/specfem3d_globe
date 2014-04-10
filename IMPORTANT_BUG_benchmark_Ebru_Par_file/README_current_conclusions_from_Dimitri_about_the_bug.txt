The three input files needed: Par_file, CMTSOLUTION and STATIONS are in this directory. (and should be moved to directory "DATA" before starting a run)

====================

runs are OK when attenuation is turned off in the Par_file; thus the bug is related at least in part to attenuation

the problem (very likely) seems to be in the mesher because when not using mcmodel=medium to compile the mesher only and create the mesh, the results are then unchanged whether we compile the solver with mcmodel=medium or not; however, if we compile the mesher with mcmodel=medium instead of without, then the seismograms are (very) different

attenuation with a 1D Earth model is OK (GJI 2002 attenuation benchmark is OK);
note from Dimitri: also confirmed by this comment from last August "All 1D examples we have tested are fine (i.e. give identical seismograms with and without -mcmodel=medium -shared-intel), thus the bug is likely in the 3D routines."

thus the bug seems to be in the mesher, and come from attenuation in 3D models (or at least in the 3D model used by Ebru, which is s362ani); a good test could be to try another 3D model (e.g. s40rts) with attenuation, to see if the bug is still there. I think I remember from last August that the answer is yes (the bug is still there) but I am only 90% sure because it was several months ago; we should redo the test

the bug does not come from option -ftz (Flush to Zero) used to compile the code, which in principle can slightly change roundoff results, because the seismograms are unchanged whether we use it or not

the bug is probably not in the MPI routines, since the 1D attenuation benchmark is OK (and the MPI routines do not know if the model is 1D or 3D), and last summer Daniel rewrote the MPI routines almost entirely but the bug is still there

the bug is not something obvious such as an array index going out of bounds, because the code does not display any error message for that benchmark when running with full range checking options (as well as stack checking, unassigned variable detection and so on)

the bug is not related to the new UNDO_ATTENUATION option, because it is already present in older versions (UNDO_ATTENUATION was added only recently), and in addition UNDO_ATTENUATION is off in the Par_file of this test

Ebru's current runs are highly affected (the example is her Par_file, with a single change: 60 min instead of 100 min for the total duration of the seismograms, just to reduce the cost of the tests)

The bug has been there for at least a year, probably more (my comment about this in the to-do list was from August 2013, back then I thought the switch to Daniel's new version at the end of Aug 2013 would fix it, but unfortunately it is still there); thus it is likely in a part of the code that did not change much when we switched to the new version

The bug is weird and leads to "random" behavior because for instance gfortran v4.6 does not give the same seismograms as gfortran v4.8, even without mcmodel=medium, and the difference is VERY big

However, it is reproducible: running the same weird case twice seems to always give the same seismograms.

All stations are affected (some more than others), thus I only use and plot the first one (station "AAE").

to all: please confirm the random behavior, but the machine I used is very reliable (European supercomputing center) and the 1D results with attenuation are always OK and fit the normal-mode exact solution perfectly, thus I am 99.99% sure it is a bug in the code rather than a problem on that machine. However, let us confirm that.

