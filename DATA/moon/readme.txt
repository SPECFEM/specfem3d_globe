---------------------------
readme
---------------------------

Moon model VPREMOON from reference:

  Garcia et al. (2011),
  Very preliminary reference Moon model,
  PEPI, 188, 96 - 113.

combined with core values from:

  Weber et al. (2011),
  Seismic Detection of the Lunar Core,
  Science, 331.

The file vpremoon.dat contains a modified table based on the model values as defined in Garcia et al. (2011),
table 6 (seismic model on the left). Since core values contain large uncertainties, their model values for the core
are left with question marks.

The modified table here in vpremoon.dat uses crustal/mantle values from Garcia et al. (2011) and
a modification for the outer fluid core and solid inner core as suggested by Weber et al. (2011).
(a solid inner core is also needed to run the current SPECFEM mesher).

The fluid outer core radius is taken from Garcia et at. (2011) which locates it at 380 +/- 40 km.
This differs from Weber et al. (2011) fluid outer core radius of 330 +/- 20 km.
The inner core radius is taken from Weber et al.'s 240 +/- 10 km.

Vp values for the fluid outer core and the inner core Vp,Vs,rho values are following Weber et al. (2011) suggestions.
Attenuation values for the cores are unknown and set at an "artificially" high Q value. Using this combination of core radius and
densities, the model overestimates the average core density and total mass of Moon - left todo for future work.

Note that Weber et al. (2011) also suggest a partial melt boundary (PMB) at radius 480 +/- 15 km,
which is not included in the VPREMOON model here.


