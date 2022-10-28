------------------------------

Moho models for Mars & Moon

------------------------------

The crustal Moho models here are calculated with Mark's ctplanet tool:

  see: https://github.com/MarkWieczorek/ctplanet

    > git clone https://github.com/MarkWieczorek/ctplanet.git
    > cd ctplanet/
    > pip install .
    > cd examples/
    > mkdir -p InSight; mkdir -p InSight/constant

  as an example for Mars, we use as a map generated with an input Moho depth estimate of 20km at the lander site.
  this corresponds to a thin crust, plan A in fig. 4 of Knapmeyer-Endrun et al. (2021, Science, 373).

  no porosity layer assumed, no density change across the dichotomy boundary.

    > python Mars-Crust-InSight.py <<EOF
      20
      0
      EOF

   this will generate various crustal maps in InSight/constant/ for different models (DWThot, EH45Tcold, LFAK,..).
   the crustal features look all very similar, main difference is the estimated Moho depth ranges
   due to different average density values assumed.

   -> we will take the model DWAK-2550 (average crustal density 2550) with a depth range ~ 2-62 km

      reason: it has a slightly thicker minimum crust and larger dynamic range than others, but no real preference here.

              the minimum crustal thickness is around 2km, which makes mesh elements very thin when we honor the Moho depths.
              these small elements will then determine the largest possible time step size for the simulations,
              making these simulations more costly than with thicker crustal elements.
              however, we could/should test all appropriate models.

   references:
     M. A. Wieczorek, A.-C. Plesa, B. Knapmeyer-Endrun, S. M. McLennan, F. Nimmo, C. Michaut,
     A. Broquet, S. Smrekar, and W. B. Banerdt (2020),
     GLOBAL CRUSTAL THICKNESS MODELING OF MARS USING INSIGHT SEISMIC CONSTRAINTS,
     51st Lunar and Planetary Science Conference (2020)

     Wieczorek, M. A., G. A. Neumann, F. Nimmo, W. S. Kiefer, G. J. Taylor, H. J. Melosh, R. J. Phillips,
     S. C. Solomon, J. C. Andrews-Hanna, S. W. Asmar, A. S. Konopliv, F. G. Lemoine, D. E. Smith,
     M. M. Watkins, J. G. Williams, M. T. Zuber (2013),
     The crust of the Moon as seen by GRAIL,
     Science, 339, 671-675, doi:10.1126/science.1231530.

     Wieczorek, M. A., M. Beuthe, A. Rivoldini, and T. Van Hoolst (2019),
     Hydrostatic interfaces in bodies with nonhydrostatic lithospheres,
     Journal of Geophysical Research: Planets, 124, doi:10.1029/2018JE005909.


In this folder, the following models are provided for testing:

* Mars:

  the thickness choices are motivated by the possible moho models suggested in Knapmeyer-Endrun et al. (2021), see figure 4:
    plan A: moho at 20km, 2-layer, upper/lower mid-crust boundary at 9km
    plan B: moho at 39km, 3-layer, upper/middle boundary at 9km, middle/lower boundary at 24km

  - Moho-Mars-DWAK-20-2550.sh      : mars model DWAK and a crustal thickness of 20 km at Lander site
  - Moho-Mars-DWAK-39-2550.sh      : mars model DWAK and a crustal thickness of 39 km at Lander site
  - Moho-Mars-DWAK-39-2550-2650.sh : mars model DWAK w/ dichotomy boundary and a crustal thickness of 39 km at Lander site

* Moon:
  - Moho-Moon-1-Constant-Density.sh : moon constant density model
  - Moho-Moon-2-Variable-Density.sh : moon variable density model

  note: Moon output files were created by a modified script to limit the output to a spherical harmonic degree lmax = 90,
        instead of lmax = 900 provided by Mark's original routine. The resolution is therefore similar to the Mars models above.

        Using the higher lmax resolution would be possible, but the creation of crustmaps becomes very slow. Furthermore, the crustmaps are
        lower resolution (1/6 degree), and such high resolutions maps are likely overconfident with respect to the uncertainties in such crustal models.


