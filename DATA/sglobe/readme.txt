------------------------------------------------
SGLOBE-rani
------------------------------------------------

3D global shear-wave isotropic and radially anisotropic model from a joint inversion for multiple data sets

SGLOBE-rani is a global radially anisotropic shear wave speed model with radial anisotropy allowed in the whole mantle.
It is based on a seismic data set of over 43M seismic surface wave (fundamental and overtones) and body wave measurements.
It models simultaneously crustal thickness and mantle structure, and its reference model is PREM.

NOTE: Kustowski et al., 2008 and Chang et al., 2014 showed that the anisotropic structure in the lowermost mantle
      retrieved from global tomographic inversions can be strongly affected by leakage effects, so we discourage
      interpreting SGLOBE-rani's anisotropic structure below ~1500 km depth.

reference:
  Chang, S.-J., A.M.G. Ferreira, J. Ritsema, H.J. van Heijst, and J.H. Woodhouse (2015),
  Joint inversion for global isotropic and radially anisotropic mantle structure including crustal thickness perturbations,
  J. Geophys. Res., 120, 4278-4300, https://doi.org/10.1002/2014JB011824.

implementation:
  The model parameterization uses spherical harmonics up to degree 35 (horizontally) and 20 vertical splines.

  P-wave velocity perturbations (dvp) taken from P12 of S20RTS,
  density perturbations (drho) scaled from Vsv perturbations (dvsv);

  The mantle model is defined between Moho and CMB;
  uses PREM as 1D reference (also for attenuation & eta-parameter)

  The isotropic version takes the Voigt average for Vs and Vp from vsv & vsh and vpv & vph, respectively.


