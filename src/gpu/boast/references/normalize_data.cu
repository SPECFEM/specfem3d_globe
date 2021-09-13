// original implementation

#ifndef NGLL3
#define NGLL3 125
#endif


__global__ void normalize_data(realw_p data_smooth,
                               realw_const_p normalisation,
                               int nker,
                               int nspec_me){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;

  // note: normalization coefficient is added nker times, thus divide by nker
  realw norm = normalisation[NGLL3*ispec + igll] / nker;

  // avoids division by zero
  if (norm < 1.e-24) norm = 1.0f;

  // normalizes smoothed kernel values
  for (int iker=0; iker<nker; iker++) data_smooth[NGLL3*nspec_me*iker + NGLL3*ispec + igll] /= norm;
}

