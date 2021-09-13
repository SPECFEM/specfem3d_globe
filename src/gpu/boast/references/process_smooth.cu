// original implementation

#ifndef NGLL3
#define NGLL3 125
#endif

__global__ void process_smooth(realw_const_p xstore_me,
                               realw_const_p ystore_me,
                               realw_const_p zstore_me,
                               realw_const_p xstore_other,
                               realw_const_p ystore_other,
                               realw_const_p zstore_other,
                               realw_const_p data_other,
                               const realw sigma_h2_inv,
                               const realw sigma_v2_inv,
                               const int iker,
                               const int nspec_me,
                               const int nspec_other,
                               const realw v_criterion,
                               const realw h_criterion,
                               realw_const_p integ_factor,
                               realw_p data_smooth,
                               realw_p normalisation){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;

  int gll_other;
  realw x_me,y_me,z_me;
  realw x_other,y_other,z_other;
  realw center_x,center_y,center_z;
  realw alpha,ratio,theta;
  realw r0_squared,r1_squared;
  realw dist_h,dist_v;
  realw val,val_gaussian;
  realw coef, normalisation_slice;
  realw dat;

  // PI squared
  const realw PI2 = 9.869604401089358f;

  __shared__ int sh_test[NGLL3];
  __shared__ realw sh_x_other[NGLL3];
  __shared__ realw sh_y_other[NGLL3];
  __shared__ realw sh_z_other[NGLL3];
  __shared__ realw sh_integ_factor[NGLL3];
  __shared__ realw sh_data[NGLL3];

  // for each reference GLL point, we can check a block of 125 neighbor elements
  // by convenience, the block size is set to the number of threads 125 of this kernel
  int n_loop = nspec_other/NGLL3 + 1;

  // reference GLL point position
  x_me = xstore_me[NGLL3*ispec + igll];
  y_me = ystore_me[NGLL3*ispec + igll];
  z_me = zstore_me[NGLL3*ispec + igll];

  __syncthreads();

  dat = 0.f;
  normalisation_slice = 0.f;

  // We test 125 spectral elements at a time
  for (int i=0; i < n_loop; i++){
    __syncthreads();

    // each thread helps to test a different element in the other slice (using the center position)
    // number of threads == NGLL3 == 125
    // for i==0: element range [0,124]
    // for i==1: element range [125,(125+124)]
    // ..
    // for i==n_loop-1: element range [NGLL3*(nloop-1),NGLL3*(nloop-1)+124]
    //                  where NGLL3*(nloop-1)+124 is equal to nspec_other (or slightly greater)
    int ispec_other = NGLL3*i + igll;

    if (ispec_other < nspec_other){
      // center position
      center_x = (xstore_other[ispec_other * NGLL3] + xstore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;
      center_y = (ystore_other[ispec_other * NGLL3] + ystore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;
      center_z = (zstore_other[ispec_other * NGLL3] + zstore_other[ispec_other * NGLL3 + (NGLL3 - 1)]) * 0.5f;

      // note: instead of distance we use distance squared to avoid too many sqrt() operations

      // Cartesian case
      // distance horizontal = (x-x0)**2 + (y-y0)**2, and vertical = (z-z0)**2
      //dist_h = (x_me - center_x)*(x_me - center_x) + (y_me - center_y)*(y_me - center_y);
      //dist_v = (z_me - center_z)*(z_me - center_z);

      // Spherical case
      // vertical distance
      r0_squared = x_me*x_me + y_me*y_me + z_me*z_me;
      r1_squared = center_x*center_x + center_y*center_y + center_z*center_z;

      // vertical distance (squared)
      // dist_v = (r1 - r0)*(r1 - r0)
      //        = r1**2 + r0**2 - 2 * alpha
      //          with alpha = sqrt( r0**2 * r1**2 ) = r0 * r1
      // this avoids using sqrt() function too often which is costly
      alpha = sqrt( r0_squared * r1_squared );
      dist_v = r1_squared + r0_squared - 2.0f * alpha;

      // epicentral distance
      // (accounting for spherical curvature)
      // calculates distance of circular segment
      // angle between r0 and r1 in radian
      // given by dot-product of two vectors
      if (alpha > 0.0f){
        ratio = (x_me*center_x + y_me*center_y + z_me*center_z) / alpha;
      } else {
        ratio = 1.0f;
      }

      // checks boundaries of ratio (due to numerical inaccuracies)
      if (ratio >= 1.0f){
        // ratio = 1.0_CUSTOM_REAL
        // -> acos(1) = 0
        // -> dist_h = 0
        dist_h = 0.0f;
      } else if (ratio <= -1.0f) {
        // ratio = -1.0_CUSTOM_REAL
        // -> acos(-1) = PI
        // -> dist_h = r1**2 * PI**2
        dist_h = r1_squared * PI2;
      } else {
        theta = acos( ratio );
        // segment length at heigth of r1 (squared)
        dist_h = r1_squared * (theta*theta);
      }
    } else {
      // artificial high values
      // (h_criterion and v_criterion are normalized in global version)
      dist_v = 99999999.f;
      dist_h = 99999999.f;
    }

    // tests if element is too far away
    sh_test[igll] = ( ispec_other >= nspec_other
                    || dist_h > h_criterion
                    || dist_v > v_criterion ) ? 1 : 0 ;


    __syncthreads();

    // loops over each spectral element tested
    for (int k=0; k < NGLL3; k++){
      __syncthreads();

      // skips element if test was true (too far away)
      if (sh_test[k]) continue ;

      // loads data from other slice to shared memory
      int ispec_test = i*NGLL3 + k;
      sh_x_other[igll] = xstore_other[ispec_test*NGLL3 + igll];
      sh_y_other[igll] = ystore_other[ispec_test*NGLL3 + igll];
      sh_z_other[igll] = zstore_other[ispec_test*NGLL3 + igll];

      sh_data[igll] = data_other[ispec_test*NGLL3 + igll];
      sh_integ_factor[igll] = integ_factor[ispec_test*NGLL3 + igll];

      __syncthreads();

      // loops over gll points
      for (int j=0; j < NGLL3; j++){
        gll_other = (igll + j) % NGLL3;

        x_other = sh_x_other[gll_other];
        y_other = sh_y_other[gll_other];
        z_other = sh_z_other[gll_other];

        // Cartesian case
        // distance horizontal = (x-x0)**2 + (y-y0)**2, and vertical = (z-z0)**2
        //dist_h = (x_me - x_other)*(x_me - x_other) + (y_me - y_other)*(y_me - y_other);
        //dist_v = (z_me - z_other)*(z_me - z_other);
        //coef = expf(- sigma_h2_inv * dist_h - sigma_v2_inv * dist_v) * sh_integ_factor[gll_other];

        // Spherical case
        // vertical distance
        r0_squared = x_me*x_me + y_me*y_me + z_me*z_me;
        r1_squared = x_other*x_other + y_other*y_other + z_other*z_other;

        // vertical distance (squared)
        // dist_v = (r1 - r0)*(r1 - r0)
        //        = r1**2 + r0**2 - 2 * alpha
        //          with alpha = sqrt( r0**2 * r1**2 ) = r0 * r1
        // this avoids using sqrt() function too often which is costly
        alpha = sqrt( r0_squared * r1_squared );
        dist_v = r1_squared + r0_squared - 2.0f * alpha;

        // epicentral distance
        // (accounting for spherical curvature)
        // calculates distance of circular segment
        // angle between r0 and r1 in radian
        // given by dot-product of two vectors
        if (alpha > 0.0f){
          ratio = (x_me*x_other + y_me*y_other + z_me*z_other) / alpha;
        } else {
          ratio = 1.0f;
        }

        // checks boundaries of ratio (due to numerical inaccuracies)
        if (ratio >= 1.0f){
          // ratio = 1.0_CUSTOM_REAL
          // -> acos(1) = 0
          // -> dist_h = 0
          dist_h = 0.0f;
        } else if (ratio <= -1.0f) {
          // ratio = -1.0_CUSTOM_REAL
          // -> acos(-1) = PI
          // -> dist_h = r1**2 * PI**2
          dist_h = r1_squared * PI2;
        } else {
          theta = acos( ratio );
          // segment length at heigth of r1 (squared)
          dist_h = r1_squared * (theta*theta);
        }

        // Gaussian function
        val = - dist_h*sigma_h2_inv - dist_v*sigma_v2_inv;

        // limits to single precision
        if (val < - 86.0f){
          // smaller than numerical precision: exp(-86) < 1.e-37
          val_gaussian = 0.0f;
        } else {
          val_gaussian = expf(val);
        }

        coef = val_gaussian * sh_integ_factor[gll_other];

        normalisation_slice = normalisation_slice + coef;
        dat += sh_data[gll_other] * coef;
      } //loop on each gll_other
    } //loop on each spec_other tested
  } //loop on each serie of 125 spec_other

  data_smooth[NGLL3*nspec_me*iker + NGLL3*ispec + igll] += dat;

  // note: normalization coefficient is added nker times
  normalisation[NGLL3*ispec + igll] += normalisation_slice;
}

