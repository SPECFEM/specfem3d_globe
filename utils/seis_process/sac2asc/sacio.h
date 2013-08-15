/*******************************************************************
*     sacio.c
* SAC I/O functions:
* read_sachead  read SAC header
* read_sac  read SAC binary data
* write_sac write SAC binary data
* wrtsac0   write 1D array as evenly-spaced SAC binary data
* wrtsac0_  fortran wrap for wrtsac0
* wrtsac2   write 2 1D arrays as XY SAC data
* wrtsac2_  fortran wrap for wrtsac2
* swab4   reverse byte order for integer/float
*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"

/***********************************************************

  read_sachead

  Description:  read binary SAC header from file.

  Author: Lupei Zhu

  Arguments:  const char *name  file name
    SACHEAD *hd   SAC header to be filled

  Return: 0 if success, -1 if failed

  Modify history:
  05/29/97  Lupei Zhu Initial coding
************************************************************/

int read_sachead(const char *name,
    SACHEAD   *hd,
    int   swap_bytes
  )
{
  FILE    *strm;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return -1;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     fclose(strm);
     return -1;
  }

if(swap_bytes)
  swab4((char *) hd, HD_SIZE);

  fclose(strm);
  return 0;

}


/***********************************************************

  read_sac

  Description:  read binary data from file. If succeed, it will return
    a float pointer to the read data array. The SAC header
    is also filled. A NULL pointer is returned if failed.

  Author: Lupei Zhu

  Arguments:  const char *name  file name
    SACHEAD *hd   SAC header to be filled

  Return: float pointer to the data array, NULL if failed

  Modify history:
  09/20/93  Lupei Zhu Initial coding
  12/05/96  Lupei Zhu adding error handling
  12/06/96  Lupei Zhu swap byte-order on PC
************************************************************/

float*  read_sac(const char *name,
    SACHEAD   *hd,
    int   swap_bytes
  )
{
  FILE    *strm;
  float   *ar;
  unsigned  sz;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

if(swap_bytes)
  swab4((char *) hd, HD_SIZE);

if(hd->npts < 0 || hd->npts > 1e+08 || hd->delta < 1.0e-15)
   {
   fprintf(stderr,"%d %e\n",hd->npts,hd->delta);

   if(swap_bytes == 1)
      swap_bytes = 0;
   else
      swap_bytes = 1;

   swab4((char *) hd, HD_SIZE);
   }

  sz = hd->npts*sizeof(float);
  if ((ar = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }

  if (fread((char *) ar, sz, 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }

  fclose(strm);

if(swap_bytes)
  swab4((char *) ar, sz);

  return ar;

}



/***********************************************************

  write_sac

  Description:  write SAC binary data.

  Author: Lupei Zhu

  Arguments:  const char *name  file name
    SACHEAD hd    SAC header
    const float *ar   pointer to the data

  Return: 0 if succeed; -1 if failed

  Modify history:
  09/20/93  Lupei Zhu Initial coding
  12/05/96  Lupei Zhu adding error handling
  12/06/96  Lupei Zhu swap byte-order on PC
************************************************************/

int write_sac(const char  *name,
    SACHEAD   hd,
    const float *ar,
    int   swap_bytes
  )
{
  FILE    *strm;
  unsigned  sz;
  float   *data;
  int   error = 0;

  sz = hd.npts*sizeof(float);
  if (hd.iftype == IXY) sz *= 2;

  if ((data = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for writing %s\n",name);
     error = 1;
  }

  if ( !error && memcpy(data, ar, sz) == NULL) {
     fprintf(stderr, "Error in copying data for writing %s\n",name);
     error = 1;
  }

if(swap_bytes)
  {
  swab4((char *) data, sz);
  swab4((char *) &hd, HD_SIZE);
  }

  if ( !error && (strm = fopen(name, "w")) == NULL ) {
     fprintf(stderr,"Error in opening file for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(&hd, sizeof(SACHEAD), 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC header for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(data, sz, 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC data for writing %s\n",name);
     error = 1;
  }

  free(data);
  fclose(strm);

  return (error==0) ? 0 : -1;

}



/*****************************************************

  swab4

  Description:  reverse byte order for float/integer

  Author: Lupei Zhu

  Arguments:  char *pt  pointer to byte array
    int    n  number of bytes

  Return: none

  Modify history:
  12/03/96  Lupei Zhu Initial coding

************************************************************/

void  swab4(  char  *pt,
    int n
  )
{
  int i;
  char temp;
  for(i=0;i<n;i+=4) {
    temp = pt[i+3];
    pt[i+3] = pt[i];
    pt[i] = temp;
    temp = pt[i+2];
    pt[i+2] = pt[i+1];
    pt[i+1] = temp;
  }
}



/***********************************************************

  wrtsac0

  Description:  write 1D array into evenly spaced SAC data.

  Author: Lupei Zhu

  Arguments:  const char *name  file name
    float dt    sampling interval
    int ns    number of points
    float b0    starting time
    float dist    distance range
    const float *ar   data array

  Return: 0 if succeed; -1 if failed

  Modify history:
  09/20/93  Lupei Zhu Initial coding
  12/05/96  Lupei Zhu adding error handling
  12/06/96  Lupei Zhu swab byte-order on PC
************************************************************/

int    wrtsac0(const char *name,
    float   dt,
    int   ns,
    float   b0,
    float   dist,
    const float *ar,
    int swap_bytes
  )
{

  SACHEAD hd = sac_null;

  hd.npts = ns;
  hd.delta = dt;
  hd.dist = dist;
  hd.b = b0;
  hd.o = 0.;
  hd.e = b0+(hd.npts-1)*hd.delta;
  hd.iztype = IO;
  hd.iftype = ITIME;
  hd.leven = TRUE;

  return write_sac(name, hd, ar,swap_bytes);

}



/***********************************************************

  wrtsac2

  Description:  write 2 arrays into XY SAC data.

  Author: Lupei Zhu

  Arguments:  const char *name  file name
    int ns    number of points
    const float *x    x data array
    const float *y    y data array

  Return: 0 succeed, -1 fail

  Modify history:
  09/20/93  Lupei Zhu Initial coding
  12/05/96  Lupei Zhu adding error handling
  12/06/96  Lupei Zhu swap byte-order on PC
************************************************************/

int wrtsac2(const char  *name,
    int   n,
    const float *x,
    const float *y,
    int   swap_bytes
  )
{
  SACHEAD hd = sac_null;
  float   *ar;
  unsigned  sz;
  int   exit_code;

  hd.npts = n;
  hd.iftype = IXY;
  hd.leven = FALSE;

  sz = n*sizeof(float);

  if ( (ar = (float *) malloc(2*sz)) == NULL ) {
     fprintf(stderr, "error in allocating memory%s\n",name);
     return -1;
  }

  if (memcpy(ar, x, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }
  if (memcpy(ar+sz, y, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }

  exit_code = write_sac(name, hd, ar,swap_bytes);

  free(ar);

  return exit_code;

}


/* write user data t[n] into sac header starting at pos (0-9)*/
void sacUdata(float *hd, int pos, float *t, int n, int type)
{
  int i;
  hd += type;
  for(i=pos;i<n && i<10;i++) hd[i] = t[i];
}


/* for fortran--write evenly-spaced data */
void    wrtsac0_(const char *name, float dt, int ns, float b0, float dist, const float *ar) {
  wrtsac0(name, dt, ns, b0, dist, ar,0);
}

/* for fortran--write x-y data */
void    wrtsac2_(const char *name, int n, const float *x, const float *y) {
  wrtsac2(name, n, x, y,0);
}

