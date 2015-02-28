#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <ninkasi_config.h>
#include <ninkasi.h>
#include <ninkasi_mathutils.h>
#include <noise.h>
  //#include <dirfile.h>
#include <readtod.h>
#include <math.h>
#ifdef _MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif


#include <fftw3.h>
#include <assert.h>


#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;

/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

}
/*--------------------------------------------------------------------------------*/
char *get_char_from_arg(charMatrix ch)
{
  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}


/*--------------------------------------------------------------------------------*/

actData get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  actData myval=(actData)myptr(0,0);
  return myval;

}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (fit_tod_noise_c, args, nargout, "Fit TOD noise.  If rotmat is non-null, rotate by it.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  assert(mytod->have_data);
  if (mytod->rotmat)
    rotate_data(mytod,'n',NULL);
  
  mytod->noise=nkFitTODNoise(mytod, MBNOISE_LINEAR_POWLAW,0.5,80,-1.5);
  for (int i=0;i<mytod->ndet;i++) {
    if (1)  //patch up bad fits here.
      ;
  }
  

  if (mytod->rotmat)
    rotate_data(mytod,'t',NULL);

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (filter_tod_noise_c, args, nargout, "Apply noise filters to TOD.  If rotmat is non-null, rotate by it.\n")
  
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  assert(mytod->have_data);
  if (mytod->rotmat)
    rotate_data(mytod,'n',NULL);
  
  filter_data_wnoise(mytod);
  
  if (mytod->rotmat)
    rotate_data(mytod,'t',NULL);
  
  return octave_value_list();
}

