#include <octave/oct.h>
#include <octave/Cell.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <mc_var.h>

#include <fftw3.h>
#include <omp.h>
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

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=(double)myptr(0,0);
  return myval;

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
char *get_char_from_ov(octave_value cch)
{
  charMatrix ch=cch.char_matrix_value();

  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}

/*--------------------------------------------------------------------------------*/
int matrix_nelem(Matrix mat)
{
  dim_vector dm=mat.dims();
  return dm(0)*dm(1);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (ccamb_hello, args, nargout, "Make sure library shows up.\n")
{
  camb_hello();

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (ccamb_test, args, nargout, "Simple trial run.\n")
{
  struct Pstruct params;
  params.omega_b=0.025;
  params.omega_d=0.11;
  params.Omega_L=0.7;
  params.Omega_k=0;

  params.tauz=0.08;
  params.n_ad=0.98;

  params.A_s=2.3;
  params.a_iso=0;
  params.b_iso=0;
  params.n_iso=0;
  params.A_tens=0;
  params.n_run=0;
  params.sz_amp=0;
  params.w_lam=0;
  params.omega_n=0;
  params.N_nu=0;
  params.z_r=10;
  params.del_z=1;
  params.w_a=0;
  params.y_he=0.24;
  params.A_p=0;
  params.A_c=0;


  int nell=(int)get_value(args(0));
  Matrix cl(nell,4);
  
  int nk=get_nk_max();
  Matrix pk(nk,1);
  Matrix kvec(nk,1);
  Matrix out_vec(22,1);
  Matrix tscale(4,1);
  //dim_vector dm(nk,4,3);
  //NDArray lrg_vec(dm);
  Matrix lrg_vec(nk,4*3);

  double *vv=cl.fortran_vec();
  int nl_flag=2;
  int reion_flag=0;
  get_spectra((double *)&params,nl_flag,reion_flag,nell,vv,vv+nell,vv+2*nell,vv+3*nell,pk.fortran_vec(),kvec.fortran_vec(),out_vec.fortran_vec(),tscale.fortran_vec(),lrg_vec.fortran_vec());

  printf("Back from get_spectra.\n");

  return octave_value(cl);
  
}


/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_nk_max, args, nargout, "Find nkmax from camb.\n")
{
  int nk_max=get_nk_max();
  return octave_value(nk_max);
}
