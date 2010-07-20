#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>


DEFUN_DLD (double_vec, args, nargout,
           "Doubling a vector.\n")
{
  int nargin = args.length ();
  if (nargin != 1)
    print_usage ();
  NDArray A = args(0).array_value ();

  printf("have %d elements.\n",A.nelem());
  printf("First element is %12.4g.\n",A.elem(0,0));

  dim_vector dd=A.dims();
  int nx=dd.elem(0);
  int ny=dd.elem(1);

#if 0
  double *vec=A.fortran_vec();
  for (int i=0;i<nx*ny;i++)
    vec[i]*=2;
#else
#if 0
  for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
      A(i,j)=2*A(i,j);
#else
  A=2*A;
#endif
#endif
  printf("First element is now %12.4g.\n",A.elem(0,0));
  return octave_value(A);
  //octave_idx_type nelem (void)                                                                                                                                                                           
  return octave_value_list ();
}







