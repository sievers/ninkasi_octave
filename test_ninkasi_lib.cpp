#include <octave/oct.h>
#include <iostream>
using namespace std;


#ifdef __cplusplus
extern "C"
{
#endif
#include <ninkasi.h>
  //#include <dirfile.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif




#include <string.h>

DEFUN_DLD (test_ninkasi_lib, args, nargout, "Testing my library.\n")
{
  octave_value_list retval;
  int nargin = args.length();
  if (nargin==0)
    return retval;
  
  charMatrix ch = args(0).char_matrix_value ();
  
  cout << "Hello.\n";
  cout << ch.row_as_string(0) << "\n";
  cout << "length is " << ch.length() << "\n";

  int nn=ch.length();

  char *s=ch.fortran_vec();

  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  printf("\n\nThis is still a test %s.\n",ss);



  mbTOD *tod=read_dirfile_tod_header(s );
  int n=tod->ndata;
  Matrix a(n,1);
  double *tmp=a.fortran_vec();
  for (int i=0;i<n;i++)
    tmp[i]=tod->az[i];
  
  
  
  printf("In octave, have %d elements.\n",n);
  
  
  

  return octave_value(a);


}       
