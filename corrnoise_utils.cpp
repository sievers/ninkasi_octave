#include <iostream>
#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#ifdef __cplusplus
extern "C"
{
#endif
  void   dpotrf_(char *uplo, int *n, double *a, int *lda, int *info, int uplolen);
  void   dpotrs_(char *uplo, int *n, int *nrhs,double *a, int *lda, double *b, int *ldb, int *info, int uplolen);

  
#ifdef __cplusplus
}  /* end extern "C" */
#endif


void invert_3x3_matrix(double *mat, double *mat2)
{

  //d=x(1,1)*x(2,2)*x(3,3)-x(1,1)*x(2,3)^2-x(1,2)^2*x(3,3)+2*x(1,2)*x(2,3)*x(1,3)-x(1,3)^2*x(2,2);  
  //double det=x[0][0]*x[1][1]*x[2][2]-x[0][0]*x[1][2]*x[1][2]-x[0][1]*x[0][1]*x[2][2]+2*x[0][1]*x[1][2]*x[0][2]-x[0][2]*x[0][2]*x[1][1];
  
  double det=mat[0*3+0]*mat[1*3+1]*mat[2*3+2]-mat[0*3+0]*mat[1*3+2]*mat[1*3+2]-mat[0*3+1]*mat[0*3+1]*mat[2*3+2]+2*mat[0*3+1]*mat[1*3+2]*mat[0*3+2]-mat[0*3+2]*mat[0*3+2]*mat[1*3+1];
  det=1.0/det;

  mat2[0*3+0]=det*(mat[1*3+1]*mat[2*3+2]-mat[1*3+2]*mat[1*3+2]);
  mat2[1*3+1]=det*(mat[0*3+0]*mat[2*3+2]-mat[0*3+2]*mat[0*3+2]);
  mat2[2*3+2]=det*(mat[1*3+1]*mat[0*3+0]-mat[0*3+1]*mat[0*3+1]);
  mat2[0*3+1]=det*(mat[1*3+2]*mat[0*3+2]-mat[0*3+1]*mat[2*3+2]);
  mat2[0*3+2]=det*(mat[0*3+1]*mat[1*3+2]-mat[1*3+1]*mat[0*3+2]);
  mat2[1*3+2]=det*(mat[0*3+1]*mat[0*3+2]-mat[1*3+2]*mat[0*3+0]);

  mat2[2*3+1]=mat2[1*3+2];
  mat2[1*3+0]=mat2[0*3+1];
  mat2[2*3+0]=mat2[0*3+2];
  
}


/*--------------------------------------------------------------------------------*/

void solve_posdef_matrix(double *a, double *b, int n, int nrhs)
{
  int info;
  char uplo='u';
  //printf("nrhs is %d\n",nrhs);
  dpotrf_(&uplo,&n,a,&n,&info,1);
  if (info) {
    printf("Non-positive-definite in solve_posdef_matrix.\n");
    return;
  }

  dpotrs_(&uplo,&n,&nrhs,a,&n,b,&n,&info,1);


}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (test_solve_posdef, args, nargout, "Multiply correlated noise w/ priors.  Args are map,prior,matrix.\n")
{

  Matrix a=args(0).matrix_value();
  Matrix b=args(1).matrix_value();
  dim_vector dm=b.dims();
  solve_posdef_matrix(a.fortran_vec(),b.fortran_vec(),dm(0),dm(1));
  return octave_value(b);
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (multiply_corrnoise_wprior_c, args, nargout, "Multiply correlated noise w/ priors.  Args are map,prior,matrix.\n")
{
  Matrix map=args(0).matrix_value();
  Matrix prior=args(1).matrix_value();
  Matrix mat=args(2).matrix_value();
  
  dim_vector dims=map.dims();
  dim_vector mat_dims=mat.dims();
  int npar=mat_dims(0);
  if (dims(0)!=mat_dims(0)) {
    fprintf(stderr,"Error in matrix vs. map size.\n");
    return octave_value_list();
  }
  dim_vector dims2=prior.dims();
  if ((dims2(0)!=dims(0))||(dims2(1)!=dims(1)))  {
    fprintf(stderr,"Error in map vs. prior size.\n");
    return octave_value_list();
  }
  
  
  
  double *vec=mat.fortran_vec();
  double *matinv=(double *)malloc(sizeof(double)*npar*npar);
  
  //  if (npar!=3)  {
  //  fprintf(stderr,"Only 3x3 supported now.\n");
  //  return octave_value_list();
  // }
    

  double *mapptr=map.fortran_vec();
  
  Matrix prod(dims);
  double *prod_ptr=prod.fortran_vec();

  double *prior_ptr=prior.fortran_vec();

  double *matcopy=(double *)malloc(sizeof(double)*npar*npar);
  int ndata=dims(1);
  for (int i=0;i<ndata;i++) {
    memcpy(matcopy,vec,sizeof(double)*npar*npar);
    
#if 1
    for (int j=0;j<npar;j++)
      matcopy[j*npar+j]+=prior_ptr[npar*i+j];
    memcpy(prod_ptr,mapptr,sizeof(double)*npar);
    solve_posdef_matrix(matcopy,prod_ptr,npar,1);
#else
    matcopy[0*npar+0]+=prior_ptr[npar*i];
    matcopy[1*npar+1]+=prior_ptr[npar*i+1];
    matcopy[2*npar+2]+=prior_ptr[npar*i+2];
    invert_3x3_matrix(matcopy,matinv);
    
    for (int j=0;j<npar;j++) {
      prod_ptr[i*npar+j]=0;
      for (int k=0;k<npar;k++)
	prod_ptr[i*npar+j]+=mapptr[i*npar+k]*matinv[j+k*npar];
    }
#endif
    
  }
  
  

  
  free(matinv);
  free(matcopy);  
  return octave_value(prod);
}

