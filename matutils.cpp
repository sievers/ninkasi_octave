#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <time.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include <omp.h>
  void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc, int uplolen, int translen);
  void dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info, int uplolen);
#include <nk_clapack.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif



//void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc, int uplolen, int translen);

//void dgemm_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc, int uplolen, int translen);


/*-----------------------------------------------------------------------------------------*/

DEFUN_DLD (dpotrs, args, nargout, "Multiply a matrix by its transposed self..\n")
{
  Matrix a=args(0).matrix_value();
  dim_vector dm=a.dims();

  Matrix vecs=args(1).matrix_value();
  dim_vector dm2=vecs.dims();
  int n=dm(0);
  int nrhs=dm2(1);
  char uplo='U';
  int info;

  dpotrs_(&uplo,&n,&nrhs,a.fortran_vec(),&n,vecs.fortran_vec(),&n,&info,1);
  return octave_value(vecs);

}
/*-----------------------------------------------------------------------------------------*/

DEFUN_DLD (dsyrk, args, nargout, "Multiply a matrix by its transposed self..\n")
{
  Matrix a=args(0).matrix_value();
  dim_vector dm=a.dims();

  dim_vector dm2=a.dims();
  dm2(0)=dm2(1);
  Matrix c(dm2);

  char uplo='u';
  char trans='t';
  int k=dm(0);
  int n=dm(1);

  double alpha=1.0;
  double beta=0.0;
  int lda=k;
  int ldc=n;

  int info=0;
  //clapack_dpotrf(uplo,n,a.fortran_vec(),n,&info);
  
  
  
  //clapack_dsyrk(uplo,trans,n,k,alpha,a.fortran_vec(),lda,beta,c.fortran_vec(),ldc);
  
  int nt;
#pragma omp parallel shared(nt)
#pragma omp single
  nt=omp_get_num_threads();




  memset(c.fortran_vec(),0,n*n*sizeof(double));
  
  if (k>nt)
    {
#pragma omp parallel shared(k,n,dm2,alpha,beta,lda,ldc,c,a,uplo,trans,nt) default(none)
      {

	double *myptr=(double *)malloc(sizeof(double)*n*n);
	memset(myptr,0,n*n*sizeof(double));
	int myid=omp_get_thread_num();
	int myk0=((k+nt-1)/nt)*myid;
	int myk1=((k+nt-1)/nt)*(myid+1)-1;
	if (myk1>=k)
	  myk1=k-1;
	int kk=1+(myk1-myk0);
	double *dptr=a.fortran_vec();
	dsyrk_(&uplo,&trans,&n,&kk,&alpha,dptr+myk0,&lda,&beta,myptr,&ldc,1,1);
	
#pragma omp critical 
	{
	  double *cptr=c.fortran_vec();
	  for (int i=0;i<n*n;i++)
	    cptr[i]+=myptr[i];
	}
	free(myptr);
	
      }
      
#if 1
      double *cptr=c.fortran_vec();
      for (int i=0;i<n;i++)
	for (int j=i;j<n;j++)
	  cptr[i*n+j]=cptr[j*n+i];
#endif
    }
  else
    dsyrk_(&uplo,&trans,&n,&k,&alpha,a.fortran_vec(),&lda,&beta,c.fortran_vec(),&ldc,1,1);
  //dsyrk_(&uplo,&trans,&n,&k,&alpha,a.fortran_vec(),&lda,&beta,c.fortran_vec(),&ldc,1,1);
  return octave_value(c);
  
}



/*--------------------------------------------------------------------------------*/



DEFUN_DLD (dsyrk2, args, nargout, "Multiply a matrix by its transposed self..\n")
{
  Matrix a=args(0).matrix_value();
  dim_vector dm=a.dims();

  dim_vector dm2=a.dims();
  dm2(1)=dm2(0);
  Matrix c(dm2);

  char uplo='u';
  char trans='n';
  int k=dm(1);
  int n=dm(0);

  double alpha=1.0;
  double beta=0.0;
  int lda=n;
  int ldc=n;

  int info=0;
  //clapack_dpotrf(uplo,n,a.fortran_vec(),n,&info);
  
  
  
  //clapack_dsyrk(uplo,trans,n,k,alpha,a.fortran_vec(),lda,beta,c.fortran_vec(),ldc);
  
  int nt;
#pragma omp parallel shared(nt)
#pragma omp single
  nt=omp_get_num_threads();


  printf("k and n are %d %d\n",k,n);

  memset(c.fortran_vec(),0,n*n*sizeof(double));
  
  if (k>nt)
    {
#pragma omp parallel shared(k,n,dm2,alpha,beta,lda,ldc,c,a,uplo,trans,nt) default(none)
      {

	double *myptr=(double *)malloc(sizeof(double)*n*n);
	memset(myptr,0,n*n*sizeof(double));
	int myid=omp_get_thread_num();
	int myk0=((k+nt-1)/nt)*myid;
	int myk1=((k+nt-1)/nt)*(myid+1)-1;
	if (myk1>=k)
	  myk1=k-1;
	int kk=1+(myk1-myk0);
	double *dptr=a.fortran_vec();
	dsyrk_(&uplo,&trans,&n,&kk,&alpha,dptr+myk0*n,&lda,&beta,myptr,&ldc,1,1);
#pragma omp critical 
	{
	  double *cptr=c.fortran_vec();
	  for (int i=0;i<n*n;i++)
	    cptr[i]+=myptr[i];
	}
	free(myptr);
	
      }
      
#if 1
      double *cptr=c.fortran_vec();
      for (int i=0;i<n;i++)
	for (int j=i;j<n;j++)
	  cptr[i*n+j]=cptr[j*n+i];
#endif
    }
  else
    dsyrk_(&uplo,&trans,&n,&k,&alpha,a.fortran_vec(),&lda,&beta,c.fortran_vec(),&ldc,1,1);
  //dsyrk_(&uplo,&trans,&n,&k,&alpha,a.fortran_vec(),&lda,&beta,c.fortran_vec(),&ldc,1,1);
  return octave_value(c);
  
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD(invsafe_c,args,nargout,"Invert a possibly singular matrix.  args are mat,[thresh].  any eigenvalues less than thresh*max are set to zero.\n")
{
  Matrix mat=args(0).matrix_value();
  double thresh=1e-8;
  if (args.length()>1) {
    Matrix tmp=args(1).matrix_value();
    thresh=tmp(0,0);
    //printf("setting threshold to %14.5g\n",thresh);
  }
  double *mm=mat.fortran_vec();
  dim_vector dm=mat.dims();
  dinvsafe(mm,dm(1),thresh);
  return octave_value(mat);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (scale_rowcol_mat, args, nargout, "Test complex stuff.\n")
{
  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();

  Matrix to_scale=args(1).matrix_value();
  dim_vector dm_scale=to_scale.dims();
  
  double *matptr=mat.fortran_vec();
  double *scaleptr=to_scale.fortran_vec();

  if ((dm_scale(0)==1)&&(dm_scale(1)==dm(1))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[j];
    return octave_value(mat);
  }


  if ((dm_scale(1)==1)&&(dm_scale(0)==dm(0))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[i];
    return octave_value(mat);
  }

  //If we got here, dimensions didn't line up.  Now try to guess.
  //guaranteed non-square if we get a match here, as we would have
  //matched previously.

  if ((dm_scale(1)==1)&&(dm_scale(0)==dm(1))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[j];
    return octave_value(mat);
  }


  if ((dm_scale(0)==1)&&(dm_scale(1)==dm(0))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[i];
    return octave_value(mat);
  }


  fprintf(stderr,"Warning - didn't find an appropriate size match in scale_rowcol_mat - %d %d vs. %d %d.\n",dm(0),dm(1),dm_scale(0),dm_scale(1));
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
#if 0
DEFUN_DLD (scale_rowcol_complex, args, nargout, "Test complex stuff.\n")
{
  ComplexMatrix mat=args(0).complex_matrix_value();
  dim_vector dm=mat.dims();

  Matrix to_scale=args(1).matrix_value();
  dim_vector dm_scale=to_scale.dims();
  
  double *matptr=mat.fortran_vec();
  double *scaleptr=to_scale.fortran_vec();

  if ((dm_scale(0)==1)&&(dm_scale(1)==dm(1))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[j];
    return octave_value(mat);
  }


  if ((dm_scale(1)==1)&&(dm_scale(0)==dm(0))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[i];
    return octave_value(mat);
  }

  //If we got here, dimensions didn't line up.  Now try to guess.
  //guaranteed non-square if we get a match here, as we would have
  //matched previously.

  if ((dm_scale(1)==1)&&(dm_scale(0)==dm(1))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[j];
    return octave_value(mat);
  }


  if ((dm_scale(0)==1)&&(dm_scale(1)==dm(0))) {
    int n=dm(0);
    int m=dm(1);
#pragma omp parallel for shared(n,m,matptr,scaleptr) default(none)
    for (int j=0;j<m;j++)
      for (int i=0;i<n;i++)
	matptr[j*n+i]*=scaleptr[i];
    return octave_value(mat);
  }


  fprintf(stderr,"Warning - didn't find an appropriate size match in scale_rowcol_mat.\n");
  return octave_value_list();

}
#endif
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (test_complex, args, nargout, "Test complex stuff.\n")
{
  ComplexMatrix mat=args(0).complex_matrix_value();

  if (args(0).is_complex_matrix()) {
    printf("doing complex.\n");
  }
  else
    return octave_value_list();


  double *vec=(double *)mat.fortran_vec();
  dim_vector dm=mat.dims();
  int numel=dm(0)*dm(1);
  for (int i=0;i<numel*2;i++) {
    printf("%5d %12.5f\n",i,vec[i]);
  }
  
  

  return octave_value(mat);
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (colmult, args, nargout, "Multiply columns of a matrix by a row vector.\n")
{
  if (args.length()<2) {
    printf("Need at least two arguments to colmult.\n");
    return octave_value_list();
  }
  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();
  Matrix vec=args(1).matrix_value();
  dim_vector dv=vec.dims();

  int nrow=dm(0);
  int ncol=dm(1);
  int nelem=dv(0)*dv(1);
  if (nelem!=nrow) {
    printf("size mismatch in colmult: %d %d\n",nrow,nelem);
    return octave_value_list();
  }
  int i,j;
  double *matptr=mat.fortran_vec();
  double *vecptr=vec.fortran_vec();
  //#pragma omp parallel for shared(matptr,vecptr,nrow,ncol) private(j) default(none)
  for (i=0;i<ncol;i++)
    for (j=0;j<nrow;j++)
      matptr[i*nrow+j]*=vecptr[i];

  return octave_value(mat);
  

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(unix_ctime,args,nargout,"Return the curent ctime.\n")
{
  time_t tt=time(NULL);
  return octave_value(tt);
}
