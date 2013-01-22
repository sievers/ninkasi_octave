#include <octave/oct.h>
#include <iostream>
#include <stdio.h>


#ifdef __cplusplus
extern "C"
{
#endif

#ifndef NO_FFTW
#include <fftw3.h>
#else
#include <fftw/fftw3.h>
#endif

#include <omp.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

/*--------------------------------------------------------------------------------*/

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=(double)myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (fft_omp_r2c_many, args, nargout, "Initialize FFTW with threads.\n")
{

  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();
  double *vec=mat.fortran_vec();
  int n=dm(0);
  int m=dm(1);
  int nn=n/2+1;
  
  ComplexMatrix matft(nn,m);
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();
  double t1=omp_get_wtime();
  int *nvec=(int *)malloc(sizeof(int)*m);
  for (int i=0;i<m;i++)
    nvec[i]=n;
  fftw_plan p=fftw_plan_many_dft_r2c(1,nvec,m,vec,NULL,1,n,cvec,NULL,1,nn,FFTW_ESTIMATE);
  fftw_execute(p);
  //printf("took %14.4e seconds.\n",omp_get_wtime()-t1);

  free(nvec);
  fftw_destroy_plan(p);
  return octave_value(matft);
				     
}
/*--------------------------------------------------------------------------------*/

octave_value fft_omp_r2c(Matrix mat)
{

  dim_vector dm=mat.dims();
  ComplexMatrix matft(dm);
  double *vec=mat.fortran_vec();
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();
  
  
  
  int n=dm(0);
  int m=dm(1);


  double t1=omp_get_wtime();


  fftw_plan p=fftw_plan_dft_r2c_1d(dm(0),vec,cvec,FFTW_MEASURE);  



#pragma omp parallel for shared(n,m,p,vec,cvec) default(none) 
  for (int i=0;i<m;i++) {
    fftw_execute_dft_r2c(p,&vec[i*n],&cvec[i*n]);
  }
  fftw_destroy_plan(p);
  double t2=omp_get_wtime();
  printf("Elapsed time is %12.4g seconds.\n",t2-t1); 

  return octave_value(matft);

}

/*--------------------------------------------------------------------------------*/

octave_value fft_omp_r2c_tight(Matrix mat)
{

  dim_vector dm=mat.dims();
  dim_vector cdm=mat.dims();
  int nn=(cdm(0)/2)+1;
  cdm(0)=nn;
  ComplexMatrix matft(cdm);
  double *vec=mat.fortran_vec();
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();

  
  
  int n=dm(0);
  int m=dm(1);
  //make sure you don't overwrite incoming data...
  double *vec2=(double *)malloc(sizeof(double)*n);

  fftw_plan p=fftw_plan_dft_r2c_1d(dm(0),vec2,cvec,FFTW_MEASURE);  
  free(vec2);
#pragma omp parallel for shared(n,nn,m,p,vec,cvec) default(none) 
  for (int i=0;i<m;i++) {
    fftw_execute_dft_r2c(p,&vec[i*n],&cvec[i*nn]);
  }
  fftw_destroy_plan(p);
  return octave_value(matft);

}
/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (fft_c2r_c, args, nargout, "Take fft along columns w/ openmp.\n") 
{

  int iseven=1;
  if (args.length()>1)
    iseven=(int)get_value(args(1));
  ComplexMatrix matft=args(0).complex_matrix_value();
  dim_vector cdm=matft.dims();
  dim_vector dm=matft.dims();
  int nn=cdm(0);
  int n;
  if (iseven)
    n=(nn-1)*2;
  else
    n=nn*2-1;
  dm(0)=n;
  int m =cdm(1);
  Matrix mat(dm);
  double *vec=mat.fortran_vec();

  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();
  int *nvec=(int *)malloc(sizeof(int)*m);
  for (int i=0;i<m;i++)
    nvec[i]=n;

  fftw_plan p=fftw_plan_many_dft_c2r(1,nvec,m,cvec,NULL,1,nn,vec,NULL,1,n,FFTW_ESTIMATE);
  double t1=omp_get_wtime();
  fftw_execute(p);
  //printf("execution time is %12.4e\n",omp_get_wtime()-t1);
  free(nvec);
  fftw_destroy_plan(p);

  double fac=1./n;
  for (int i=0;i<n*m;i++) {
    vec[i]*=fac;
  }
  return octave_value(mat);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (fft_r2c_many, args, nargout, "Take fft along columns w/ openmp.\n") 
{
  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();
  dim_vector dm2=dm;
  int nn=(dm(0)/2)+1;
  dm2(0)=nn;
  ComplexMatrix matft(dm2);

  double *vec=mat.fortran_vec();
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();
  

  int ntrans=dm(1);
  int len=dm(0);
  
  int *n=(int *)malloc(sizeof(int)*ntrans);
  for (int i=0;i<ntrans;i++)
    n[i]=len;
  fftw_plan plan=fftw_plan_many_dft_r2c(1,n,ntrans,vec,NULL,1,len,cvec,NULL,1,nn,FFTW_MEASURE);
  double t1=omp_get_wtime();
  fftw_execute(plan);
  double t2=omp_get_wtime();
  printf("Took %8.5f seconds to execute plan.\n",t2-t1);
  fftw_destroy_plan(plan);
  free(n);
  return octave_value(matft);

}
/*--------------------------------------------------------------------------------*/

octave_value fft_omp_c2r_tight(ComplexMatrix matft,int iseven)
{

  dim_vector cdm=matft.dims();
  dim_vector dm=matft.dims();
  int nn=cdm(0);
  int n;
  if (iseven)
    n=(nn-1)*2;
  else
    n=nn*2-1;
  dm(0)=n;
  int m =cdm(1);
  Matrix mat(dm);
  double *vec=mat.fortran_vec();
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();

  fftw_complex *cvec_tmp=(fftw_complex *)malloc(sizeof(fftw_complex)*nn);  
  fftw_plan p=fftw_plan_dft_c2r_1d(dm(0),cvec_tmp,vec,FFTW_MEASURE);  
  free(cvec_tmp);

#pragma omp parallel for shared(n,nn,m,p,vec,cvec) default(none) 
  for (int i=0;i<m;i++) {
    fftw_execute_dft_c2r(p,&cvec[i*nn],&vec[i*n]);
  }
  fftw_destroy_plan(p);
  return octave_value(mat);

}

/*--------------------------------------------------------------------------------*/

octave_value fft_omp_c2c(ComplexMatrix mat)
{

  dim_vector dm=mat.dims();
  ComplexMatrix matft(dm);
  fftw_complex *vec=(fftw_complex *)mat.fortran_vec();
  fftw_complex *cvec=(fftw_complex *)matft.fortran_vec();


  int n=dm(0);
  int m=dm(1);
  
  fftw_plan p=fftw_plan_dft_1d(n,vec,cvec,FFTW_FORWARD,FFTW_MEASURE);
#pragma omp parallel for shared(n,m,p,vec,cvec) default(none) 
  for (int i=0;i<m;i++) {
    fftw_execute_dft(p,&vec[i*n],&cvec[i*n]);
  }
  fftw_destroy_plan(p);
  return octave_value(matft);

}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (fftw_init_threaded, args, nargout, "Initialize FFTW with threads.\n")
{
  fftw_init_threads();
  int nthread;
  if (args.length()>0)
    nthread=(int)get_value(args(0));
  else {
#pragma omp parallel
#pragma omp single
    nthread=omp_get_num_threads();
  }
  
  printf("planning with %d thread.\n",nthread);
  fftw_plan_with_nthreads(nthread);
  
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (fft_r2c, args, nargout, "Take fft along columns w/ openmp.\n") 
{
  //fprintf(stderr,"Taking tight fft.\n");
  return fft_omp_r2c_tight(args(0).matrix_value());
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (fft_c2r_c_old, args, nargout, "Take fft along columns w/ openmp.\n") 
{
  int iseven=1;
  if (args.length()>1)
    iseven=(int)get_value(args(1));
  
  //fprintf(stderr,"Taking tight fft.\n");
  return fft_omp_c2r_tight(args(0).complex_matrix_value(),iseven);
  //return fft_omp_c2r_many(args(0).complex_matrix_value(),iseven);

}
/*--------------------------------------------------------------------------------*/

//DEFUN_DLD (fft_c2r, args, nargout, "Take fft along columns w/ openmp.\n") 
//{
// return fft_omp_r2c(args(0).matrix_value());
//}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (fft_omp, args, nargout, "Take fft along columns w/ openmp.\n")

{

  if (args(0).is_complex_matrix()) 
    return fft_omp_c2c(args(0).complex_matrix_value());
  else
    return fft_omp_r2c(args(0).matrix_value());
  
  return octave_value_list();  //never get here

		       
}
