#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>


#include <healpix_cxx/healpix_base.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/arr.h>
#include <healpix_cxx/alm.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/pointing.h>

#ifdef __cplusplus
extern "C"
{
#endif
  //#include <chealpix.h>
#include <omp.h>
#define PI_OVER_TWO 1.5707963267948966


  //void map2alm_(int *nsmax, int *nlmax, int *nmmax, 
  //call map2alm*( nsmax, nlmax, nmmax, map_TQU, alm_TGC, zbounds, w8ring_TQU [, plm] )

  void map2alm_sc_d_(int *nsmax, int *nlmax, int *nmmax, double *map, double *alm, double *zbounds, double *w8ring);
  void __alm_tools_MOD_map2alm_sc_d(int *nsmax, int *nlmax, int *nmmax, double *map, double *alm, double *zbounds, double *w8ring);


#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;


/*--------------------------------------------------------------------------------*/
double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=(double)myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (map2alm_simple_c, args, nargout, "Transform a map into alms.\n")
{
  Matrix map=args(0).matrix_value();
  dim_vector dm=map.dims();
  int nelem=dm(0)*dm(1);
  int nside=sqrt(nelem/12);
  printf("nside is %d\n",nside);

  int lmax=nside;
  int mmax=lmax;
  dim_vector dm_out(2);
  dm(0)=lmax;
  dm(1)=2*mmax+1;
  ComplexMatrix alms(dm_out);
  double *mapptr=(double *)map.fortran_vec();
  double *almptr=(double *)alms.fortran_vec();
  __alm_tools_MOD_map2alm_sc_d(&nside,&lmax,&mmax,mapptr,almptr,NULL,NULL);
  return octave_value(alms);
 
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (test_healpix_arr, args, nargout, "Test healpix c++ array functionality.\n")
{
  Matrix map=args(0).matrix_value();
  dim_vector dm=map.dims();
  long nelem=dm(0)*dm(1);
  printf("have %d elements.\n",nelem);
  double *myptr=map.fortran_vec();
  arr<double> x(myptr,nelem);
  printf("allocated array.\n");
  double min,max;
  //x->minmax(&min,&max);
  printf("first element is %12.4f\n",x[nelem-1]);
  Healpix_Map <double> crap(x,RING);
  //Healpix_Map <double> crap();
  double val=2;
  //crap.fill(val);
  printf("average value is %12.4f\n",crap.average());
  Alm < xcomplex < double > > my_alms(12,12);
  map2alm(crap,my_alms,x,false);
  const arr <xcomplex < double > > vals=my_alms.Alms();
  //copyToPtr
  long nalm=vals.size();
  printf("have %d alms\n",nalm);
  ComplexColumnVector myans(nalm);
  complex<double> *ans_ptr=myans.fortran_vec();
  void *beg=(void *)vals.begin();
  memcpy(ans_ptr,beg,nalm*sizeof(complex<double>));
  //vals.copyToPtr(ans_ptr);
  return octave_value(myans);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (map2alm_intensity_c, args, nargout, "Transform healpix ring map into alms.\n")
{

  double t1=omp_get_wtime();

  Matrix map=args(0).matrix_value();
  int lmax=(int)get_value(args(1));
  int mmax=(int)get_value(args(2));
  Matrix weights=args(3).matrix_value();
  dim_vector dm=map.dims();
  long nelem=dm(0)*dm(1);
  dim_vector dm2=weights.dims();
  long nelem2=dm(0)*dm(1);
  if (nelem2!=nelem) {
    printf("Error in map2alm_intensity_c - mismatch in weights vs. map size %d vs. %d.\n",nelem2,nelem);
    return octave_value_list();    
  }
  
  
  double *myptr=map.fortran_vec();
  double *wtptr=weights.fortran_vec();
  //printf("did octave stuff at %12.6f\n",omp_get_wtime()-t1);
  arr<double> x(myptr,nelem);
  arr<double> wt(wtptr,nelem);

  Healpix_Map <double> crap(x,RING);
  Alm < xcomplex < double > > my_alms(lmax,mmax);
  //printf("did healpix allocs at %12.6f\n",omp_get_wtime()-t1);
  map2alm(crap,my_alms,wt,false);
  //printf("did alm transform at %12.6f\n",omp_get_wtime()-t1);
  //double t2=omp_get_wtime();
  //printf("map2alm took %12.4f seconds.\n",t2-t1);
  const arr <xcomplex < double > > vals=my_alms.Alms();
  //copyToPtr
  long nalm=vals.size();

  ComplexColumnVector myans(nalm);
  complex<double> *ans_ptr=myans.fortran_vec();
  void *beg=(void *)vals.begin();
  memcpy(ans_ptr,beg,nalm*sizeof(complex<double>));
  //printf("finishing up at %12.6f\n",omp_get_wtime()-t1);
  return octave_value(myans);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (pix2ang_c,args,nargout,"Interface to healpix c++ pix2ang.\n")
{
  int nside=(int)get_value(args(0));
  int32NDArray pix=args(1).array_value();
  dim_vector dm=pix.dims();
  
  int nelem=1;
  for (int i=0;i<dm.length();i++)
    nelem*=dm(i);
  
  printf("array has %d dimensions and %d members.\n",dm.length(),nelem);
  NDArray phi(dm);
  NDArray theta(dm);
  
  pointing tmp2(0.5,1.5);
  printf("pointing is %12.6f %12.6f\n",tmp2.theta,tmp2.phi);
  int *pixptr=(int *)pix.fortran_vec();
  double *phiptr=phi.fortran_vec();
  double *thetaptr=phi.fortran_vec();
  

  //T_Healpix_Base < int > base1 (nside,RING,nside);
  const int zero=0;
  //T_Healpix_Base < int > base1 (nside,RING,zero);
  T_Healpix_Base < int > base1 (nside,RING);

  for (int i=0;i<nelem;i++) {
    pointing tmp=base1.pix2ang(pixptr[i]);
    phiptr[i]=tmp.phi;
    thetaptr[i]=tmp.theta;
  }

  octave_value_list retval;
  retval(0)=theta;
  retval(1)=phi;
  return retval;

}
