#include <octave/oct.h>
#include <octave/ov-struct.h>

#include <iostream>
using namespace std;

#include <mpi.h>
#ifdef __cplusplus
extern "C"
{
#endif
  
#include <ninkasi_config.h>
#include <ninkasi.h>
#include <dirfile.h>
#include <slalib.h>
#include <actpol/actpol.h>
#include <assert.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif




/*--------------------------------------------------------------------------------*/

actData get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  actData myval=(actData)myptr(0,0);
  return myval;

}



/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

}

 /*--------------------------------------------------------------------------------*/

DEFUN_DLD ( ACTpolArray_init, args, nargout, "Initialize and ACTpol pointing array.  Args are (x,y,angle,[freq])\n")
{
  
  ColumnVector x=args(0).column_vector_value();
  ColumnVector y=args(1).column_vector_value();
  ColumnVector a=args(2).column_vector_value();
  double *xx=x.fortran_vec();
  double *yy=y.fortran_vec();
  double *aa=a.fortran_vec();



  double freq_ghz=148;
  if (args.length()>=4)
    freq_ghz=get_value(args(3));

  int nhorns=x.length();
  assert(y.length()==nhorns);
  assert(a.length()==nhorns);
  printf("have %d detectors.\n",nhorns);
  
  double xtot=0,ytot=0;
  for (int i=0;i<nhorns;i++){
    xtot+=xx[i];
    ytot+=yy[i];
  }
  double xcent=xtot/((double)nhorns);
  double ycent=ytot/((double)nhorns);
  

  
  ACTpolArray *array = ACTpolArray_alloc(nhorns);
  printf("Hi, made it here.  Centers are %14.6f %14.6f\n",xcent,ycent);
  ACTpolArray_init(array, freq_ghz, xcent,ycent);
  
  for (int i = 0; i < nhorns; i++) {
    ACTpolFeedhorn_init(array->horn+i, xx[i], yy[i], aa[i]);
  }


  int64NDArray myptr(1);
  myptr(0)=(long)array;

  return octave_value(myptr);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD ( ACTpolPointingEvaluate, args, nargout, "Given an ACTpol structure, get the pointing.  Args are (tod,detarray,[weather, currently unsupported]\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  ACTpolArray *array=(ACTpolArray *)get_pointer(args(1));
  if (tod->ndet!=array->nhorns) {
    printf("Error in ACTpolPointingEvaluate.  Mismatch in horns, tod has %d, array has %d\n",tod->ndet,array->nhorns);
    return octave_value_list();
  }
  assert(tod->ndet==array->nhorns);  //make sure we agree on how many detectors, should keep detectors on same horn separate
  
  Matrix ra_mat(tod->ndata,tod->ndet);
  //Matrix ra_mat(tod->ndet,tod->ndata);
  Matrix dec_mat(tod->ndata,tod->ndet);
  //Matrix dec_mat(tod->ndet,tod->ndata);
  double *ra=ra_mat.fortran_vec();
  double *dec=dec_mat.fortran_vec();
  ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
  ACTpolArrayCoords_init(coords);

  double azmin,azmax,altmin,altmax;
  azmin=tod->az[0];
  azmax=tod->az[0];
  altmin=tod->alt[0];
  altmax=tod->alt[0];
  for (int j=1;j<tod->ndata;j++) {
    if (tod->az[j]<azmin)
      azmin=tod->az[j];
    if (tod->az[j]>azmax)
      azmax=tod->az[j];
    if (tod->alt[j]<altmin)
      altmin=tod->alt[j];
    if (tod->alt[j]>altmax)
      altmax=tod->alt[j];
  }
  printf("alt/az limits are %14.6f %14.6f %14.6f %14.6f\n",altmin,altmax,azmin,azmax);

  
  ACTpolState *state = ACTpolState_alloc();
  ACTpolState_init(state);
  

  ACTpolScan scan;
  ACTpolScan_init(&scan, 0.5*(altmin+altmax), 0.5*(azmin+azmax), 0.5*(azmax-azmin));

  
  ACTpolWeather weather;
  ACTpolWeather_default(&weather);
  
  ACTpolArrayCoords_update_refraction(coords, &scan, &weather);
  
  


  
  double *raptr=ra_mat.fortran_vec();
  double *decptr=dec_mat.fortran_vec();

  int ndet=tod->ndet;
  for (int isamp=0;isamp<tod->ndata;isamp++) {
    double ctime=tod->ctime+((double)isamp)*(tod->deltat);
    ACTpolState_update(state, ctime, tod->alt[isamp],tod->az[isamp]);
    ACTpolArrayCoords_update_fast(coords, state);
    double *ra=raptr+isamp*(long)ndet;
    double *dec=decptr+isamp*(long)ndet;
    for (int k=0;k<ndet;k++) {
      ACTpolFeedhornCoords *fc = coords->horn + k;
      ra[k]=fc->ra;
      dec[k]=fc->sindec;
      //ra[k]=coords->horn[k].ra;
      //dec[k]=coords->horn[k].sindec;
      
    }
  }
  
  octave_value_list retval;
  retval(0)=ra_mat;
  retval(1)=dec_mat;
  return retval;
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (FindTodPixLimits, args, nargout, "Given a tod and a map, find the range of indices spanned by the tod in the map.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  MAP *map=(MAP *)get_pointer(args(1));
  int imin,imax;
  find_map_index_limits(map,tod,&imin,&imax);
  octave_value_list retval;
  retval(0)=imin;
  retval(1)=imax;
  return retval;
  
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (tod2map_actpol, args, nargout, "Project a tod into a map with actpoley goodness.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  MAP *map=(MAP *)get_pointer(args(1));
  Matrix ipiv=args(2).matrix_value();
  dim_vector dm=ipiv.dims();
  int nelem=dm(0)*dm(1);
  double *vv=ipiv.fortran_vec();
  int *ii=(int *)malloc(sizeof(int)*nelem);
  for (int i=0;i<nelem;i++)
    ii[i]=vv[i];
  tod2map_actpol(map, tod, ii);
  free(ii);

  return octave_value_list();
}
