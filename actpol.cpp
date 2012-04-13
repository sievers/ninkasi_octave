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

DEFUN_DLD ( ACTpolArray_init, args, nargout, "Initialize an ACTpol pointing array.  Args are (x,y,angle,[freq])\n")
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
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(initialize_actpol_pointing,args,nargout,"Initialize a TOD with actpol pointing.  Args are tod,dx,dy,(angle),freq,dpiv\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));

  ColumnVector x=args(1).column_vector_value();
  actData *dx=x.fortran_vec();
  ColumnVector y=args(2).column_vector_value();
  actData *dy=y.fortran_vec();
  ColumnVector a=args(3).column_vector_value();
  actData *angle;
  if (a.length()==tod->ndet)
    angle=a.fortran_vec();
  else
    angle=NULL;
  actData freq=get_value(args(4));
  int dpiv=(int)get_value(args(5));
  initialize_actpol_pointing(tod,dx,dy,angle,freq,dpiv);
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(precalc_actpol_pointing,args,nargout,"Precalculated stuff for an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  precalc_actpol_pointing(tod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(precalc_actpol_pointing_free,args,nargout,"Free precalculated stuff from an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  precalc_actpol_pointing_free(tod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD(get_detector_radec_actpol_c, args, nargout, "Get RA/Dec of a detector.\n")
{
  int nargin = args.length();
  if (nargin<2) {
    printf("Only have %d arguments in get_detector_radec_c\n",nargin);
    return octave_value_list();
  }

  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int det=(int)get_value(args(1));
  PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
  
  get_radec_one_det_actpol(tod,det,scratch);
  
  
  Matrix radec(tod->ndata,2);
  double *dat=radec.fortran_vec();
  for (int i=0;i<tod->ndata;i++) {
    dat[i]=scratch->ra[i];
    dat[i+tod->ndata]=scratch->dec[i];
  }
  destroy_pointing_fit_scratch(scratch);

  return octave_value(radec);

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_radec_from_altaz_actpol_c,args,nargout,"Convert az, el, and ctime to ra/dec and sin/cos(2gamma).\n")
{
  int nargin = args.length();
  if (nargin<3) {
    printf("Only have %d arguments in get_radec_from_altaz_actpol_c, need at least 3.\n",nargin);
    return octave_value_list();
  }
  Matrix azmat=args(0).matrix_value();
  Matrix elmat=args(1).matrix_value();
  Matrix tmat=args(2).matrix_value();
  double dx=0;
  if (nargin>3)
    dx=get_value(args(3));
  double dy=0;
  if (nargin>4)
    dy=get_value(args(4));

  double theta=0;
  if (nargin>5)
    theta=get_value(args(5));


  double *az=azmat.fortran_vec();
  double *el=elmat.fortran_vec();
  double *tvec=tmat.fortran_vec();
  dim_vector dm=azmat.dims();
  int nelem=dm(0)*dm(1);
  double azmin=az[0];
  double azmax=az[0];
  
  Matrix ramat(dm);
  Matrix decmat(dm);
  Matrix sin2gammamat(dm);
  Matrix cos2gammamat(dm);
  double *ra=ramat.fortran_vec();
  double *dec=decmat.fortran_vec();
  double *sin2gamma=sin2gammamat.fortran_vec();
  double *cos2gamma=cos2gammamat.fortran_vec();
  

  for (int i=0;i<nelem;i++) {
    if (az[i]<azmin)
      azmin=az[i];
    if (az[i]>azmax)
      azmax=az[i];
  }
  double az_cent=(azmin+azmax)/2;
  double az_throw=(azmin-azmax)/2;


  ACTpolPointingFit *fit=(ACTpolPointingFit *)calloc(1,sizeof(ACTpolPointingFit));
  assert(fit);
  int nhorns=1;
  ACTpolArray *array = ACTpolArray_alloc(nhorns);
  assert(array);

  double xcent=0.0;
  double ycent=0.0;
  double freq=148.0;
  ACTpolArray_init(array, freq, xcent,ycent);
  
  ACTpolFeedhorn_init(array->horn,dx,dy,theta);
  ACTpolWeather weather;
  ACTpolWeather_default(&weather);


  ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
  ACTpolArrayCoords_init(coords);

  ACTpolState *state = ACTpolState_alloc();
  ACTpolState_init(state);

  
  ACTpolScan scan;
  ACTpolScan_init(&scan, el[0], az_cent,az_throw);
  ACTpolArrayCoords_update_refraction(coords, &scan, &weather);
  for (int i = 0; i < nelem; i++)
    {
      ACTpolState_update(state, tvec[i],el[i],az[i]);
      ACTpolArrayCoords_update(coords, state);
      ACTpolFeedhornCoords *fc = coords->horn;
      ra[i]=fc->ra;
      dec[i]=fc->dec;
      sin2gamma[i]=fc->sin2gamma;
      cos2gamma[i]=fc->cos2gamma;
    }
  octave_value_list retval;
  retval(0)=ramat;
  retval(1)=decmat;
  retval(2)=sin2gammamat;
  retval(3)=cos2gammamat;
  return retval;

}
/*--------------------------------------------------------------------------------*/
