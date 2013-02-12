#include <octave/oct.h>
#include <octave/ov-struct.h>

#include <iostream>
using namespace std;

#include <mpi.h>
#ifdef __cplusplus
extern "C"
{
#endif

#define ACTPOL  
#include <ninkasi_config.h>
#include <ninkasi.h>
#include <noise.h>
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
DEFUN_DLD(get_tod_actpol_pointing_offsets_c,args,nargout,"Return detector pointing offsets.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  ACTpolPointingFit *fit=tod->actpol_pointing;
  Matrix offsets(tod->ndet,2);
  double *offs=offsets.fortran_vec();
  for (int i=0;i<tod->ndet;i++) {
    //offs[2*i]=fit->dx[i];
    //offs[2*i+1]=fit->dy[i];
    offs[i]=fit->dx[i];
    offs[tod->ndet+i]=fit->dy[i];

  }
  return octave_value(offsets);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(get_detector_offsets_actpol,args,nargout,"Pull the detector offsets stored in a TOD with actpol pointing.  Arg is tod.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (tod->actpol_pointing==NULL) {
    fprintf(stderr,"don't have actpol pointing in this tod, cannot pull detector offsets.\n");
    return octave_value_list();
  }
  int ndet=tod->ndet;
  Matrix dx(ndet,1);
  Matrix dy(ndet,1);
  Matrix theta(ndet,1);

  double *dxptr=dx.fortran_vec();
  double *dyptr=dy.fortran_vec();
  double *thptr=theta.fortran_vec();
  ACTpolPointingFit *fit=tod->actpol_pointing;
  
  for (int i=0;i<ndet;i++) {
    dxptr[i]=fit->dx[i];
    dyptr[i]=fit->dy[i];
    thptr[i]=fit->theta[i];
  }
  octave_value_list retval;
  retval(0)=dx;
  retval(1)=dy;
  retval(2)=theta;
  return retval;
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
  if (tod->actpol_pointing)
    update_actpol_pointing(tod,dx,dy,angle,freq,dpiv);
  else
    initialize_actpol_pointing(tod,dx,dy,angle,freq,dpiv);
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(precalc_actpol_pointing_exact,args,nargout,"Precalculated stuff for an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  int op_flag=NINKASI_DO_RADEC|NINKASI_DO_TWOGAMMA;
  if (args.length()>1)
    op_flag=(int)get_value(args(1));
  precalc_actpol_pointing_exact(tod,op_flag);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(find_tod_radec_lims_actpol_pointing_exact_c,args,nargout,"Find the TOD ra/dec limits without storing the full pointing.  Arguments are tod,(wrap).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  actData rawrap=3.14159265;
  if (args.length()>1)
    rawrap=get_value(args(1));
  find_tod_radec_lims_actpol_pointing_exact(tod,rawrap);
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

  double zero=0;
  double *dx=&zero,*dy=&zero;
  Matrix dxmat,dymat;
  int ndet=1;
  
  if (nargin>3) {
    dxmat=args(3).matrix_value();
    dx=dxmat.fortran_vec();
    dim_vector dm=dxmat.dims();
    ndet=dm(0)*dm(1);
  }
  if (nargin>4) {
    dymat=args(4).matrix_value();
    dy=dymat.fortran_vec();
    dim_vector dm=dymat.dims();
    assert(dm(0)*dm(1)==ndet);
  }

  double *theta;
  Matrix theta_mat;
  int free_theta=0;
  if (nargin>5) {
    theta_mat=args(5).matrix_value();
    dim_vector dm=theta_mat.dims();
    if (dm(0)*dm(1)==1) {
      theta=(double *)malloc(sizeof(double)*ndet);
      free_theta=1;
      for (int i=0;i<ndet;i++)
	theta[i]=theta_mat(1,1);
    }
    else {
      assert(dm(0)*dm(1)==ndet);
      theta=theta_mat.fortran_vec();
    }
  }
  else {
    theta=(double *)malloc(sizeof(double)*ndet);
    free_theta=1;
    for (int i=0;i<ndet;i++)
      theta[i]=0;
  }
  
  
  double *az=azmat.fortran_vec();
  double *el=elmat.fortran_vec();
  double *tvec=tmat.fortran_vec();
  dim_vector dm=azmat.dims();
  int nelem=dm(0)*dm(1);
  double azmin=az[0];
  double azmax=az[0];
  
  Matrix ramat(nelem,ndet);
  Matrix decmat(nelem,ndet);
  Matrix sin2gammamat(nelem,ndet);
  Matrix cos2gammamat(nelem,ndet);
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
  //Fix this - JLS 23 June 2012.  Note to Nolta - for constant az, this makes no difference.
  double az_cent=(azmin+azmax)/2;
  double az_throw=(azmin-azmax)/2;


  ACTpolPointingFit *fit=(ACTpolPointingFit *)calloc(1,sizeof(ACTpolPointingFit));
  assert(fit);
  int nhorns=ndet;
  ACTpolArray *array = ACTpolArray_alloc(nhorns);
  assert(array);

  double xcent=0.0;
  double ycent=0.0;
  double freq=148.0;
  ACTpolArray_init(array, freq, xcent,ycent);
  for (int i=0;i<nhorns;i++) {
    ACTpolFeedhorn_init(&(array->horn[i]),dx[i],dy[i],theta[i]);
  }

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
      for (int j=0;j<nhorns;j++) {
	ACTpolFeedhornCoords *fc = &(coords->horn[j]);
	ra[i+j*nelem]=fc->ra;
	dec[i+j*nelem]=fc->dec;
	sin2gamma[i+j*nelem]=fc->sin2gamma;
	cos2gamma[i+j*nelem]=fc->cos2gamma;
      }
    }
  if (free_theta==1)
    free(theta);
  octave_value_list retval;
  retval(0)=ramat;
  retval(1)=decmat;
  retval(2)=sin2gammamat;
  retval(3)=cos2gammamat;
  return retval;

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_radec_from_altaz_actpol_c_old,args,nargout,"Convert az, el, and ctime to ra/dec and sin/cos(2gamma).\n")
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


//void get_radec_from_altaz_actpol_c(double *az, double *el, double *tvec, double *dx, double *dy, double *theta, double *ra, double *dec, int nhorns, int nt);


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_radec_from_altaz_actpol_test,args,nargout,"Convert az, el, and ctime to ra/dec and sin/cos(2gamma).\n")
{
  int nargin = args.length();
  if (nargin<3) {
    printf("Only have %d arguments in get_radec_from_altaz_actpol_c, need at least 3.\n",nargin);
    return octave_value_list();
  }
  Matrix azmat=args(0).matrix_value();
  Matrix elmat=args(1).matrix_value();
  Matrix tmat=args(2).matrix_value();

  double zero=0;
  double *dx=&zero,*dy=&zero;
  Matrix dxmat,dymat;
  int ndet=1;
  
  if (nargin>3) {
    dxmat=args(3).matrix_value();
    dx=dxmat.fortran_vec();
    dim_vector dm=dxmat.dims();
    ndet=dm(0)*dm(1);
  }
  if (nargin>4) {
    dymat=args(4).matrix_value();
    dy=dymat.fortran_vec();
    dim_vector dm=dymat.dims();
    assert(dm(0)*dm(1)==ndet);
  }

  double *theta;
  Matrix theta_mat;
  int free_theta=0;
  if (nargin>5) {
    theta_mat=args(5).matrix_value();
    dim_vector dm=theta_mat.dims();
    if (dm(0)*dm(1)==1) {
      theta=(double *)malloc(sizeof(double)*ndet);
      free_theta=1;
      for (int i=0;i<ndet;i++)
	theta[i]=theta_mat(1,1);
    }
    else {
      assert(dm(0)*dm(1)==ndet);
      theta=theta_mat.fortran_vec();
    }
  }
  else {
    theta=(double *)malloc(sizeof(double)*ndet);
    free_theta=1;
    for (int i=0;i<ndet;i++)
      theta[i]=0;
  }
  
  
  double *az=azmat.fortran_vec();
  double *el=elmat.fortran_vec();
  double *tvec=tmat.fortran_vec();
  dim_vector dm=azmat.dims();
  int nelem=dm(0)*dm(1);
  double azmin=az[0];
  double azmax=az[0];
  
  Matrix ramat(nelem,ndet);
  Matrix decmat(nelem,ndet);
  //Matrix sin2gammamat(nelem,ndet);
  //Matrix cos2gammamat(nelem,ndet);
  double *ra=ramat.fortran_vec();
  double *dec=decmat.fortran_vec();
  //double *sin2gamma=sin2gammamat.fortran_vec();
  //double *cos2gamma=cos2gammamat.fortran_vec(); 

  get_radec_from_altaz_actpol_c(az, el, tvec,dx,dy, theta, ra, dec,ndet,nelem);

#if 0

  for (int i=0;i<nelem;i++) {
    if (az[i]<azmin)
      azmin=az[i];
    if (az[i]>azmax)
      azmax=az[i];
  }
  //Fix this - JLS 23 June 2012.  Note to Nolta - for constant az, this makes no difference.
  double az_cent=(azmin+azmax)/2;
  double az_throw=(azmin-azmax)/2;


  ACTpolPointingFit *fit=(ACTpolPointingFit *)calloc(1,sizeof(ACTpolPointingFit));
  assert(fit);
  int nhorns=ndet;
  ACTpolArray *array = ACTpolArray_alloc(nhorns);
  assert(array);

  double xcent=0.0;
  double ycent=0.0;
  double freq=148.0;
  ACTpolArray_init(array, freq, xcent,ycent);
  for (int i=0;i<nhorns;i++) {
    ACTpolFeedhorn_init(&(array->horn[i]),dx[i],dy[i],theta[i]);
  }

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
      for (int j=0;j<nhorns;j++) {
	ACTpolFeedhornCoords *fc = &(coords->horn[j]);
	ra[i+j*nelem]=fc->ra;
	dec[i+j*nelem]=fc->dec;
	sin2gamma[i+j*nelem]=fc->sin2gamma;
	cos2gamma[i+j*nelem]=fc->cos2gamma;
      }
    }

#endif
  if (free_theta==1)
    free(theta);
  octave_value_list retval;
  retval(0)=ramat;
  retval(1)=decmat;
  //retval(2)=sin2gammamat;
  // retval(3)=cos2gammamat;
  return retval;

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_altaz_from_radec_ctime_c,args,nargout,"Convert ra/dec/ctime to alt/az.\n")
{
  double ra=get_value(args(0));
  double dec=get_value(args(1));
  double ct=get_value(args(2));
  double alt,az;


  int myerr=actpol_radec_to_crude_altaz(ct,ra,dec,&alt,&az);
  //assert(1==0);  //this call is missing on scinet...

  //printf("error code is %d\n",myerr);
  if (myerr) {
    fprintf(stderr,"Error in get_altaz_from_radec_ctime_c.\n");
    return octave_value_list();
  }
  octave_value_list retval;
  retval(0)=alt;
  retval(1)=az;
  return retval;
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_map_poltag,args,nargout,"Return the polarization tag for a submap inside a polarized map.  Args are (map,which_pol).\n")
{
  MAP *map=(MAP *)get_pointer(args(0));  
  int which_map=(int)get_value(args(1));
  if ((which_map<=0)||(which_map>get_npol_in_map(map))) {
    fprintf(stderr,"Requested tag %d is out of allowed range in map.  Use something between 1 and %d\n",which_map,get_npol_in_map(map));
    return octave_value_list();
  }
  int icur=0;
  for (int i=0;i<MAX_NPOL;i++) {
    if (map->pol_state[i])
      icur++;
    if (icur==which_map) {
      if (icur==1)
	return octave_value('I');
      if (icur==2)
	return octave_value('Q');
      if (icur==3)
	return octave_value('U');
      if (icur==4)
	return octave_value('V');
    }
  }
  fprintf(stderr,"Managed to not find a polarization.  Very odd...\n");
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(fit_hwp_poly_to_data_c,args,nargout,"Fit sins/cosines/(low-order) polynomials to tod data.  Args are (tod, nsin, npoly).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if ((tod->data==NULL) ||(tod->hwp==NULL)) {
    fprintf(stderr,"TOD isn't sufficiently populated in fit_hwp_poly_to_data_c.\n");
    return octave_value_list();
  }
  int npoly=1;
  int nsine=20;
  int nargin=args.length();
  if (nargin>1)
    nsine=(int)get_value(args(1));
  if (nargin>2)
    npoly=(int)get_value(args(2));
  if (1) {
    remove_hwp_poly_from_data(tod,nsine,npoly);
    return octave_value_list();
  }


  int nparam=2*nsine+npoly;
  actData **fitp=matrix(tod->ndet,nparam);
  memset(fitp[0],0,sizeof(actData)*tod->ndet*nparam);

  fit_hwp_poly_to_data(tod, nsine,npoly,fitp,NULL);
  Matrix myfitp(nparam,tod->ndet);
  //printf("doing memcpy now.\n");
  memcpy(myfitp.fortran_vec(),fitp[0],nparam*tod->ndet*sizeof(actData));
  //printf("finished.\n");
  free(fitp[0]);
  free(fitp);
  return octave_value(myfitp);


}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_sincos_mat,args,nargout,"Fill a matrix with sin/cos(n*hwp).  Args are (tod,n_sin).\n")
{
  int nterm=(int)get_value(args(1));
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  Matrix mat(2*nterm,tod->ndata);
  actData *vec=mat.fortran_vec();
  fill_tod_sin_cos_vec(tod, nterm, vec);
  return octave_value(mat);

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (test_ninkasi_linfit,args,nargout,"Test out ninkasi linear fitting.\n")
{
  int ndata=5;
  int nparam=3;
  int ncol=2;
  actData **crap=matrix(ndata,nparam);
  for (int i=0;i<ndata;i++) {
    crap[i][0]=1;
    for (int j=1;j<nparam;j++)
      crap[i][j]=crap[i][j-1]*i;
  }  
  printf("crap block is:\n %16.6e %16.6e %16.6e\n %16.6e %16.6e %16.6e\n %16.6e %16.6e %16.6e\n",crap[0][0],crap[0][1],crap[0][2],crap[1][0],crap[1][1],crap[1][2],crap[2][0],crap[2][1],crap[2][2]);
  actData **dat=matrix(ncol,ndata);
  for (int i=0;i<ncol;i++)
    for (int j=0;j<ndata;j++)
      dat[i][j]=(1-2*i)*j+(2*i);
  for (int i=0;i<ndata;i++) {
    for (int j=0;j<ncol;j++)
      printf("%4.0f ",dat[j][i]);
    printf("\n");
  }


  


  //linfit_many_vecs(NULL,crap,ndata,1,nparam);
  actData **fitp=matrix(ncol,nparam);
  linfit_many_vecs(dat,crap,ndata,ncol,nparam,fitp);
  for (int i=0;i<nparam;i++) {
    for (int j=0;j<ncol;j++) 
      printf("%14.6f ",fitp[j][i]);
    printf("\n");
  }
  return octave_value_list();

  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(get_demodulated_hwp_data_c,args,nargout,"Get (complex) I/QU data with pol data demodulated by the HWP.  Original data will be overwritten.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  actData hwp_freq=2.5;
  if (args.length()>1)
    hwp_freq=get_value(args(1));
  int nmode=get_demodulated_hwp_data(tod,hwp_freq,NULL,NULL);
  //printf("allocating.\n");
  ComplexMatrix cmat(2*nmode-1,tod->ndet);
  ComplexMatrix mat(nmode,tod->ndet);
  Complex *cvec=cmat.fortran_vec();
  Complex *vec=mat.fortran_vec();

  actComplex **mycmat=(actComplex **)malloc(sizeof(actComplex *)*tod->ndet);
  actComplex **mymat=(actComplex **)malloc(sizeof(actComplex *)*tod->ndet);


  void *crap=&(cvec[0]);
  actComplex *cvecptr=(actComplex *)crap;
  crap=&(vec[0]);
  actComplex *vecptr=(actComplex *)crap;

  //printf("got my stuff allocated.\n");

  for (int i=0;i<tod->ndet;i++) {
    mycmat[i]=cvecptr+i*(2*nmode-1);
    mymat[i]=vecptr+i*(nmode);
  }
  //printf("pointers are set.\n");
  for (int i=0;i<tod->ndet;i++) {
    memset(mycmat[i],0,(2*nmode-1)*sizeof(actComplex));
    memset(mymat[i],0,nmode*sizeof(actComplex));
    }
  //printf("cleared out some of my stuff.\n");
  get_demodulated_hwp_data(tod,hwp_freq,mymat,mycmat);
  octave_value_list retval;
  retval(0)=octave_value(mat);
  retval(1)=octave_value(cmat);
  
  return octave_value(retval);
  
}
