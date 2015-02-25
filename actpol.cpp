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
#include <omp.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

#define ACTPOL_NEW


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
#ifdef ACTPOL_NEW
  ACTpolArrayCoords_init(coords,ACTPOL_COORDSYS_RA_DEC);
#else
  ACTpolArrayCoords_init(coords);
#endif

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
  //printf("alt/az limits are %14.6f %14.6f %14.6f %14.6f\n",altmin,altmax,azmin,azmax);

  
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
#ifdef ACTPOL_NEW
    ACTpolArrayCoords_update(coords, state);
#else
    ACTpolArrayCoords_update_fast(coords, state);
#endif
    double *ra=raptr+isamp*(long)ndet;
    double *dec=decptr+isamp*(long)ndet;
    for (int k=0;k<ndet;k++) {
      ACTpolFeedhornCoords *fc = coords->horn + k;
#ifdef ACTPOL_NEW
      ra[k]=fc->a;
      //dec[k]=acos(fc->b);
      dec[k]=(fc->b);
#else
      ra[k]=fc->ra;
      dec[k]=fc->sindec;
#endif
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

  if (args.length()==1)    {
    ACTpolPointingFit *fit=tod->actpol_pointing;
    if (!fit) {
      fprintf(stderr,"Unable to update pointing as it doesn't exist, and you haven't given me any more information.\n");
      return octave_value_list();

    }
    update_actpol_pointing(tod,NULL,NULL,NULL,150,1);
    return octave_value_list();
    
  }

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
DEFUN_DLD(destroy_actpol_pointing,args,nargout,"Destroy an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (tod->actpol_pointing) {
    ACTpolPointingFit *fit=tod->actpol_pointing;
    free(fit->ipiv);
    ACTpolArray_free(fit->array);
    free(fit->dx);
    free(fit->dy);
    free(fit->theta);
    free(fit);
    tod->actpol_pointing=NULL;
  }
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(destroy_actpol_pointing_fit,args,nargout,"Destroy an actpol pointing fit structure.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (tod->ACTPol_pointing_fit) {
    ACTPolPointingFit2 *fit=tod->ACTPol_pointing_fit;
    free(fit->az_scale);
    free(fit->alt_scale);
    free(fit->t_scale);
    free(fit->ra_fitp[0][0]);
    free(fit->ra_fitp[0]);
    free(fit->ra_fitp);
    free(fit->dec_fitp[0][0]);
    free(fit->dec_fitp[0]);
    free(fit->dec_fitp);
    free(fit);
    tod->ACTPol_pointing_fit=NULL;
  }
  else
    fprintf(stderr,"Warning - fit not found in destroy_actpol_pointing.\n");
  return octave_value_list();


}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(initialize_actpol_pointing_fit_c,args,nargout,"Set up an actpol pointing fit.  Args are tod, (az,alt,t) scaled vecs, ra fit parameters, dec fit parameters.\n")
{

  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  ColumnVector az=args(1).column_vector_value();
  ColumnVector alt=args(2).column_vector_value();
  ColumnVector t=args(3).column_vector_value();
  Matrix ra_fit=args(4).matrix_value();
  Matrix dec_fit=args(5).matrix_value();

  actData *azptr=az.fortran_vec();
  actData *altptr=alt.fortran_vec();
  actData *tptr=t.fortran_vec();
  actData *raptr=ra_fit.fortran_vec();
  actData *decptr=dec_fit.fortran_vec();

  initialize_actpol_pointing_fit(tod,azptr,altptr,tptr,raptr,decptr);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(evaluate_actpol_pointing_fit,args,nargout,"Evaluate an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  evaluate_actpol_pointing_fit(tod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (does_tod_have_twogamma_fit,args,nargout,"Check to see if a tod already has a twogamma fit saved.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (!tod->actpol_pointing) 
    return octave_value(0);
  if (tod->actpol_pointing->gamma_az_sin_coeffs) 
    return octave_value(1);
  return octave_value(0);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(set_tod_twogamma_fit_c,args,nargout,"Install a saved twogamma fit to az.\n")
{
  if (args.length()<3) {
    printf("need at least 3 args to set_tod_twogamma_fit_c.\n");
    return octave_value_list();
  }
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (!tod->actpol_pointing) {
    printf("actpol_pointing not saved in TOD.\n");
    return octave_value_list();
  }

  if (!tod->actpol_pointing) {
    printf("missing pointing structure in set_tod_twogamma_fit_c, aborting.\n");
    return octave_value_list();
  }
  if (tod->actpol_pointing->gamma_az_sin_coeffs) {
    printf("found twogamma fit in set_tod_twogamma_fit_c already, so skipping.\n");
    return octave_value_list();
  }

  Matrix sinmat=args(1).matrix_value();
  Matrix cosmat=args(2).matrix_value();
  dim_vector dm_sin=sinmat.dims();
  dim_vector dm_cos=cosmat.dims();
  //printf("dims are %d %d, %d %d\n",dm_sin(0),dm_sin(1),dm_cos(0),dm_cos(1));
  if ((dm_sin(0)!=dm_cos(0))||(dm_sin(1)!=dm_cos(1))||(dm_sin(1)!=tod->ndet)) {
    printf("dimesion mismatch in set_tod_twogamma_fit_c.\n");
    return octave_value_list();
  }
  ACTpolPointingFit *pfit=tod->actpol_pointing;
  int nterm=dm_sin(0)-1;
  pfit->gamma_az_sin_coeffs=matrix(tod->ndet,nterm);
  pfit->gamma_az_cos_coeffs=matrix(tod->ndet,nterm);
  pfit->gamma_ctime_sin_coeffs=(actData *)malloc(sizeof(actData)*tod->ndet);
  pfit->gamma_ctime_cos_coeffs=(actData *)malloc(sizeof(actData)*tod->ndet);
  pfit->n_gamma_az_coeffs=nterm;
  for (int det=0;det<tod->ndet;det++) {
    for (int j=0;j<nterm;j++) {
      pfit->gamma_az_sin_coeffs[det][j]=sinmat(nterm-1-j,det);
      pfit->gamma_az_cos_coeffs[det][j]=cosmat(nterm-1-j,det);
    }
    pfit->gamma_ctime_sin_coeffs[det]=sinmat(nterm,det);
    pfit->gamma_ctime_cos_coeffs[det]=cosmat(nterm,det);
  }
  //  for (int i=0;i<nterm;i++) {
  // printf("first detector fits are %d %14.6g %14.6g\n",i,pfit->gamma_az_sin_coeffs[0][i],pfit->gamma_az_cos_coeffs[0][i]);
  //}
  
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(precalc_actpol_pointing_exact_subsampled_c,args,nargout,"Calculate a subsample of the array pointing.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  int fac=(int)get_value(args(1));
  int len=(tod->ndata+fac-1)/fac;
  //printf("len is %d %d %d\n",len,tod->ndata,fac);
  
  dim_vector dm(len,tod->ndet);
  Matrix ra(dm);
  Matrix dec(dm);
  Matrix twogamma(dm);

  double *ravec=ra.fortran_vec();
  double *decvec=dec.fortran_vec();
  double *twogammavec=twogamma.fortran_vec();
  
  double **ramat=(double **)malloc(sizeof(double *)*tod->ndata);
  double **decmat=(double **)malloc(sizeof(double *)*tod->ndata);
  double **twogammamat=(double **)malloc(sizeof(double *)*tod->ndata);

  

  for (int i=0;i<tod->ndet;i++) {
    ramat[i]=ravec+i*len;
    decmat[i]=decvec+i*len;
    twogammamat[i]=twogammavec+i*len;
  }


  precalc_actpol_pointing_exact_subsampled(tod,fac,ramat,decmat,twogammamat);
  free(ramat);
  free(decmat);
  free(twogammamat);
  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;
  retval(2)=twogamma;

  return retval;


}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD(precalc_actpol_pointing_exact,args,nargout,"Precalculated stuff for an actpol pointing fit.\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (tod->ACTPol_pointing_fit) {
    evaluate_actpol_pointing_fit(tod);
    return octave_value_list();
  }
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
#ifdef ACTPOL_NEW
  ACTpolArrayCoords_init(coords,ACTPOL_COORDSYS_RA_DEC);
#else
  ACTpolArrayCoords_init(coords);
#endif

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

#ifdef ACTPOL_NEW
	ra[i+j*nelem]=fc->a;
	//dec[i+j*nelem]=acos(fc->b);
	dec[i+j*nelem]=(fc->b);
#else
	ra[i+j*nelem]=fc->ra;
	dec[i+j*nelem]=fc->dec;
#endif
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
#ifdef ACTPOL_NEW
  ACTpolArrayCoords_init(coords,ACTPOL_COORDSYS_RA_DEC);
#else
  ACTpolArrayCoords_init(coords);
#endif

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
#ifdef ACTPOL_NEW
      ra[i]=fc->a;
      //dec[i]=acos(fc->b);
      dec[i]=(fc->b);
#else
      ra[i]=fc->ra;
      dec[i]=fc->dec;
#endif
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
#ifdef ACTPOL_NEW
  ACTpolArrayCoords_init(coords,ACTPOL_COORDSYS_RA_DEC);
#else
  ACTpolArrayCoords_init(coords);
#endif

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
#ifdef ACTPOL_NEW
	ra[i+j*nelem]=fc->a;
	//dec[i+j*nelem]=acos(fc->b);
	dec[i+j*nelem]=(fc->b);
#else
	ra[i+j*nelem]=fc->ra;
	dec[i+j*nelem]=fc->dec;
#endif
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
#if 0
  int myerr=actpol_radec_to_crude_altaz(ct,ra,dec,&alt,&az);
#else
  int myerr=99999;
  assert(1==0);  //this call is missing on scinet...  
#endif

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
	return octave_value("QQ");
      if (icur==5)
	return octave_value("QU");
      if (icur==6)
	return octave_value("UU");
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
  if (nargout==0) {
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
DEFUN_DLD(fit_hwp_az_poly_to_data_c,args,nargout,"Fit sins/cosines/(low-order) az/polynomials to tod data.  Args are (tod, nsin, naz,npoly).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  //if ((tod->data==NULL) ||(tod->hwp==NULL)) {
  if ((tod->data==NULL)) {
    fprintf(stderr,"TOD isn't sufficiently populated in fit_hwp_poly_to_data_c.\n");
    return octave_value_list();
  }
  int npoly=1;
  int naz=1;
  int nsine=20;
  int nargin=args.length();
  if (nargin>1)
    nsine=(int)get_value(args(1));
  if (nargin>2)
    naz=(int)get_value(args(2));
  if (nargin>3)
    npoly=(int)get_value(args(3));
  if (nargout==0) {
    remove_hwp_az_poly_from_data(tod,nsine,naz,npoly);
    return octave_value_list();
  }


  int nparam=2*nsine+naz+npoly;
  actData **fitp=matrix(tod->ndet,nparam);
  memset(fitp[0],0,sizeof(actData)*tod->ndet*nparam);

  fit_hwp_az_poly_to_data(tod, nsine,naz,npoly,fitp,NULL);
  Matrix myfitp(nparam,tod->ndet);
  //printf("doing memcpy now.\n");
  memcpy(myfitp.fortran_vec(),fitp[0],nparam*tod->ndet*sizeof(actData));
  //printf("finished.\n");
  free(fitp[0]);
  free(fitp);
  return octave_value(myfitp);


}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (print_detector_offsets_actpol,args,nargout,"Print the detector offsets, angle, etc.  Args are (tod,row,col).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  int myrow=(int)get_value(args(1));
  int mycol=(int)get_value(args(2));
  int det=tod->dets[myrow][mycol];
  printf("detector number is %d\n",det);
  printf("Offsets are %g %g with angle %g\n",tod->actpol_pointing->dx[det] ,tod->actpol_pointing->dy[det], tod->actpol_pointing->theta[det]);
  if (tod->actpol_pointing->gamma_ctime_cos_coeffs) {
    printf("ctime coefficients are %g %g\n",tod->actpol_pointing->gamma_ctime_cos_coeffs[det],tod->actpol_pointing->gamma_ctime_sin_coeffs[det]);
    for (int i=0;i<tod->actpol_pointing->n_gamma_az_coeffs;i++) {
      printf("fit coeffs %d are %g %g\n",i,tod->actpol_pointing->gamma_az_cos_coeffs[det][i],tod->actpol_pointing->gamma_az_sin_coeffs[det][i]);
    }
  }
  return octave_value_list();
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
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(remodulate_hwp_data_c,args,nargout,"Turn demodulated T/pol data back into a timestream and write into the TOD.  Args are (tod,T,pol,[hwp frequency]) as made by get_demodulated_hwp_data_c.\n")
{

  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  ComplexMatrix mat=args(1).complex_matrix_value();
  ComplexMatrix cmat=args(2).complex_matrix_value();


  actData hwp_freq=2.5;
  if (args.length()>3)
    hwp_freq=get_value(args(3));

  
  int nmode=get_demodulated_hwp_data(tod,hwp_freq,NULL,NULL);

  
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
  //printf("cleared out some of my stuff.\n");
  remodulate_hwp_data(tod,hwp_freq,mymat,mycmat);
  octave_value_list retval;
  //retval(0)=octave_value(mat);
  //retval(1)=octave_value(cmat);
  
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (init_demod_data_c, args, nargout, "Set up a demod data structure and return it.  Args are tod,freqs,[hwp_freq,bottom_freq,top_freq]\n")
{
  
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  Matrix freqs(2,1);
  freqs(0,0)=2;
  freqs(1,0)=4;
  if (args.length()>=2)
    freqs=args(1).matrix_value();
  dim_vector dm=freqs.dims();
  int nfreq=dm(0)*dm(1);
  actData hwp_freq=0;
  if (args.length()>=3)
    hwp_freq=get_value(args(2));
  actData low_pass=2;
  if (args.length()>=4)
    low_pass=get_value(args(3));
  actData high_pass=0;
  if (args.length()>=5)
    high_pass=get_value(args(4));

  DemodData *demod=init_demod_data(tod,hwp_freq,low_pass,0,high_pass,0);
  set_demod_freqs(demod,freqs.fortran_vec(),nfreq);

  int64NDArray myptr(1);

  myptr(0)=(long)demod;
  tod->demod=demod;
  
  return octave_value(myptr);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD(demodulate_data_c,args,nargout,"Demodulate data using HWP and a demod structure.  Args are (tod,demod).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  if (demod==NULL)
    demod=(DemodData *)get_pointer(args(1));
  demodulate_data(tod,demod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD(remodulate_data_c,args,nargout,"Remodulate data using HWP and a demod structure.  Args are (tod,demod).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  if (demod==NULL)
    demod=(DemodData *)get_pointer(args(1));
  remodulate_data(tod,demod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(write_demodulated_detector_text,args,nargout,"Write a text file with demodulated data in it.\n")
{
  DemodData *demod=(DemodData *)get_pointer(args(0));
  int det=0;
  if (args.length()>1)
    det=(int)get_value(args(1));

  char fname[512];
  //sprintf(fname,"demodulate_streams_r%02dc%02d.txt",tod->rows[det],tod->cols[det]);
  FILE *outfile=fopen("demodulated_streams.txt","w");
  //FILE *outfile=fopen(fname,"w");
  int ncol=get_demod_nchannel(demod);
  for (int i=0;i<demod->nmode;i++) {
    for (int j=0;j<ncol;j++)  {
      //fprintf(outfile,"%16.7e %16.7e   ",creal(demod->data[j+det*ncol][i]),cimag(demod->data[j+det*ncol][i]));      
      fprintf(outfile,"%16.7e %16.7e   ",(demod->data[j+det*ncol][i][0]),(demod->data[j+det*ncol][i][1]));      
    }
    fprintf(outfile,"\n");
  }
  fclose(outfile);

  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(return_data_from_demod_c,args,nargout,"Return the demodulated data stored in a demod C structure.  Args are (tod,demod).\n")
{
  if (args.length()<1) {
    printf("need at least (tod) inputs in return_data_from_demod_c.\n");
    return octave_value_list();
  }
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  if ((demod==NULL)||(args.length()>1))
    demod=(DemodData *)get_pointer(args(1));
  if (demod->data==NULL) {
    printf("do not have demodulated data stored in return_data_from_demod_c.\n");
    return octave_value_list();
  }
  int ncol=get_demod_nchannel(demod);
  ComplexMatrix cmat(demod->nmode,ncol*tod->ndet);
  //memcpy(cmat.fortran_vec(),&(demod->data[0][0][0]),ncol*demod->nmode*tod->ndet*sizeof(actComplex));
  memcpy(cmat.fortran_vec(),(demod->data[0]),ncol*demod->nmode*tod->ndet*sizeof(actComplex));
  return octave_value(cmat);

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(push_data_to_demod_c,args,nargout,"Push demodulated data into the demod structure.  Args are (tod,data).\n")
{
  if (args.length()<2) {
    printf("need at least (tod,data) inputs in push_data_to_demod_c.\n");
    return octave_value_list();
  }
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  if (demod==NULL) {
    fprintf(stderr,"Demod data not present in push_data_to_demod_c.\n");
    return octave_value_list();
  }
  
  ComplexMatrix cmat=args(1).complex_matrix_value();
  if (demod->data==NULL) {
    printf("do not have demodulated data stored in return_data_from_demod_c.\n");
    return octave_value_list();
  }
  int ncol=get_demod_nchannel(demod);
  //ComplexMatrix cmat(demod->nmode,ncol*tod->ndet);
  //memcpy(cmat.fortran_vec(),&(demod->data[0][0][0]),ncol*demod->nmode*tod->ndet*sizeof(actComplex));
  memcpy(&(demod->data[0][0][0]),cmat.fortran_vec(),ncol*demod->nmode*tod->ndet*sizeof(actComplex));
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(free_demod_c,args,nargout,"Free data stored in a demod structure.  Args are (tod).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  free_demod_data(tod->demod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(destroy_demod_c,args,nargout,"Destroy data stored in a demod structure.  Args are (tod).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  DemodData *demod=tod->demod;
  destroy_demod_data(tod->demod);
  tod->demod=NULL;
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(does_tod_have_demod,args,nargout,"Report if a TOD has a demod struture in it.  Arg is (tod).\n")
{
  mbTOD *tod=(mbTOD *)get_pointer(args(0));
  if (tod->demod)
    return octave_value(1);
  else return octave_value(0);
  
  return octave_value_list();  //never get here.
}
