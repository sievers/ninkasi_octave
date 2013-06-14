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
char *get_char_from_arg(charMatrix ch)
{
  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}




/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (allocate_tod_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{
  mbTOD *mytod=(mbTOD *)malloc(sizeof(mbTOD));
  memset(mytod,0,sizeof(mbTOD));
  int64NDArray myptr(1);
  
  myptr(0)=(long)mytod;
  //long *ptr=myptr.fortran_vec();
  //ptr[0]=(long)mytod;
  return octave_value(myptr);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_dt_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData dt=get_value(args(1));
  mytod->deltat=dt;
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_timevec_c, args, nargout, "Store a full time vector, useful for unevenly sampled data.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix tvec=args(1).matrix_value();
  if (mytod->dt)
    free(mytod->dt);
  mytod->dt=(double *)malloc(sizeof(double)*mytod->ndata);
  double *tt=tvec.fortran_vec();
  for (int i=0;i<mytod->ndata;i++)
    mytod->dt[i]=tt[i];
  printf("first time sample is %16.7e %16.7e %d\n",mytod->dt[0],tt[0],mytod->ndata);
  mytod->ctime=tt[0];

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_ndata_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int ndata=(int)get_value(args(1));
  mytod->ndata=ndata;
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_pointing_saved, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix ra=args(1).matrix_value();
  Matrix dec=args(2).matrix_value();
  
  dim_vector dm_ra=ra.dims();
  dim_vector dm_dec=dec.dims();

  printf("dm_ra is %d %d\n",dm_ra(0),dm_ra(1));

  mytod->ra_saved=matrix(mytod->ndet,mytod->ndata);
  mytod->dec_saved=matrix(mytod->ndet,mytod->ndata);

  memcpy(mytod->ra_saved[0],ra.fortran_vec(),mytod->ndet*mytod->ndata*sizeof(actData));
  memcpy(mytod->dec_saved[0],dec.fortran_vec(),mytod->ndet*mytod->ndata*sizeof(actData));

  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_tod_pointing_saved_branch,args,nargout,"Set the branch cut on a TOD pointing.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->ra_saved) {
    printf("Do not have ra saved in set_tod_pointing_saved_branch.\n");
    return octave_value_list();
  }
  double val=get_value(args(1));
  for (int i=0;i<mytod->ndet;i++)
    for (int j=0;j<mytod->ndata;j++)
      if (mytod->ra_saved[i][j]<val)
	mytod->ra_saved[i][j]+=2*3.141592653589793;
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (free_tod_ra_saved, args, nargout, "Free cached RA.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->ra_saved) {
    free(mytod->ra_saved[0]);
    free(mytod->ra_saved);
    mytod->ra_saved=NULL;
  }
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (free_tod_dec_saved, args, nargout, "Free cached DEC.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->dec_saved) {
    free(mytod->dec_saved[0]);
    free(mytod->dec_saved);
    mytod->dec_saved=NULL;
  }
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (free_tod_pointing_saved, args, nargout, "Read a TOD header, including pointing info etc.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->ra_saved) {
    free(mytod->ra_saved[0]);
    free(mytod->ra_saved);
    mytod->ra_saved=NULL;
  }
  if (mytod->dec_saved) {
    free(mytod->dec_saved[0]);
    free(mytod->dec_saved);
    mytod->dec_saved=NULL;
  }
  
 

#ifdef ACTPOL
  int free_twogamma=1;
  if (args.length()>1) {
    free_twogamma=(int)get_value(args(1));
  }
  if (free_twogamma) {
    if (mytod->twogamma_saved) {
      free(mytod->twogamma_saved[0]);
      free(mytod->twogamma_saved);
      mytod->twogamma_saved=NULL;
    }
  }
  else
    printf("Skipping two_gamma free.\n");
  
  
#endif
  
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (free_tod_2gamma_saved,args,nargou,"Free 2gamma in a tod.\n")
{
#ifdef ACTPOL
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->twogamma_saved) {
    free(mytod->twogamma_saved[0]);
    free(mytod->twogamma_saved);
    mytod->twogamma_saved=NULL;
  }
#endif
  return octave_value_list();
  
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (free_tod_data_saved, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->data_saved) {
    free(mytod->data_saved[0]);
    free(mytod->data_saved);
    mytod->data_saved=NULL;
  }
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_data_saved, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix data=args(1).matrix_value();
  
  dim_vector dm=data.dims();


  mytod->data_saved=matrix(mytod->ndet,mytod->ndata);

  memcpy(mytod->data_saved[0],data.fortran_vec(),mytod->ndet*mytod->ndata*sizeof(actData));

  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (update_tod_data_saved, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix data=args(1).matrix_value();
  
  dim_vector dm=data.dims();


  //mytod->data_saved=matrix(mytod->ndet,mytod->ndata);

  memcpy(mytod->data_saved[0],data.fortran_vec(),mytod->ndet*mytod->ndata*sizeof(actData));

  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_rowcol_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix rows=args(1).matrix_value();
  Matrix cols=args(2).matrix_value();
  dim_vector dm=rows.dims();

  
  int ndet=dm(0)*dm(1);
  mytod->ndet=ndet;
  printf("setting ndet to %d\n",ndet);
  mytod->rows=(int *)malloc(sizeof(int)*ndet);
  mytod->cols=(int *)malloc(sizeof(int)*ndet);
  double *rowptr=rows.fortran_vec();
  double *colptr=cols.fortran_vec();
  int maxrow=0;
  int maxcol=0;

  for (int i=0;i<ndet;i++) {
    mytod->rows[i]=rowptr[i];
    mytod->cols[i]=colptr[i];
    if (mytod->rows[i]>maxrow)
      maxrow=mytod->rows[i];
    if (mytod->cols[i]>maxcol)
      maxcol=mytod->cols[i];
  }
  mytod->nrow=maxrow+1;
  mytod->ncol=maxcol+1;
  printf("setting nrow/ncol to %d %d\n",mytod->nrow,mytod->ncol);

  mytod->dets=imatrix(mytod->nrow,mytod->ncol);
  
  for (int i=0;i<mytod->ndet;i++)
    mytod->dets[mytod->rows[i]][mytod->cols[i]]=i;
  
  
  return octave_value_list();

}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (alloc_tod_cuts_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_radec_lims_c, args, nargout, "Return the tod pointing limits.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (args.length()>4) {
    double ramin=get_value(args(1));
    double ramax=get_value(args(2));
    double decmin=get_value(args(3));
    double decmax=get_value(args(4));
    mytod->ramin=ramin;
    mytod->ramax=ramax;
    mytod->decmin=decmin;
    mytod->decmax=decmax;
    return octave_value_list();
  }

  double ramin,ramax,decmin,decmax;
  //printf("getting ready to allocate.\n");
#pragma omp parallel shared(mytod,ramin,ramax,decmin,decmax)  default(none) 
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(mytod);
    //printf("allocated.\n");
    ramin=1000;
    ramax=-1000;
    decmin=1000;
    decmax=-1000;
    double myramin=ramin;
    double myramax=ramax;
    double mydecmin=decmin;
    double mydecmax=decmax;
    //printf("checking ndet.\n");
    //printf("it is %d\n",mytod->ndet);
#pragma omp for
    for (int det=0;det<mytod->ndet;det++) {      
      get_radec_from_altaz_fit_1det_coarse(mytod,det,scratch);
      for (int i=0;i<mytod->ndata;i++) {
	if (scratch->ra[i]<myramin)
	  myramin=scratch->ra[i];
	if (scratch->dec[i]<mydecmin)
	  mydecmin=scratch->dec[i];
	if (scratch->ra[i]>myramax)
	  myramax=scratch->ra[i];
	if (scratch->dec[i]>mydecmax)
	  mydecmax=scratch->dec[i];	
      }
    }
    destroy_pointing_fit_scratch(scratch);

#pragma omp critical
    {
      if (myramin<ramin)
	ramin=myramin;
      if (myramax>ramax)
	ramax=myramax;
      if (mydecmin<decmin)
	decmin=mydecmin;
      if (mydecmax>decmax)
	decmax=mydecmax;

    }
  }
  //printf("ramin/ramax decmin/decmax are %14.6f %14.6f %14.6f %14.6f\n",ramin,ramax,decmin,decmax);
  mytod->ramin=ramin;
  mytod->ramax=ramax;
  mytod->decmin=decmin;
  mytod->decmax=decmax;

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
#ifdef ACTPOL
DEFUN_DLD (set_tod_hwp_angle_c,args,nargout,"Set the HWP angle of a TOD.  Args are (tod,hwp).\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix hwp_mat=args(1).matrix_value();
  double *hwp=hwp_mat.fortran_vec();
  
  tod->hwp=(actData *)malloc(sizeof(actData)*tod->ndata);
  for (int i=0;i<tod->ndata;i++)
    tod->hwp[i]=hwp[i];
  
  return octave_value_list();  

}
#endif

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (cuts_extend_c, args, nargout, "Extend a cut.  args are (tod,start,stop,row,col).\n")
{
  if (args.length() <5) {
    fprintf(stderr,"You need to pass at least 5 arguments to cuts_extend_c.\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int start=(int)get_value(args(1));
  int stop=(int)get_value(args(2));
  int row=(int)get_value(args(3));
  int col=(int)get_value(args(4));

  mbCutsExtend(mytod->cuts,start,stop,row,col);

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_altaz_c, args, nargout, "Write alt/az into a TOD.  args are (tod,alt,az).\n")
{
  int nargin = args.length();
  if (nargin<3) {
    printf("not enough arguments to set_tod_altaz_c.  Need at least 3.\n");
    return octave_value_list();
  }
  int force_alloc=0;
  if (nargin>3) {
    int val=(int)get_value(args(3));
    if (val)
      force_alloc=1;
  }

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix alt=args(1).matrix_value();
  Matrix az=args(2).matrix_value();

  if ((mytod->alt==NULL)||(force_alloc))
    mytod->alt=(actData *)malloc(sizeof(actData)*mytod->ndata);
  if ((mytod->az==NULL)||(force_alloc))
    mytod->az=(actData *)malloc(sizeof(actData)*mytod->ndata);

  memcpy(mytod->alt,alt.fortran_vec(),mytod->ndata*sizeof(actData));
  memcpy(mytod->az,az.fortran_vec(),mytod->ndata*sizeof(actData));

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_tod_filename_c,args,nargout,"Set a TODs dirfile name.  args are (tod,fname).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  charMatrix fname=args(1).char_matrix_value();
  
  mytod->dirfile=get_char_from_arg(fname);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(free_tod_altaz_c,args,nargout,"Free the alt/az from a TOD.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (tod->alt) {
    free(tod->alt);
    tod->alt=NULL;
  }
  if (tod->az) {
    free(tod->az);
    tod->az=NULL;
  }
  return octave_value_list();
    
    
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(free_tod_timevec_c,args,nargout,"Free the time from a TOD.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (tod->dt) {
    free(tod->dt);
    tod->dt=NULL;
  }
  return octave_value_list();    
    
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(free_tod_hwp_c,args,nargout,"Free the HWP from a TOD.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (tod->hwp) {
    free(tod->hwp);
    tod->hwp=NULL;
  }
  return octave_value_list();    
    
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD(free_tod_pointer_c,args,nargout,"Free a TOD pointer.  Anything still saved in here will be lost!\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (tod)
    free(tod);
  return octave_value_list();
}
