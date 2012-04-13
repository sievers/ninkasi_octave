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
  mytod->dt=(float *)malloc(sizeof(float)*mytod->ndata);
  double *tt=tvec.fortran_vec();
  for (int i=0;i<mytod->ndata;i++)
    mytod->dt[i]=tt[i];
  
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


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (cuts_extend_c, args, nargout, "Extend a cut.  args are (tod,start,stop,row,col).\n")
{
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
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix alt=args(1).matrix_value();
  Matrix az=args(2).matrix_value();
  mytod->alt=(actData *)malloc(sizeof(actData)*mytod->ndata);
  mytod->az=(actData *)malloc(sizeof(actData)*mytod->ndata);

  memcpy(mytod->alt,alt.fortran_vec(),mytod->ndata*sizeof(actData));
  memcpy(mytod->az,az.fortran_vec(),mytod->ndata*sizeof(actData));

  return octave_value_list();
}
