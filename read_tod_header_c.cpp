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
#include <getdata.h>
#include <dirfile.h>
#include <readtod.h>
#include <slalib.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif




#include <string.h>


char *get_char_from_arg(charMatrix ch);


/*--------------------------------------------------------------------------------*/


actData get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  actData myval=(actData)myptr(0,0);
  return myval;

}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_median_altaz, args, nargout, "Get median alt/az of a tod name.\n")
{
  
  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();
  
  
  
  
  char *myfroot=get_char_from_arg(args(0).char_matrix_value());  
  
  int decimate=0;
  if (nargin>1) {
    decimate=get_value(args(1));
  }
  
  mbTOD *tod=read_dirfile_tod_header_decimate(myfroot,decimate);
  
  //printf("read tod header.\n");
  
  double azmed=(double)sselect(tod->ndata/2,tod->ndata-1,tod->az);  
  double altmed=(double)sselect(tod->ndata/2,tod->ndata-1,tod->alt);  

  double ctime_med=(double)(tod->ctime+tod->deltat*(tod->ndata/2));

  //printf("got medians.\n");
  

  
  free(tod->alt);
  free(tod->az);
  free(tod->rows);
  free(tod->cols);
  free(tod);
  
  
  octave_value_list retval;
  retval(0)=azmed;
  retval(1)=altmed;
  retval(2)=ctime_med;
  return retval;

}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (alt_az_ctime2act_ra_dec_c, args, nargout, "Turn an alt, az, and ctime for ACT into an ra/dec")
{
  actData alt=get_value(args(0));
  actData az=get_value(args(1));
  actData ctime=get_value(args(2));
  Site site;
  ACTSite(&site);
  actData ra,dec;
  if (act_observed_altaz_to_mean_radec(&site,150.0,1,&ctime,&alt,&az, &ra, &dec)) {
    fprintf(stderr,"Error in alt_az_ctime2act_ra_dec_c.\n");
    return octave_value_list();
  }


  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;
  return retval;
  
    
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (read_tod_header_nopoint_c, args, nargout, "Read a TOD header, but w/out pointing.\n")
{
  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();

  char *myfroot=get_char_from_arg(args(0).char_matrix_value());
  int decimate=0;
  if (nargin>1) {
    decimate=get_value(args(1));
  }
  mbTOD *mytod=read_dirfile_tod_header_decimate(myfroot,decimate);
  mprintf(stdout,"file %s had %d detectors and %d data elements, rows and cols are %d %d.  dt=%12.4e\n",myfroot,mytod->ndet,mytod->ndata,mytod->nrow, mytod->ncol,mytod->deltat);
  mytod->pointingOffset=NULL;
  
  mprintf(stdout,"allocating cuts.\n");
  mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol);

  int64NDArray myptr(1);
  long mytod_asint=(long)mytod;
  myptr(0)=mytod_asint;


  octave_value_list retval;
  retval(0)=myptr;
  
  return retval;
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (read_tod_header_abs_nopoint_c, args, nargout, "Read a TOD header, but w/out pointing.\n")
{
  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();

  char *myfroot=get_char_from_arg(args(0).char_matrix_value());
  int decimate=0;
  if (nargin>1) {
    decimate=get_value(args(1));
  }
  mbTOD *mytod=read_dirfile_tod_header_abs(myfroot);
  mprintf(stdout,"file %s had %d detectors and %d data elements, rows and cols are %d %d.  dt=%12.4e\n",myfroot,mytod->ndet,mytod->ndata,mytod->nrow, mytod->ncol,mytod->deltat);
  mytod->pointingOffset=NULL;
  
  mprintf(stdout,"allocating cuts.\n");
  mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol);

  int64NDArray myptr(1);
  long mytod_asint=(long)mytod;
  myptr(0)=mytod_asint;


  octave_value_list retval;
  retval(0)=myptr;
  
  return retval;
}
/*--------------------------------------------------------------------------------*/





DEFUN_DLD (read_tod_header_c, args, nargout, "Read a TOD header, including pointing info etc.\n")
{
  
  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();




  char *myfroot=get_char_from_arg(args(0).char_matrix_value());  

  int decimate=0;
  if (nargin>2) {
    decimate=get_value(args(2));
  }

  mbTOD *mytod=read_dirfile_tod_header_decimate(myfroot,decimate);

  mprintf(stdout,"file %s had %d detectors and %d data elements, rows and cols are %d %d.  dt=%12.4e\n",myfroot,mytod->ndet,mytod->ndata,mytod->nrow, mytod->ncol,mytod->deltat);
  if (args(1).length()>0) { 
    if (args(1).is_sq_string ()) {
      char *pointing_file=get_char_from_arg(args(1).char_matrix_value());
      mprintf(stdout,"reading pointing offsets from %s\n",pointing_file);
      
      mytod->pointingOffset=nkReadPointingOffset(pointing_file);
      
      cut_mispointed_detectors(mytod);
      free(pointing_file);
    }
    else {
      Octave_map poff = args(1).map_value ();
      Octave_map::const_iterator myalt = poff.seek ("dalt");
      octave_value tmpalt =  poff.contents(myalt)(0);
      Octave_map::const_iterator myaz = poff.seek ("daz");
      octave_value tmpaz =  poff.contents(myaz)(0);
      Matrix dalt=tmpalt.matrix_value();
      Matrix daz=tmpaz.matrix_value();
      dim_vector dm=dalt.dims();
      dim_vector dm2=daz.dims();


      Octave_map::const_iterator myparams = poff.seek ("fitparams");
      octave_value tmppar =  poff.contents(myparams)(0);
      Matrix fitpar=tmppar.matrix_value();
      dim_vector dd_par=fitpar.dims();


      assert(dm(0)==dm2(0));
      assert(dm(1)==dm2(1));
      int fittype=0;  //maybe need to fix this later.
      int nrow=dm(0);
      int ncol=dm(1);
      mbPointingOffset *offset = nkPointingOffsetAlloc(nrow,ncol,fittype);
      assert(dd_par(0)*dd_par(1)==offset->nparam);
      double *ptr=fitpar.fortran_vec();
      for (int i=0;i<offset->nparam;i++)
	offset->fit[i]=ptr[i];
      
      for (int i=0;i<nrow;i++)
	for (int j=0;j<ncol;j++) {
	  offset->offsetAlt[i][j]=dalt(i,j);
	  offset->offsetAzCosAlt[i][j]=daz(i,j);
	}
      mytod->pointingOffset=offset;
    }
  }
  else
    mytod->pointingOffset=NULL;
  
  mprintf(stdout,"allocating cuts.\n");
  mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol);
  cut_mispointed_detectors(mytod);
  
  assign_tod_ra_dec(mytod);
  //pivot finding lives here.
  find_pointing_pivots(mytod,0.5);
  find_tod_radec_lims(mytod); 

#ifdef ACTDATA_DOUBLE
  //printf("reading double %ld.\n",sizeof(actData));
#else
  //printf("reading single %ld\n",sizeof(actData));
#endif

  //#ifdef ACTDATA_DOUBLE
  //printf("Ra/dec lims are %14.6lf %14.6lf %14.6lf %14.6lf\n",mytod->ramin,mytod->ramax,mytod->decmin,mytod->decmax);
  //#else
  mprintf(stdout,"Ra/dec lims are %14.6f %14.6f %14.6f %14.6f\n",mytod->ramin,mytod->ramax,mytod->decmin,mytod->decmax);
  //#endif
  int64NDArray myptr(1);
  long mytod_asint=(long)mytod;
  myptr(0)=mytod_asint;
  Matrix lims(4,1);
  lims(0,0)=mytod->ramin;
  lims(1,0)=mytod->ramax;
  lims(2,0)=mytod->decmin;
  lims(3,0)=mytod->decmax;

  //Array retlen(2,3,4,5);



  octave_value_list retval;
  retval(0)=myptr;
  retval(1)=lims;
  
  return retval;
}       


/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
char *get_char_from_arg(charMatrix ch)
{
  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}





#if 0  
  
  cout << "Hello.\n";
  cout << ch.row_as_string(0) << "\n";
  cout << "length is " << ch.length() << "\n";

  int nn=ch.length();
  char *s=ch.fortran_vec();

  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  printf("\n\nThis is still a test %s.\n",ss);



  mbTOD *tod=read_dirfile_tod_header(s );
  int n=tod->ndata;
  Matrix a(n,1);
  double *tmp=a.fortran_vec();
  for (int i=0;i<n;i++)
    tmp[i]=tod->az[i];
  
  
  
  printf("In octave, have %d elements.\n",n);
  
  
  

  return octave_value(a);
#endif

/*--------------------------------------------------------------------------------*/
template <typename T> ColumnVector vec2col(T *vec_in, int n)
{
  ColumnVector vec(n);
  double *vv=vec.fortran_vec();
  for (int i=0;i<n;i++)
    vv[i]=(double)vec_in[i];

  return vec;
    
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (read_dirfile_channel, args, nargout, "Read a dirfile channel.  Allowable types are usUSfd\n")
{

  if (args.length()<3) {
    fprintf(stderr,"Error - need at least three arguments (filename, fieldname, type) in read_dirfile_channel.\n");
    return octave_value_list();
  }

  charMatrix filename=args(0).char_matrix_value();
  char *fname=get_char_from_arg(filename);
  //printf("Filename is .%s.\n",fname);

  int n=0,status=0; 

  charMatrix fieldname=args(1).char_matrix_value();
  char *ffname=get_char_from_arg(fieldname);
  charMatrix ntype=args(2).char_matrix_value();  
  char mytype=ntype(0);
  //printf("mytype is %c\n",mytype);
  

  void *data=dirfile_read_channel_direct(mytype,fname,ffname,&n);
  free(ffname);
  free(fname);
  if (data==NULL) {
    return octave_value_list();
  }
  
  ColumnVector dd;
  switch(mytype) {
  case 's':
    dd=vec2col((int16_t *)data,n);
    break;
  case 'u':
    dd=vec2col((uint16_t *)data,n);
    break;
  case 'f':
    dd=vec2col((float *)data,n);
    break;
  case 'd':
    dd=vec2col((double *)data,n);
    break;
  case 'U':
    dd=vec2col((uint32_t *)data,n);
    break;
  case 'S':
    dd=vec2col((int32_t *)data,n);
    break;
    
  default:
    fprintf(stderr,"type not recognized in switch.  A bit odd...\n");
  }

  return octave_value(dd);
  
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (read_many_dirfile_channels_c, args, nargout, "Read many dirfile channels.\n")
{

  if (args.length()<2) {
    fprintf(stderr,"Error - need at least three arguments (filename, fieldname) in read_dirfile_channel.\n");
    return octave_value_list();
  }

  int nfields=args.length();

  octave_value_list retval;

  charMatrix filename=args(0).char_matrix_value();  
  char *fname=get_char_from_arg(filename);
  
  int status;
  struct FormatType *format = GetFormat( fname, NULL, &status );
  if (format==NULL) {
    fprintf(stderr,"Warning - problem reading format from file .%s.\n",fname);
    return octave_value_list();
  }
  
  for (int i=1;i<nfields;i++) {
    charMatrix channame=args(i).char_matrix_value();
    char *cname=get_char_from_arg(channame);
    int nsamp=0;
    double *crap=dirfile_read_double_channel(format,cname,&nsamp);
    ColumnVector vec(nsamp);
    double *vv=vec.fortran_vec();
    memcpy(vv,crap,nsamp*sizeof(double));

    free(crap);
    free(cname);
    retval(i-1)=vec;    
  }

  GetDataClose(format);
  free(format);
  
   
  return retval;
}

/*--------------------------------------------------------------------------------*/
