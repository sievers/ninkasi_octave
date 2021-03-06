#include <octave/oct.h>
#include <octave/Cell.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <getdata_v2.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <getdata.h>


  //#include "/home/nolta/local/include/getdata.h"

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

DEFUN_DLD (init_getdata_file_c, args, nargout, "Open a compressed dirfile.\n")
{
  
  char *myfroot=get_char_from_arg(args(0).char_matrix_value());
  const DIRFILE *file=gd_open(myfroot,0);
  int myerr=gd_error(file);
  if (myerr) {
    fprintf(stderr,"Error opening dirfile %s for reading.\n",myfroot);
    return octave_value_list();
  }
  int64NDArray myptr(1);
  myptr(0)=(long)file;
  return octave_value(myptr);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (close_getdata_file_c, args, nargout, "Open a compressed dirfile.\n")
{
  DIRFILE *file=( DIRFILE *)get_pointer(args(0));
  gd_close(file);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (getdata_double_channel_c,args,nargout,"Read a dirfile channel as a double.  Args are (dirfile,channel name).\n")
{
  DIRFILE *file=( DIRFILE *)get_pointer(args(0));
  char *myfield=get_char_from_arg(args(1).char_matrix_value());
  int nframes=gd_nframes(file);
  //printf("nframes is %d\n",nframes);

  int spf=gd_spf(file,myfield);
  int myerr=gd_error(file);
  if (myerr!=GD_E_OK) {
    fprintf(stderr,"Error reading samples per frame on channel %s\n",myfield);
    return octave_value(myerr);
  }
  //printf("samples per frame is %d\n",spf);

  Matrix data(spf*(nframes+1),1);
  double *dataptr=data.fortran_vec();
  //return octave_value(data);
  int nread=gd_getdata(file,myfield,0,0,nframes,spf,GD_FLOAT64,dataptr);
  //int nread=gd_getdata(file,myfield,0,0,nframes-1,spf,GD_FLOAT64,dataptr);
  //printf("Made it here.\n");
  myerr=gd_error(file);
  if (myerr!=GD_E_OK) {
    fprintf(stderr,"Error reading data from channel %s\n",myfield);
    return octave_value(myerr);
  }
  if (nread<spf*nframes) {
    fprintf(stderr,"Only read %d samples, expected %d on dirfile channel %s.  You may have some extra zeros.\n",nread,spf*nframes,myfield);
  }
  gd_raw_close(file,myfield);

  octave_value_list retval;
  retval(0)=data;
  retval(1)=nread;
  return retval;
  
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (getdata_double_channel_piece_c,args,nargout,"Read a piece of a dirfile channel as a double.  Args are (dirfile,channel name,first_frame,first_sample,num_frames,num_samples).\n")
{
  DIRFILE *file=( DIRFILE *)get_pointer(args(0));
  char *myfield=get_char_from_arg(args(1).char_matrix_value());
  long first_frame=(long)get_value(args(2));
  long first_sample=(long)get_value(args(3));
  long num_frames=(long)get_value(args(4));
  long num_samples=(long)get_value(args(5));

  long nframes=gd_nframes(file);
  //printf("nframes is %d\n",nframes);



  long spf=gd_spf(file,myfield);
  int myerr=gd_error(file);
  if (myerr!=GD_E_OK) {
    fprintf(stderr,"Error reading samples per frame on channel %s\n",myfield);
    return octave_value(myerr);
  }
  //printf("samples per frame is %d\n",spf);
  
  long ndata=gd_getdata(file,myfield,first_frame,first_sample,num_frames,num_samples,GD_NULL,NULL);

  Matrix data(ndata,1);
  double *dataptr=data.fortran_vec();
  //return octave_value(data);
  int nread=gd_getdata(file,myfield,first_frame,first_sample,num_frames,num_samples,GD_FLOAT64,dataptr);
  //int nread=gd_getdata(file,myfield,0,0,nframes-1,spf,GD_FLOAT64,dataptr);
  //printf("Made it here.\n");
  myerr=gd_error(file);
  if (myerr!=GD_E_OK) {
    fprintf(stderr,"Error reading data from channel %s\n",myfield);
    return octave_value(myerr);
  }
  //if (nread<spf*nframes) {
  //  fprintf(stderr,"Only read %d samples, expected %d on dirfile channel %s.  You may have some extra zeros.\n",nread,spf*nframes,myfield);
  //}
  gd_raw_close(file,myfield);

  octave_value_list retval;
  retval(0)=data;
  retval(1)=nread;
  return retval;
  
}
