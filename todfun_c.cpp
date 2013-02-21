#include <octave/oct.h>
#include <octave/Cell.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <ninkasi_config.h>
#include <ninkasi.h>
#include <mbCuts.h>
#include <ninkasi_mathutils.h>
#include <noise.h>
#include <dirfile.h>
#include <readtod.h>
#include <math.h>
#ifdef _MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif

#ifndef NO_FFTW
#include <fftw3.h>
#else
#include <fftw/fftw3.h>
#endif

#include <omp.h>

void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;



void dgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc)
{
  dgemm_(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}


/*--------------------------------------------------------------------------------*/
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

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
char *get_char_from_ov(octave_value cch)
{
  charMatrix ch=cch.char_matrix_value();

  int nn=ch.length();
  char *s=ch.fortran_vec();
  char *ss=strndup(s,nn+1);
  ss[nn]='\0';
  return ss;
}


/*--------------------------------------------------------------------------------*/
double **matrix2cmat(double *vec, dim_vector dm)

{
  //double *vec=mat.fortran_vec();
  //dim_vector dm=mat.dims();
  int ny=dm(0);
  int nx=dm(1);
  double **mm=(double **)malloc(sizeof(double *)*nx);
  for (int i=0;i<nx;i++)
    mm[i]=vec+i*ny;
  return mm;

}
/*--------------------------------------------------------------------------------*/
int matrix_nelem(Matrix mat)
{
  dim_vector dm=mat.dims();
  return dm(0)*dm(1);
}
/*--------------------------------------------------------------------------------*/

actData get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  actData myval=(actData)myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (test_omp_c, args, nargout, "Read TOD data into memory.\n")
{

  if (args.length()>0) {
    int n=(int)get_value(args(0));
    printf("Setting threads to %d\n",n);
    omp_set_num_threads(n);
  }


#pragma omp parallel
  printf("I am %d of %d\n",omp_get_thread_num(),omp_get_num_threads());

  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_dist_c, args, nargout, "Find the closest approach a TOD makes to an RA/Dec.\n")
{
  if (args.length()<3) {
    printf("Error - need at least three arguments to get_tod_dist_c.\n");
    return octave_value_list();
  }
  actData ra=get_value(args(0));
  actData dec=get_value(args(1));
  mbTOD  *mytod=(mbTOD *)get_pointer(args(2));

  return octave_value(how_far_am_i_from_radec_radians(ra,dec,mytod));
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_dist_vec_c, args, nargout, "Find the closest approach a TOD makes to an RA/Dec vector.\n")
{
  if (args.length()<3) {
    printf("Error - need at least three arguments to get_tod_dist_c.\n");
    return octave_value_list();
  }
  Matrix ra=args(0).matrix_value();
  Matrix dec=args(1).matrix_value();
  double *ra_vec=ra.fortran_vec();
  double *dec_vec=dec.fortran_vec();
  dim_vector dm=ra.dims();
  Matrix dists(dm);
  double *dists_vec=dists.fortran_vec();

  int npt=dm(1)*dm(0);
  mbTOD  *mytod=(mbTOD *)get_pointer(args(2));
  double *vec=how_far_am_i_from_radec_radians_vec(ra_vec,dec_vec,npt,mytod);
  for (int i=0;i<npt;i++)
    dists_vec[i]=vec[i];
  free(vec);


  return octave_value(dists);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (read_tod_data, args, nargout, "Read TOD data into memory.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data)
    allocate_tod_storage(mytod);
  if (mytod->data_saved) 
    memcpy(mytod->data[0],mytod->data_saved[0],mytod->ndet*mytod->ndata*sizeof(actData));
  else {
    //printf("reading tod data.\n");
    read_tod_data(mytod);
  }
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (reverse_tod_data, args, nargout, "Read TOD data into memory.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  reverse_tod_data(mytod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_ndata, args, nargout, "How many data elements in a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  return octave_value(mytod->ndata);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_dt, args, nargout, "How long is a sample in a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  return octave_value(mytod->deltat);
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_tod_tvec_c, args, nargout, "How long is a sample in a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->dt) {
    dim_vector dm(mytod->ndata,1);
    Matrix tvec(dm);
    for (int i=0;i<mytod->ndata;i++)
      tvec(i)=mytod->dt[i];
    return octave_value(tvec);
  }
  else
    return octave_value(0);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_ctime_c, args, nargout, "How long is a sample in a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  return octave_value(mytod->ctime);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_name_c, args, nargout, "What is a TOD's name?.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  return octave_value(mytod->dirfile);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (cut_detector_c, args, nargout, "Cut a detector.  args are (tod,row,col).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int row=(int)get_value(args(1));
  int col=(int)get_value(args(2));
  mbCutsSetAlwaysCut(mytod->cuts,row,col);
  
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (read_cuts_c, args, nargout, "What is a TOD's name?.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  char *myfroot=get_char_from_arg(args(1).char_matrix_value());  
  //mytod->cuts=mbCutsAlloc(mytod->nrow,mytod->ncol); //should already be allocated
  printf("reading cuts from %s\n",myfroot);
  mbCutsRead(mytod->cuts, myfroot);
  printf("cutting mispointed detectors.\n");
  cut_mispointed_detectors(mytod);
  printf("cut.\n");
  
  for (int i=0;i<mytod->decimate;i++)
    mbCutsDecimate(mytod->cuts);
  
  free(myfroot);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (write_cuts_c, args, nargout, "Write cuts to disk.  Args are tod,fname.\n")
{
  if (args.length() !=2) {
    printf("error - need exactly 2 inputs in write_cuts_c, not %d\n",args.length());
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  char *fname=get_char_from_arg(args(1).char_matrix_value());

  if (mytod->cuts)
    mbCutsWrite(mytod->cuts,fname);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (merge_cuts_c, args, nargout, "Merge cuts, write 'em to disk.\n")
{
  if (args.length() !=3) {
    printf("error - need exactly 3 inputs, not %d\n",args.length());
    return octave_value_list();
  }
  int n=args(0).length();
  assert(args(1).length()==n);
  assert(args(2).length()==n);


  Cell c1 = args (0).cell_value ();
  Cell c2 = args (1).cell_value ();
  Cell cc = args (2).cell_value ();



#pragma omp parallel for shared(c1,c2,cc,n) default(none) schedule(dynamic,2)
  for (int i=0;i<n;i++) {
    char *s1=get_char_from_ov(c1.elem (i));
    char *s2=get_char_from_ov(c2.elem (i));
    char *ss=get_char_from_ov(cc.elem (i));

    printf("going to merge %s and %s into %s\n",s1,s2,ss);
    mbCuts **cuts_in=(mbCuts **)malloc(sizeof(mbCuts *)*2);
    cuts_in[0]=mbCutsAlloc(33,32);
    cuts_in[1]=mbCutsAlloc(33,32);
    mbCutsRead(cuts_in[0], s1);
    mbCutsRead(cuts_in[1], s2);
    mbCuts *cuts=mbCutsAllocFromCuts((const mbCuts **)cuts_in,2);
    
    if (cuts)
      mbCutsWrite(cuts,ss);
    

    //huge memory leak.
    //CutsFree(cuts_in[0]);
    //CutsFree(cuts_in[1]);
    //CutsFree(cuts);
    

    free(s1);
    free(s2);
    free(ss);
  }
  
  return octave_value_list();

}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (print_tod_uncut_regions, args, nargout, "Print out the uncut regions of a detector in a tod.\n")
{
  
  if (args.length()<3) {
    fprintf(stderr,"Error - usage is (tod,row,column)\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int myrow=(int)get_value(args(1));
  int mycol=(int)get_value(args(2));
  int do_cuts=0;
  if (args.length()>3)
    do_cuts=(int)get_value(args(3));
  
  if (mbCutsIsAlwaysCut(mytod->cuts,myrow,mycol)) {
    printf("detector %d %d is completely cut.\n",myrow,mycol);
    return octave_value_list();
  }
  
#if 1
  mbUncut *uncuts=NULL;
  if (do_cuts==1) {
    if (mytod->cuts_as_uncuts==NULL) {
      fprintf(stderr,"Do not have cuts_as_uncuts filled in print_tod_uncut_regions.\n");
      return octave_value_list();
    }    
    uncuts=mytod->cuts_as_uncuts[myrow][mycol];
  }
  if (do_cuts==2) {
    if (mytod->cuts_as_vec==NULL) {
      fprintf(stderr,"Do not have cuts_as_uncuts filled in print_tod_uncut_regions.\n");
      return octave_value_list();
    }    
    uncuts=mytod->cuts_as_vec[myrow][mycol];
  }
  if (do_cuts==0)
    uncuts=mytod->uncuts[myrow][mycol];
  if (uncuts==NULL) {
    printf("do_cuts value not recognized in print_tod_uncut_regions, value is %d\n",do_cuts);
    return octave_value_list();
  }

  if (nargout==0)
    printf("Have %d regions.\n",uncuts->nregions);
  Matrix uncut(uncuts->nregions,2);
  for (int i=0;i<uncuts->nregions;i++) {
    if (nargout==0)
      printf("using %5d to %5d.\n",uncuts->indexFirst[i],uncuts->indexLast[i]);
    uncut(i,0)=uncuts->indexFirst[i];
    uncut(i,1)=uncuts->indexLast[i];
  }  
#else
  if (nargout==0)
    printf("Have %d regions.\n",mytod->uncuts[myrow][mycol]->nregions);
  Matrix uncut(mytod->uncuts[myrow][mycol]->nregions,2);
  for (int i=0;i<mytod->uncuts[myrow][mycol]->nregions;i++) {
    if (nargout==0)
      printf("using %5d to %5d.\n",mytod->uncuts[myrow][mycol]->indexFirst[i],mytod->uncuts[myrow][mycol]->indexLast[i]);
    uncut(i,0)=mytod->uncuts[myrow][mycol]->indexFirst[i];
    uncut(i,1)=mytod->uncuts[myrow][mycol]->indexLast[i];
  }
#endif  
  return octave_value(uncut);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_tod_window_c, args, nargout, "Cut the ends of a TOD.  Arguments are (tod, cut_length (in seconds)).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData tcut=get_value(args(1));
  set_tod_window(mytod,tcut);
  return octave_value(mytod->n_to_window);
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (cut_tod_ends_c, args, nargout, "Cut the ends of a TOD.  Arguments are (tod, cut_length (in seconds)).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData tcut=get_value(args(1));
  cut_tod_ends(mytod,tcut);
  //int nsamp=tcut/mytod->deltat;
  //mbCutsExtendGlobal(mytod->cuts,0,nsamp);
  //mbCutsExtendGlobal(mytod->cuts,mytod->ndata-nsamp,mytod->ndata);
  return octave_value(mytod->n_to_window);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (write_tod_data_c, args, nargout, "Write out the data in a TOD. \n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  char *myfroot=get_char_from_arg(args(1).char_matrix_value());
  
  FILE *outfile=fopen(myfroot,"w");
  if (!outfile) {
    fprintf(stderr,"Error writing to %s in write_data_c.\n",myfroot);
    return octave_value_list();
    
  }
  fwrite(&(mytod->ndet),1,sizeof(int),outfile);
  fwrite(&(mytod->ndata),1,sizeof(int),outfile);
  fwrite(mytod->rows,mytod->ndet,sizeof(int),outfile);
  fwrite(mytod->cols,mytod->ndet,sizeof(int),outfile);
  fwrite(mytod->data[0],mytod->ndet*mytod->ndata,sizeof(actData),outfile);
  fclose(outfile);
  free(myfroot);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (window_data_c, args, nargout, "Window the ends of a TOD. \n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int unwindow=0;
  if (args.length()>1)
    unwindow=(int)get_value(args(1));
  if (unwindow)
    unwindow_data(mytod);
  else
    window_data(mytod);
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (highpass_tod_c, args, nargout, "Highpass a TOD.  Args are (tod,nu_low,nu_high). \n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData t1=get_value(args(1));
  actData t2=get_value(args(2));
  highpass_tod(mytod,t1,t2);

  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (gapfill_data_c, args, nargout, "Window the ends of a TOD. \n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  //printf("tod->data and tod->ndet are %d %d\n",mytod->ndata,mytod->ndet);
  if (!mytod->have_data) {
    printf("Warning - tried to gapfill nonexistent data.\n");
    return octave_value_list();
  }
  
  fill_gaps_stupid(mytod);
  
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (cut_tod_global_c, args, nargout, "Do a global cut on a TOD.  Args are (tod, min_samp,max_samp).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->cuts==NULL) {
    fprintf(stderr,"Warning in cut_tod_global_c - TOD does not have cuts.\n");
    return octave_value_list();
  }
  int first=(int)get_value(args(1));
  int last=(int)get_value(args(2));
  if (first<0)
    first=0;
  if (last>=mytod->ndata)
    last=mytod->ndata-1;

  mbCutsExtendGlobal(mytod->cuts,first,last);

  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (setup_cut_fit_params, args, nargout, "Initialize storage for cuts parameterized with polynomials.\n")
{
  if (args.length()<2) {
    printf("need two arguments in setup_cut_fit_params - (tod,map of gap length to # of parameters)\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix vec=args(1).matrix_value();
  double *vv=vec.fortran_vec();
  dim_vector dm=vec.dims();
  int nelem=dm(0)*dm(1);
  int *ivec=(int *)malloc(sizeof(int)*nelem);
  for (int i=0;i<nelem;i++)
    ivec[i]=vv[i];
  mytod->cuts_fit_params=setup_cut_fit_params(mytod,ivec);
  free(ivec);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (are_cutvecs_fitparams, args, nargout, "See if cutvecs are fit params.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->cuts_fit_params)
    return octave_value(true);
  else
    return octave_value(false);
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (initialize_cutvecs_fitparams_precon, args, nargout, "Setup fitparams preconditionerss.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  setup_cutsfits_precon(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (apply_cutvecs_fitparams_precon, args, nargout, "Setup fitparams preconditionerss.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix vec_in=args(1).matrix_value();
  dim_vector dm=vec_in.dims();
  //printf("have %d elements in cutvecs params.\n",dm(0)*dm(1));
  //printf("tod has %d detectors.\n",mytod->ndet);
  Matrix vec_out(dm);
  actData *vv_in=vec_in.fortran_vec();
  actData *vv_out=vec_out.fortran_vec();
  memset(vv_out,0,sizeof(actData)*dm(0)*dm(1));
  apply_cutfits_precon(mytod,vv_in,vv_out);

  return octave_value(vec_out);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_numel_cut_c, args, nargout, "Find out how many elements have been cut from a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int ncut=get_numel_cut(mytod);
  return octave_value(ncut);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_cuts_statistics, args, nargout,"Find statistics of cuts lengths in a TOD.\n")
{
  mbTOD *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->cuts_as_uncuts==NULL) {
    fprintf(stderr,"cuts_as_uncut not set on in get_cuts_statistics.\n");
    return octave_value_list();    
  }
  Matrix vec(mytod->ndata,1);
  actData *vv=vec.fortran_vec();
  memset(vv,0,mytod->ndata*sizeof(actData));
  for (int det=0;det<mytod->ndet;det++) {
    if (mytod->cuts_fit_params) {
      mbCutFitParams *cut=mytod->cuts_fit_params[mytod->rows[det]][mytod->cols[det]];
      for (int i=0;i<cut->nregions;i++) {
	vv[cut->nparams[i]-1]++;
      }
    }
    else {
      mbUncut *cut=mytod->cuts_as_uncuts[mytod->rows[det]][mytod->cols[det]];
      for (int i=0;i<cut->nregions;i++) {
	int nelem=cut->indexLast[i]-cut->indexFirst[i];
	assert(nelem>0);
	vv[nelem-1]++;
      }
    }
  }
  return octave_value(vec);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (tod2cutvec_c, args, nargout, "Put cut data from a TOD into a vector.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int ncut=get_numel_cut(mytod);
  Matrix vec(ncut,1);
  memset(vec.fortran_vec(),0,sizeof(actData)*ncut);

  if (tod2cutvec(mytod,vec.fortran_vec())) {
    printf("had a problem in tod2cutvec.\n");
    return octave_value_list();
  }
  return octave_value(vec);

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (cutvec2tod_c, args, nargout, "Put cut data into a TOD.\n")
{
  if (args.length()<2) {
    printf("error - need at least two inputs to cutvec2tod_c.\n");
    return octave_value_list();
  }
  
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix vec=args(1).matrix_value();  

  dim_vector dm=vec.dims();
  int ncut=get_numel_cut(mytod);
  if (dm(1)*dm(0)<ncut) {
    fprintf(stderr,"Error - not enough elements in cutvec2tod_c.  Expected at least %d, got %d on tod %s\n",ncut, dm(1)*dm(2),mytod->dirfile);
    return octave_value_list();
  }
  
  if (cutvec2tod(mytod,vec.fortran_vec())) {
    printf("had a problem in tod2cutvec.\n");
  }
  return octave_value_list();
  
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_tod_cuts_global_indices_c, args, nargout, "Set up the global indexing to the cut regions in a tod.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  set_global_indexed_cuts(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_tod_cut_regions_c, args, nargout, "Set up the cut regions in a tod like the uncuts.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  get_tod_cut_regions(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_tod_uncut_regions_c, args, nargout, "Set up the uncut regions in a tod.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int do_cached=0;
  if (args.length()>1)
    do_cached=(int)get_value(args(1));
  mbUncut ***tmp=mytod->uncuts;
  
  get_tod_uncut_regions(mytod);
  if (do_cached) {
    mytod->uncuts_for_interp=mytod->uncuts;
    mytod->uncuts=tmp;
  }
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (initialize_tod_cut_params_c, args, nargout, "Set up the legendre fit parameters for gap filling.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix mat=args(1).matrix_value();
  actData *vec=mat.fortran_vec();
  dim_vector dm=mat.dims();
  if (dm(0)*dm(1)<mytod->ndata) {
    printf("Binning vector is not long enough in initialize_tod_cut_params_c.  Make sure it's the length of a timestream.\n");
    return octave_value_list();
  }
  int *ivec=(int *)malloc(sizeof(int)*mytod->ndata);
  for (int i=0;i<mytod->ndata;i++)
    ivec[i]=vec[i];
  mytod->cuts_fit_params=setup_cut_fit_params(mytod,ivec);
  free(ivec);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_tod_kept_regions_c, args, nargout, "Set up the kept regions in a tod.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  //get_tod_uncut_regions(mytod);
  mytod->kept_data=get_uncut_regions(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (reverse_tod_uncut_regions_c, args, nargout, "Set up the uncut regions in a tod.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  reverse_tod_uncut_regions(mytod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (detrend_data_c_c, args, nargout, "Set up the uncut regions in a tod.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  detrend_data(mytod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod2map, args, nargout, "Project tod into a map.  Args are (tod,map)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  MAP *mymap=(MAP *)get_pointer(args(1));

  //omp_set_num_threads(8); 
  
  tod2map(mymap,mytod,NULL);

  return octave_value_list();  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD(tod2polmap,args,nargout," Project a tod into a polmap.  Args are (tod,map)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  MAP *mymap=(MAP *)get_pointer(args(1));
  
  tod2polmap(mymap,mytod);

  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod_times_map, args, nargout, "Multiply and sum a tod by a map.  Args are (tod,map)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  MAP *mymap=(MAP *)get_pointer(args(1));

  //omp_set_num_threads(8); 
  
  return octave_value(tod_times_map(mymap,mytod,NULL));

}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (map2tod, args, nargout, "Project map into a tod.  Args are (map,tod).  Optionally, scale factor\n")
{
  MAP *mymap=(MAP *)get_pointer(args(0));
  mbTOD  *mytod=(mbTOD *)get_pointer(args(1));
  if (args.length()>2) {
    actData scale_fac=get_value(args(2));
    
    map2tod(mymap,mytod,(PARAMS *)&scale_fac);
  }
  else
    map2tod(mymap,mytod,NULL);
  return octave_value_list();  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (polmap2tod, args, nargout, "Project polarizedmap into a tod.  Args are (map,tod).\n")
{
  MAP *mymap=(MAP *)get_pointer(args(0));
  mbTOD  *mytod=(mbTOD *)get_pointer(args(1));
  polmap2tod(mymap,mytod);
  return octave_value_list();  
}

/*--------------------------------------------------------------------------------*/
int how_many_dets_are_kept_c(const mbTOD *tod)
{
  int i;
  int ngood=0;
  for (i=0;i<tod->ndet;i++) {
    if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i))) {
      ngood++;
    }
  }
  return ngood;
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (how_many_dets_are_kept, args, nargout, "Find out how many detectors are kept.  Arg is (tod)\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int ngood=how_many_dets_are_kept_c(tod);
  return octave_value(ngood);  
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (is_det_always_cut, args, nargout, "Find out how many detectors are kept.  Arg is (tod)\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int row=(int)get_value(args(1));
  int col=(int)get_value(args(2));

  return octave_value(mbCutsIsAlwaysCut(tod->cuts,row,col));
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (purge_cut_detectors, args, nargout, "Get rid of detectors that are cut.  Arg is (tod)\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  purge_cut_detectors(tod);
  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/
/*
DEFUN_DLD (how_many_dets_are_kept, args, nargout, "Find out how many detectors are kept.  Arg is (tod)\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int i;
  int ngood=0;
  for (i=0;i<tod->ndet;i++) {
    if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i))) {
      ngood++;
    }
  }
  printf("keeping %d out of %d\n",ngood,tod->ndet);
  return octave_value_list();  
}

*/
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (sum_tod_data, args, nargout, "Sum up total abs(data).  Arg is (tod)\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (!tod->have_data) {
    printf("Do not have data.\n");
    return octave_value(0);
  }
  int i,j;

  double tot=0;
  for (i=0;i<tod->ndet;i++) {
    if ((!mbCutsIsAlwaysCut(tod->cuts,tod->rows[i],tod->cols[i]))&&(is_det_listed(tod,NULL,i))) {
      for (j=0;j<tod->ndata;j++)
	tot+=fabs(tod->data[i][j]);
    }
  }
  return octave_value(tot);
  //return octave_value_list();  
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (allocate_tod_storage_c, args, nargout, "Allocate storage for a tod.  Arg is (tod)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  allocate_tod_storage(mytod);
  return octave_value_list();  
}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (free_tod_storage_c, args, nargout, "Free storage for a tod.  Arg is (tod)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  if (mytod->have_data)
    free_tod_storage(mytod);
  else
    printf("Warning - attempted to free already empty TOD.  Skipping.\n");
  return octave_value_list();  
}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (assign_tod_value_c, args, nargout, "Assign a value to a tod.  Args are (tod,value)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  NDArray vv=args(1).array_value();
  actData val=vv(0,0);
  assign_tod_value(mytod,val);
  return octave_value_list();  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (tod_hits_source_c, args, nargout, "Check to see if a tod might hit a source Args are (ra,dec,dist,tod)\n")
{
  actData ra=get_value(args(0));
  actData dec=get_value(args(1));
  actData r=get_value(args(2));
  mbTOD  *mytod=(mbTOD *)get_pointer(args(3));
  
  int val=tod_hits_source(ra,dec,r,mytod);
  return octave_value(val);

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_src2tod, args, nargout, "Add a source into a tod.  Args are (tod,ra,dec,amp,beam vector, dtheta)\n")
//void add_src2tod(mbTOD *tod, actData ra, actData dec, actData src_amp, const actData *beam, actData dtheta, int nbeam, int oversamp)

{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  actData ra=get_value(args(1));
  actData dec=get_value(args(2));
  actData amp=get_value(args(3));
  Matrix beam=args(4).matrix_value();
  actData dtheta=get_value(args(5));
  dim_vector dm=beam.dims();
  int nbeam=dm(0)*dm(1);
  actData *beamvec=beam.fortran_vec();
  int oversamp=1;
  if (args.length()>6)
    oversamp=(int)get_value(args(6));
  add_src2tod(mytod,ra,dec,amp,beamvec,dtheta,nbeam,oversamp);
  return octave_value_list();  
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_srcvec2tod, args, nargout, "Add a source into a tod.  Args are (tod,ra,dec,amp,beam vector, dtheta)\n")
//void add_src2tod(mbTOD *tod, actData ra, actData dec, actData src_amp, const actData *beam, actData dtheta, int nbeam, int oversamp)

{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  Matrix ram=args(1).matrix_value();
  Matrix decm=args(2).matrix_value();
  Matrix ampm=args(3).matrix_value();
  Matrix beam=args(4).matrix_value();
  actData *ra=ram.fortran_vec();
  actData *dec=decm.fortran_vec();
  actData *amp=ampm.fortran_vec();
  dim_vector srcdm=ram.dims();
  int nsrc=srcdm(0)*srcdm(1);

  actData dtheta=get_value(args(5));
  dim_vector dm=beam.dims();
  int nbeam=dm(0)*dm(1);
  actData *beamvec=beam.fortran_vec();
  int oversamp=1;
  if (args.length()>6)
    oversamp=(int)get_value(args(6));
  add_srcvec2tod(mytod,ra,dec,amp,nsrc,beamvec,dtheta,nbeam,oversamp);
  return octave_value_list();  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod2srcvec, args, nargout, "Project a tod into source vecs.  Args are (tod,ra,dec,beam vector, dtheta,[oversamp])\n")
//void add_src2tod(mbTOD *tod, actData ra, actData dec, actData src_amp, const actData *beam, actData dtheta, int nbeam, int oversamp)

{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  Matrix ram=args(1).matrix_value();
  Matrix decm=args(2).matrix_value();
  Matrix beam=args(3).matrix_value();
  actData *ra=ram.fortran_vec();
  actData *dec=decm.fortran_vec();
  dim_vector srcdm=ram.dims();
  int nsrc=srcdm(0)*srcdm(1);

  actData dtheta=get_value(args(4));
  dim_vector dm=beam.dims();
  int nbeam=dm(0)*dm(1);
  actData *beamvec=beam.fortran_vec();
  int oversamp=1;
  if (args.length()>5)
    oversamp=(int)get_value(args(5));
  
  Matrix amps(srcdm);
  memset(amps.fortran_vec(),0,sizeof(double)*nsrc);
  tod2srcvec(amps.fortran_vec(),mytod,ra,dec,nsrc,beamvec,dtheta,nbeam,oversamp);
  return octave_value(amps);  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_vector_to_tod_data, args, nargout, "Adds a vector to each detector in tod data.  Optionally, scale vector by a third agrument.   Args are (tod,vector,[ampt])\n")
{
  if (args.length()<2) {
    fprintf(stderr,"Need at least 2 args to add_vector_to_tod_data.\n");
    return octave_value_list();
  }
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (!tod->have_data) {
    fprintf(stderr,"Error in add_matrix_to_tod_data - data not allocated in tod.\n");
    return octave_value_list();    
  }
  Matrix vec=args(1).matrix_value();
  double *vv=vec.fortran_vec();
  if (args.length()==2) {
    for (int i=0;i<tod->ndet;i++)
      for (int j=0;j<tod->ndata;j++)
	tod->data[i][j]+=vv[j];
  }
  if (args.length()==3) {
    Matrix scale_facs=args(3).matrix_value();
    double *ss=scale_facs.fortran_vec();
    for (int i=0;i<tod->ndet;i++)
      for (int j=0;j<tod->ndata;j++)
	tod->data[i][j]+=vv[j]*ss[i];
  }
  
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_matrix_to_tod_data, args, nargout, "Adds a matrix to tod data.  Args are (tod,matrix,[value])\n")

{

  if (args.length()<2) {
    fprintf(stderr,"Need at least 2 args to add_matrix_to_tod_data.\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  if (!mytod->have_data) {
    fprintf(stderr,"Error in add_matrix_to_tod_data - data not allocated in tod.\n");
    return octave_value_list();

  }
  Matrix mat=args(1).matrix_value();
  dim_vector dm=mat.dims();
  
  double fac=1.0;
  if (args.length()>2)
    fac=get_value(args(2));
  if ((dm(1)!=mytod->ndet)||(dm(0)!=mytod->ndata)) {
    fprintf(stderr,"Error in add_matrix_to_tod_data - size mismatch.  Got %d %d, expected %d %d new\n",dm(0),dm(1),mytod->ndata,mytod->ndet);
    return octave_value_list();
  }

  double *vec=mat.fortran_vec();
  for (int i=0;i<mytod->ndet;i++)
    for (int j=0;j<mytod->ndata;j++) {
      mytod->data[i][j]+=fac*vec[i*mytod->ndata+j];
    }
  return octave_value_list();  
}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (read_dirfile_tod_data_from_rowcol_list_c, args, nargout, "Get the data from a list of rows & columns into octave.   Does not check cuts, so can be used for dark detectors.  Args are (tod,row,col)\n")
{
  if (args.length()!=3) {
    fprintf(stderr,"Error - expected 3 arguments in read_dirfile_tod_data_from_rowcol_list_c, got %d.\n",args.length());
    return octave_value_list();
  }
    
  assert(args.length()==3);
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix in_rows=args(1).matrix_value();
  Matrix in_cols=args(2).matrix_value();
  dim_vector dm_rows=in_rows.dims();
  dim_vector dm_cols=in_cols.dims();
  int ndet=dm_rows(0)*dm_rows(1);
  assert(dm_cols(0)*dm_cols(1)==ndet);

  int *rows=(int *)malloc(sizeof(int)*ndet);
  int *cols=(int *)malloc(sizeof(int)*ndet);
  double *rvec=in_rows.fortran_vec();
  double *cvec=in_cols.fortran_vec();
  for (int i=0;i<ndet;i++) {
    rows[i]=(int)rvec[i];
    cols[i]=(int)cvec[i];
  }
  dim_vector dm(tod->ndata,ndet);
  Matrix data(dm);
#ifdef ACTDATA_DOUBLE
  actData **dd=(actData **)malloc(sizeof(actData *)*tod->ndet);
  actData *vec=data.fortran_vec();
  for (int i=0;i<ndet;i++)
    dd[i]=vec+i*tod->ndata;
  read_dirfile_tod_data_from_rowcol_list(tod,rows,cols,ndet,dd);
#else
  fprintf(stderr,"Single precision not yet implemented.\n");
#endif

  free(dd);
  free(rows);
  free(cols);
  return octave_value(data);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_altaz, args, nargout, "Return the alt/az of the boresight for a tod.  Arg is (tod)\n")
{
  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  Matrix alt(mytod->ndata,1);
  Matrix az(mytod->ndata,1);
  double *altp=alt.fortran_vec();
  double *azp=az.fortran_vec();
  for (int i=0;i<mytod->ndata;i++) {
    altp[i]=mytod->alt[i];
    azp[i]=mytod->az[i];
  }
  octave_value_list retval;
  retval(0)=octave_value(alt);
  retval(1)=octave_value(az);
  return retval;
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_hwp, args, nargout, "Return the hwp angle of a tod.  Arg is (tod)\n")
{
  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  if (!mytod->hwp) {
    fprintf(stderr,"Error in get_tod_hwp - HWP not found in TOD.\n");
    return octave_value_list();
  }
  Matrix hwp(mytod->ndata,1);
  double *hwpp=hwp.fortran_vec();
  for (int i=0;i<mytod->ndata;i++) {
    hwpp[i]=mytod->hwp[i];
  }
  return octave_value(hwp);
}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_data, args, nargout, "Get the data from a tod into octave.  Only return non-cut dets.  Arg is (tod)\n")
{

  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  

  if (!mytod->have_data) {
    printf("Tod does not have data loaded.\n");
    return octave_value_list();
  }
  int ndet= how_many_dets_are_kept_c(mytod);
  Matrix dat(mytod->ndata,ndet);
  double *dd=dat.fortran_vec();
  int icur=0;
  int j;
  for (int i=0;i<mytod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
      for (j=0;j<mytod->ndata;j++)
	dd[icur*mytod->ndata+j]=mytod->data[i][j];
      icur++;
    }
  }
  
  return octave_value(dat);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_2gamma, args, nargout, "Get 2gamma from a tod into octave.  Only return non-cut dets.  Arg is (tod)\n")
{

  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));

  if (!mytod->twogamma_saved) {
    printf("Tod does not have twogamma loaded.\n");
    return octave_value_list();
  }
  int ndet= how_many_dets_are_kept_c(mytod);
  Matrix dat(mytod->ndata,ndet);
  double *dd=dat.fortran_vec();
  int icur=0;
  int j;
  for (int i=0;i<mytod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
      for (j=0;j<mytod->ndata;j++)
        dd[icur*mytod->ndata+j]=mytod->twogamma_saved[i][j];
      icur++;
    }
  }
  
  return octave_value(dat);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_pixellization, args, nargout, "Get the pixellization from a tod into octave.  Only return non-cut dets.  Arg is (tod)\n")
{

  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));

  if (!mytod->pixelization_saved) {
    printf("Tod does not have pixellization loaded.\n");
    return octave_value_list();
  }
  int ndet= how_many_dets_are_kept_c(mytod);
  dim_vector dm(mytod->ndata,ndet);
  int32NDArray dat(dm);
  int *dd=(int *)dat.fortran_vec();
  int icur=0;
  int j;
  for (int i=0;i<mytod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
      for (j=0;j<mytod->ndata;j++)
        dd[icur*mytod->ndata+j]=mytod->pixelization_saved[i][j];
      icur++;
    }
  }
  
  return octave_value(dat);
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (apply_calib_facs_c, args, nargout, "Apply calibration factors to a TOD.  args are (tod,facs).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data) {
    fprintf(stderr,"Data must be present in apply_calib_facs_c.\n");
    return octave_value_list();
  }
  if (mytod->calib_facs_saved) {
    printf("Using saved calibration factors.\n");
    for (int i=0;i<mytod->ndet;i++) {
      actData fac=mytod->calib_facs_saved[mytod->rows[i]][mytod->cols[i]];
      for (int j=0;j<mytod->ndata;j++)
	mytod->data[i][j]*=fac;
    }						
    return octave_value_list();
  }

  Matrix facs=args(1).matrix_value();
  dim_vector dm=facs.dims();
  if ((dm(0)!=mytod->ndet)||(dm(1)!=1)) {
    fprintf(stderr,"Error in apply_calib_facs_c - unexpected input size.\n");
    return octave_value_list();
  }

  //printf("about to calibate on %d detectors.\n",mytod->ndet);
  
#pragma omp parallel for shared(mytod,facs) default(none)
  for (int i=0;i<mytod->ndet;i++) {
    actData fac=facs(i,0);
    for (int j=0;j<mytod->ndata;j++)
      mytod->data[i][j]*=fac;
  }
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod_has_calib_facs_c, args, nargout, "Does a TOD have calibration factors saved in it?.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->calib_facs_saved)
    return octave_value(1);
  else
    return octave_value(0);
  return octave_value_list(); //never get here.
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (save_calib_facs_c, args, nargout, "Save calibration factors into a TOD.  args are (tod,facs).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix facs=args(1).matrix_value();
  dim_vector dm=facs.dims();
  //printf("dimensions are %d %d, vs. %d\n",dm(0),dm(1),mytod->ndet);
  if ((dm(0)!=mytod->ndet)||(dm(1)!=1)) {
    fprintf(stderr,"Error in save_calib_facs_c - unexpected input size.\n");
    return octave_value_list();
  }
  

  if (mytod->calib_facs_saved==NULL) 
    mytod->calib_facs_saved=matrix(mytod->nrow,mytod->ncol);

  for (int i=0;i<mytod->ndet;i++) 
    mytod->calib_facs_saved[mytod->rows[i]][mytod->cols[i]]=facs(i,0);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_tod_row_c, args, nargout, "Get the used rows from a tod into octave.  Only return non-cut dets.  Arg is (tod)\n")
{
  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  int ndet=how_many_dets_are_kept_c(mytod);
  Matrix rows(ndet,1);
  double *rr=rows.fortran_vec();
  int icur=0;
  for (int i=0;i<mytod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
      assert(icur<ndet);
      rr[icur]=mytod->rows[i];
      icur++;
    }
  }
  
  if (icur != ndet) {
    fprintf(stderr,"Warning - cuts screwiness in get_tod_row_c.\n");
  }

  return octave_value(rows);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_col_c, args, nargout, "Get the used cols from a tod into octave.  Only return non-cut dets.  Arg is (tod)\n")
{
  const mbTOD  *mytod=(mbTOD *)get_pointer(args(0));  
  int ndet=how_many_dets_are_kept_c(mytod);
  Matrix cols(ndet,1);
  double *cc=cols.fortran_vec();
  int icur=0;
  for (int i=0;i<mytod->ndet;i++) {
    if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
      assert(icur<ndet);
      cc[icur]=mytod->cols[i];
      icur++;
    }
  }
  
  if (icur != ndet) {
    fprintf(stderr,"Warning - cuts screwiness in get_tod_col_c.\n");
  }

  return octave_value(cols);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (assign_bad_timestreams, args, nargout, "Push a matrix of bad timestreams into a TOD.   Args are (data,tod)\n")
{
  Matrix dat=args(0).matrix_value();
  mbTOD *mytod=(mbTOD *)get_pointer(args(1));
  
  dim_vector dm=dat.dims();
  if (dm(0)!= mytod->ndata) {
    fprintf(stderr,"Problem in assign_bad_timestreams - mismatch in length, expected %d, got %d\n",mytod->ndata,dm(0));
    return octave_value_list();
  }

  
  BadTimestreams *timestreams=(BadTimestreams *)malloc(sizeof(BadTimestreams));
  mytod->bad_timestreams=timestreams;

  mytod->bad_timestreams->nstream=dm(1);
  mytod->bad_timestreams->n=mytod->ndata;
  mytod->bad_timestreams->timestreams=matrix(dm(1),dm(0));
  //printf("nstream is %d\n",mytod->bad_timestreams->nstream);
  double *vec=dat.fortran_vec();
  memcpy(mytod->bad_timestreams->timestreams[0],vec,dm(1)*dm(0)*sizeof(double));
  
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (pull_bad_timestreams, args, nargout, "Pull the bad timestreams out of a TOD, return the vectors.  Arg is (tod)\n")
{
  mbTOD *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->bad_timestreams==NULL) {
    fprintf(stderr,"Problem in pull_bad_timestreams - TOD is missing bad timestreams!\n");
    return octave_value_list();    
  }
  
  dim_vector dm(2);
  dm(0)=mytod->bad_timestreams->n;
  dm(1)=mytod->bad_timestreams->nstream;
  Matrix vecs(dm);
  double *vec=vecs.fortran_vec();
  memcpy(vec,mytod->bad_timestreams->timestreams[0],dm(0)*dm(1)*sizeof(double));
  
  return octave_value(vecs);
  
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (push_tod_altaz, args, nargout, "Push octave alt and az into a TOD.   Args are (data,alt,az)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix alt_mat=args(1).matrix_value();
  Matrix az_mat=args(2).matrix_value();
  double *alt=alt_mat.fortran_vec();
  double *az=az_mat.fortran_vec();

  for (int i=0;i<mytod->ndata;i++) {
    mytod->alt[i]=alt[i];
    mytod->az[i]=az[i];
  }
  
  
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (push_tod_data, args, nargout, "Push an octave matrix into a TOD.   Args are (data,tod)\n")
{
  
  NDArray dat=args(0).array_value();
  mbTOD  *mytod=(mbTOD *)get_pointer(args(1));  
  dim_vector dm=dat.dims();
  if (mytod->ndata!=dm.elem(0)) {
    fprintf(stderr,"Error in push_tod_data.  Size mismatch %d %d\n",dm.elem(0),mytod->ndata);
    return octave_value_list();
  }
  
  
  int ndet= how_many_dets_are_kept_c(mytod);
  if (ndet!=dm.elem(1)) {
    fprintf(stderr,"Error in push_tod_data.  ndet mismatch %d %d\n",dm.elem(1),ndet);
    return octave_value_list();        
  }
  
  if (!mytod->have_data) {
    allocate_tod_storage(mytod);
  }
  
  double *dd=dat.fortran_vec();
  
  int icur=0;
  for (int i=0;i<ndet;i++) {
    while ((mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[icur],mytod->cols[icur]))&&(icur<mytod->ndet))
      icur++;
    for (int j=0;j<mytod->ndata;j++)
      mytod->data[icur][j]=dd[i*mytod->ndata+j];
    icur++;
  }
  
  return octave_value(dat);

}


/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod_mean_c, args, nargout, "Multiply a row & column matrix, and add into a TOD.   Args are (tod,col_vecs,row_vecs)\n")
{  
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data) {
    fprintf(stderr,"Error - no data present in TOD when mean requested!\n");
    return octave_value_list();
  }
  actData *tot=(actData *)malloc(sizeof(actData)*mytod->ndata);
  memset(tot,0,sizeof(actData)*mytod->ndata);
#pragma omp parallel  shared(tot,mytod) default(none)
  {
    actData *mytot=(actData *)malloc(sizeof(actData)*mytod->ndata);
    memset(mytot,0,sizeof(actData)*mytod->ndata);
#pragma omp for
    for (int i=0;i<mytod->ndet;i++) {
      if (!mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[i],mytod->cols[i])) {
	for (int j=0;j<mytod->ndata;j++)
	  mytot[j]+=mytod->data[i][j];
      }
    }
#pragma omp critical
    for (int i=0;i<mytod->ndata;i++)
      tot[i]+=mytot[i];
    free(mytot);
  }
  actData ngood=(actData)how_many_dets_are_kept_c(mytod);
  ColumnVector tt(mytod->ndata);
  double *ptr=tt.fortran_vec();
  for (int i=0;i<mytod->ndata;i++)
    ptr[i]=tot[i]/ngood;
  free(tot); 
  return octave_value(tt);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (vecs2tod_blas, args, nargout, "Multiply a row & column matrix, and add into a TOD using blas.   Args are (tod,col_vecs,row_vecs)\n")
{

  mbTOD  *tod=(mbTOD *)get_pointer(args(0));

  Matrix cols=args(1).matrix_value();
  dim_vector dm_cols=cols.dims();
  Matrix rows=args(2).matrix_value();
  dim_vector dm_rows=rows.dims();
  double *colsptr=cols.fortran_vec();
  double *rowsptr=rows.fortran_vec();

  int nvec=dm_cols(1);
  assert(dm_cols(0)==tod->ndata);
  assert(dm_rows(1)==tod->ndet);
  assert(dm_cols(1)==dm_rows(0));

  //printf("dims on cols are %d %d, and on rows are %d %d\n",dm_cols(0),dm_cols(1),dm_rows(0),dm_rows(1));
  dgemm('N','N',tod->ndata,tod->ndet,nvec,1.0,colsptr,tod->ndata,rowsptr,nvec,1.0,tod->data[0],tod->ndata);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (timestreams2tod_blas, args, nargout, "Multiply a row & column matrix, and add into a TOD using blas.  Timestreams should be stored in TOD. Args are (tod,rows)\n")
{

  mbTOD  *tod=(mbTOD *)get_pointer(args(0));

  Matrix rows=args(1).matrix_value();
  dim_vector dm_rows=rows.dims();
  double *rowsptr=rows.fortran_vec();

  assert(dm_rows(0)==tod->bad_timestreams->nstream);
  assert(dm_rows(1)==tod->ndet);
  int nvec=tod->bad_timestreams->nstream;

  //printf("dims on cols are %d %d, and on rows are %d %d\n",dm_cols(0),dm_cols(1),dm_rows(0),dm_rows(1));
  dgemm('N','N',tod->ndata,tod->ndet,nvec,1.0,tod->bad_timestreams->timestreams[0],tod->ndata,rowsptr,nvec,1.0,tod->data[0],tod->ndata);
  return octave_value_list();
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (vecs2tod, args, nargout, "Multiply a row & column matrix, and add into a TOD.   Args are (tod,col_vecs,row_vecs)\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));

  Matrix cols=args(1).matrix_value();
  dim_vector dm_cols=cols.dims();
  Matrix rows=args(2).matrix_value();
  dim_vector dm_rows=rows.dims();

  if (dm_cols.elem(1)!=dm_rows.elem(0)) {
    printf("Inner dimension mismatch in vecs2tod %d %d\n",dm_cols.elem(1),dm_rows.elem(0));
    return octave_value_list();
  }
  
  if (dm_cols.elem(0)!=mytod->ndata) {
    printf("Length mismatch in vecs2tod %d %d\n",dm_cols.elem(0),mytod->ndata);
    return octave_value_list();
  }
  
  int ngood=how_many_dets_are_kept_c(mytod);
  if (dm_rows.elem(1)!=  ngood) {
    printf("Good detector mismatch in vecs2tod %d %d\n",dm_rows.elem(1),ngood);
    return octave_value_list();
  }


  
  if (!mytod->have_data) {
    allocate_tod_storage(mytod);
  }
  
  
  int nvec=dm_cols.elem(1);
  actData **mycols=matrix(mytod->ndata,nvec);
  actData **myrows=matrix(nvec,mytod->ndet);
  memset(mycols[0],0,mytod->ndata*nvec*sizeof(actData));
  memset(myrows[0],0,mytod->ndet*nvec*sizeof(actData));

  int icur=0;
  double *dd=rows.fortran_vec();
  for (int i=0;i<ngood;i++) {    
    while ((mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[icur],mytod->cols[icur]))&&(icur<mytod->ndet))
      icur++;
    for (int j=0;j<nvec;j++) {
      myrows[j][icur]=dd[j+i*nvec];
    }
    icur++;
  }
  
  double *cc=cols.fortran_vec();
  for (int i=0;i<mytod->ndata;i++) 
    for (int j=0;j<nvec;j++) {
      mycols[i][j]=cc[i+j*mytod->ndata];
    }
  
  
  int nndet=mytod->ndet;
  int nndata=mytod->ndata;

#if 1   //this should eventually be a blas call, but was seg faulting...
        //Maybe not - blas overhead for very small nvec appears to be large
  if (nvec>1) {
    printf("doing multiple.\n");
#pragma omp parallel for shared(mytod,mycols,myrows,nvec,nndata,nndet) default(none)
    for (int i=0;i<nndet;i++) {
      for (int j=0;j<nndata;j++) 
	for (int k=0;k<nvec;k++) {
	  mytod->data[i][j]+=mycols[j][k]*myrows[k][i];
	}
    }
  }
  else
    {

#if 1
      add_outer_product(myrows[0],nndet,mycols[0],nndata,mytod->data);
#else
      printf("doing single.\n");

#pragma omp parallel for shared(mytod,mycols,myrows,nvec,nndata,nndet) default(none)
      for (int i=0;i<nndet;i++)	  	
	for (int j=0;j<nndata;j++) 
	  mytod->data[i][j]+=mycols[0][j]*myrows[0][i];
#endif
    }

#else
  cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, mytod->ndata,mytod->ndet,nvec,1.0,mycols[0],nvec,myrows[0],mytod->ndet,1.0,mytod->data[0],mytod->ndata);
#endif
  free(mycols[0]);
  free(mycols);
  free(myrows[0]);
  free(myrows);
  
  return octave_value_list();

}



/*--------------------------------------------------------------------------------*/

DEFUN_DLD (clear_cut_data, args, nargout, "Set values of cut data to be equal to zero.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  clear_cut_data(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod2vecs_blas, args, nargout, "Multiply a set of rows by a TOD, and put it into a column matrix.\nAssumes cut regions have been zeroed.   Args are (tod,row_vecs)\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  if (!mytod->have_data) {
    printf("no data in tod!\n");
    return octave_value_list();
  }
  
  
  Matrix rows=args(1).matrix_value();
  dim_vector dm_rows=rows.dims();
  int nvec=dm_rows.elem(0);
  assert(dm_rows.elem(1)==mytod->ndet);  //no cut detectors!
  dim_vector dm_out(mytod->ndata,nvec);
  Matrix vecs(dm_out);
  double *vecs_ptr=vecs.fortran_vec();
  double *rows_ptr=rows.fortran_vec();
  dgemm('N','T',mytod->ndata,nvec,mytod->ndet,1.0,mytod->data[0],mytod->ndata,rows_ptr,nvec,0.0,vecs_ptr,mytod->ndata);

  return octave_value(vecs);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod2timestreams_blas, args, nargout, "Multiply a set of timestreams by a TOD, and put it into a row matrix.\nAssumes cut regions have been zeroed.   Arg is (tod)\n")
{
  
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  if (!mytod->have_data) {
    printf("no data in tod!\n");
    return octave_value_list();
  }
  
  dim_vector dm(2);
  int nvec=mytod->bad_timestreams->nstream;
  dm(0)=nvec;
  dm(1)=mytod->ndet;
  Matrix rows(dm);
  double *rowptr=rows.fortran_vec();
  //printf("calling blas.\n");
  dgemm('T','N',nvec,mytod->ndet,mytod->ndata,1.0,mytod->bad_timestreams->timestreams[0],mytod->ndata,mytod->data[0],mytod->ndata,0.0,rowptr,nvec);
  //printf("back from blas.\n");
  
#if 0 
  Matrix rows=args(1).matrix_value();
  dim_vector dm_rows=rows.dims();
  int nvec=dm_rows.elem(0);
  assert(dm_rows.elem(1)==mytod->ndet);  //no cut detectors!
  dim_vector dm_out(mytod->ndata,nvec);
  Matrix vecs(dm_out);
  double *vecs_ptr=vecs.fortran_vec();
  double *rows_ptr=rows.fortran_vec();
  dgemm('N','T',mytod->ndata,nvec,mytod->ndet,1.0,mytod->data[0],mytod->ndata,rows_ptr,nvec,0.0,vecs_ptr,mytod->ndata);
#endif
  return octave_value(rows);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (tod2vecs, args, nargout, "Multiply a row by a TOD, and put it into a column matrix.   Args are (tod,row_vecs)\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));

  if (!mytod->have_data) {
    printf("no data in tod!\n");
    return octave_value_list();
  }


  Matrix rows=args(1).matrix_value();
  dim_vector dm_rows=rows.dims();
  
  int ngood=how_many_dets_are_kept_c(mytod);
  if (dm_rows.elem(1)!=  ngood) {
    printf("Good detector mismatch in vecs2tod %d %d\n",dm_rows.elem(1),ngood);
    return octave_value_list();
  }


  int nvec=dm_rows.elem(0);
  
  Matrix cols(mytod->ndata,nvec);
  
  
  actData **mycols=matrix(mytod->ndata,nvec);
  actData **myrows=matrix(nvec,mytod->ndet);
  memset(mycols[0],0,mytod->ndata*nvec*sizeof(actData));
  memset(myrows[0],0,mytod->ndet*nvec*sizeof(actData));

  int icur=0; 
  double *dd=rows.fortran_vec();
  for (int i=0;i<ngood;i++) {    
    while ((mbCutsIsAlwaysCut(mytod->cuts,mytod->rows[icur],mytod->cols[icur]))&&(icur<mytod->ndet))
      icur++;
    for (int j=0;j<nvec;j++) {
      myrows[j][icur]=dd[j+i*nvec];
    }
    icur++;
  }
  
  
  int nndet=mytod->ndet;
  int nndata=mytod->ndata;


#pragma omp parallel shared(mytod,myrows,mycols,nvec,nndet,nndata) default(none)
  {
    actData **ompcols=matrix(mytod->ndata,nvec);
    memset(ompcols[0],0,mytod->ndata*nvec*sizeof(actData));
    if (nvec>1) {
#pragma omp for
      for (int i=0;i<nndet;i++) {
	for (int j=0;j<nndata;j++)
	  for (int k=0;k<nvec;k++)
	    ompcols[j][k]+=myrows[k][i]*mytod->data[i][j];
      }
    }
    else {
#pragma omp for
      for (int i=0;i<nndet;i++)
	for (int j=0;j<nndata;j++)
	  ompcols[0][j]+=mytod->data[i][j]*myrows[0][i];
    }
#pragma omp critical
    {
      if (nvec>1) {
	for (int i=0;i<nndata;i++)
	  for (int k=0;k<nvec;k++)
	    mycols[i][k]+=ompcols[i][k];
      }
      else {
	for (int i=0;i<nndata;i++)
	  mycols[0][i]+=ompcols[0][i];      
      }
    }
    free(ompcols[0]);
    free(ompcols);
  }

  dd=cols.fortran_vec();
  for (int i=0;i<nndata;i++)
    for (int k=0;k<nvec;k++)
      dd[i*nvec+k]=mycols[i][k];
  
  
  free(mycols[0]);
  free(mycols);
  free(myrows[0]);
  free(myrows);
  
  return octave_value(cols);

}



/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_band_fft, args, nargout, "Get the fft of a tod inside of a band.  Args are (tod,nu1,nu2)\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data) {
    printf("No data present in tod.\n");
    return octave_value_list();
  }
  
  int ndet=mytod->ndet;
  int ndata=mytod->ndata;

  actData dnu=1.0/(mytod->deltat*mytod->ndata);

  actData  nu1,nu2;
  if (args.length()>2)
    nu2=get_value(args(2));
  else
    nu2=dnu*ndata/2;
  if (args.length()>1)
    nu1=get_value(args(1));
  else
    nu1=dnu;
  if (nu2>dnu*ndata/2) {
    nu2=dnu*ndata/2;
    printf("nu2 overly large, setting to %12.4f\n",nu2);
  }
  int i1=(int)(nu1/dnu);
  int i2=(int)(nu2/dnu);

  int nfreq=(i2-i1+1);
      dim_vector dm(nfreq,ndet);
  ComplexMatrix mat(dm);
  memset(mat.fortran_vec(),0,sizeof(Complex)*dm(0)*dm(1));
  
  
  actComplex *mycvec=(actComplex *)malloc(sizeof(actComplex)*ndata);
  actData *myvec=(actData *)malloc(sizeof(actData)*ndata);
  
  act_fftw_plan p=act_fftw_plan_dft_r2c_1d(ndata,myvec,mycvec,FFTW_ESTIMATE);
  
  free(mycvec);
  free(myvec);
  
  Complex *ptr=mat.fortran_vec();    
  double *pptr=(double *)ptr;
#pragma omp parallel private(mycvec) shared(mytod,pptr,ndet,ndata,p,i1,i2,nfreq) default(none)
  {
    mycvec=(actComplex *)malloc(sizeof(actComplex)*ndata);
    actData *mypptr=(double *)mycvec;
#pragma omp for    
    for (int i=0;i<ndet;i++) {
      memset(mypptr,0,sizeof(actComplex)*ndata);
      fftw_execute_dft_r2c(p,mytod->data[i],mycvec);	
      for (int k=0;k<2*nfreq;k++) {
	pptr[i*nfreq*2 +k]=mypptr[2*i1+k];
      }
    }    
    free(mycvec);    
  }  
  act_fftw_destroy_plan(p);


  if (nargout>1) {
    dm(0)=nfreq;
    dm(1)=1;
    Matrix nu(dm);
    for (int i=0;i<nfreq;i++)
      nu(i)=dnu*(i+i1);
    octave_value_list retval;
    retval(0)=octave_value(mat);
    retval(1)=octave_value(nu);
    return retval;
  }
  else
    return octave_value(mat);
  


}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (rotate_tod_data_c, args, nargout, "Rotate the data by a matrix.  If none specified, use internal correlation matrix.  Args are (tod, trans, [mat]))\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  char trans;
  if (args.length()==1)
    trans='n';
  else {
    char *cc=get_char_from_arg(args(1).char_matrix_value());
    trans=cc[0];
    free(cc);
  }
  
  if (args.length()==3) {
#ifndef ACTDATA_DOUBLE  //more work if doing single precision.
    assert(1==0);
#endif
    Matrix mat=args(2).matrix_value();
    rotate_data(mytod,trans,mat.fortran_vec());
  }
  else
    rotate_data(mytod,trans,NULL);

  return octave_value_list();

  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_data_correlation_c, args, nargout, "Get the correlation matrix of a tod..  Arg is (tod))\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data) {
    fprintf(stderr,"Data not present in get_data_correlation.\n");
    return octave_value_list();
  }
  dim_vector dm(2);
  dm(0)=mytod->ndet;
  dm(1)=mytod->ndet;
  Matrix corrs(dm);
  
  get_data_corrs(mytod);
#ifdef ACTDATA_DOUBLE
  //printf("doing double copy.\n");
  memcpy(corrs.fortran_vec(),mytod->corrs[0],sizeof(double)*mytod->ndet*mytod->ndet);
#else
  double *ptr=corrs.fortran_vec();
  for (int i=0;i<mytod->ndet;i++)
    for (int j=i;j<mytod->ndet;j++)
      prt[i*mytod->ndet+j]=mytod->corrs[i][j];
#endif
  free(mytod->corrs[0]);
  free(mytod->corrs);
  mytod->corrs=NULL;


  return octave_value(corrs);
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (assign_rotmat_c, args, nargout, "Assign a rotation matrix to a tod.  Args are (tod,rotmat))\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix rotmat=args(1).matrix_value();
  dim_vector dm=rotmat.dims();
  if ((dm(0)!=mytod->ndet) ||(dm(1)!=mytod->ndet)) {
    mprintf(stdout,"error in assign_rotmat_c.  Size mismatch.\n");
    mytod->rotmat=NULL;
    return octave_value_list();
  }
  mytod->rotmat=matrix(mytod->ndet,mytod->ndet);
  double *mm=rotmat.fortran_vec();
  for (int i=0;i<mytod->ndet*mytod->ndet;i++) {
    mytod->rotmat[0][i]=mm[i];
  }
  
  return octave_value_list();
  
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (scale_tod_band_noise_c, args, nargout, "Set the noise in a TOD.  Args are (tod,white,knee,powlaw))\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int which_band=(int)get_value(args(1));
  actData fac;
  if (args.length()>2)
    fac=get_value(args(2));
  else
    fac=0;
  scale_banded_noise_band(mytod,which_band,fac);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_noise_c, args, nargout, "Set the noise in a TOD.  Args are (tod,white,knee,powlaw))\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  actData white=get_value(args(1));
  actData knee=get_value(args(2));
  actData powlaw=get_value(args(3));

  set_tod_noise(mytod,white,knee,powlaw);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_noise_to_tod, args, nargout, "Add noise to a TOD.  Arg is (tod).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  add_noise_to_tod_new(mytod);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (add_noise_to_tod_gaussian, args, nargout, "Add noise to a TOD.  Arg is (tod).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  add_noise_to_tod_gaussian(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (createFFTWplans_c, args, nargout, "Setup the plans for fft's.  Arg is (tod).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  createFFTWplans1TOD(mytod);
  return octave_value_list();

  
}



/*--------------------------------------------------------------------------------*/
DEFUN_DLD (deconvolve_tod_time_constants_c, args, nargout, "Apply the time constants to a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  nkDeconvolveTimeConstants(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (reconvolve_tod_time_constants_c, args, nargout, "Apply the time constants to a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  nkReconvolveTimeConstants(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (assign_tod_time_constants_c, args, nargout, "Assign the time constants to a TOD.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix taus=args(1).matrix_value();
  dim_vector dm=taus.dims();
  int nrow=dm(0);
  int ncol=dm(1);
  mytod->time_constants=matrix(nrow,ncol);
  for (int i=0;i<nrow;i++)
    for (int j=0;j<nrow;j++)
      mytod->time_constants[i][j]=taus(i,j);
  return octave_value_list();
  
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (apply_real_filter_to_data, args, nargout, "Apply a real-valued filter to the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix filt=args(1).matrix_value();
  actData *ptr=filt.fortran_vec();
  apply_real_filter_to_data(mytod,ptr);
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (apply_complex_filter_to_data, args, nargout, "Apply a real-valued filter to the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  ComplexMatrix filt=args(1).complex_matrix_value();
  actComplex *ptr=(actComplex *)filt.fortran_vec();
  apply_complex_filter_to_data(mytod,ptr);
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (debutterworth_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  nkDeButterworth(mytod);
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (rebutterworth_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  nkReButterworth(mytod);
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_data_fft_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actComplex **data_fft=fft_all_data(mytod);
  
  dim_vector dm(get_nn(mytod->ndata),mytod->ndet);
  ComplexMatrix datft(dm);
  memcpy(datft.fortran_vec(),data_fft[0],get_nn(mytod->ndata)*mytod->ndet*sizeof(actComplex));
  free(data_fft[0]);
  free(data_fft);

  return octave_value(datft);

}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (push_data_fft_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  ComplexMatrix mat=args(1).complex_matrix_value();
  
  if (!mytod->have_data) {
    fprintf(stderr,"allocating storage for TODs.\n");
    allocate_tod_storage(mytod);
  }
  
  actComplex *vec=(actComplex *)mat.fortran_vec();
  dim_vector dm=mat.dims();
  printf("dms are %d %d\n",dm(0),dm(1));
  actComplex **mm=(actComplex **)malloc(sizeof(actComplex *)*dm(0));
  int dd=dm(1);
  int nd=dm(0);
  for (int i=0;i<nd;i++) {
    mm[i]=vec+dd*i;
  }
  ifft_all_data(mytod,mm);
  free(mm);
  return octave_value_list();

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (check_data_fft_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actComplex **data_fft=fft_all_data(mytod);
  memset(mytod->data[0],0,sizeof(actData)*mytod->ndet*mytod->ndata);
  ifft_all_data(mytod,data_fft);
  free(data_fft[0]);
  free(data_fft);

  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (allocate_tod_noise_bands_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix bands=args(1).matrix_value();
  //if (mytod->noise)
  // printf("have noise in mbTOD\n");
  //else
  //printf("do not have noise in mbTOD\n");


  int nband=bands.length()-1;

  double *b=bands.fortran_vec();

  //printf("have %d bands.\n",nband);
  //for (int i=0;i<nband;i++)
  //printf("band %d goes from %12.4f to %12.4f\n",i,b[i],b[i+1]);

  allocate_tod_noise_bands(mytod,b,nband);

  return octave_value_list();
}


/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_simple_banded_noise_model_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix do_rots_octave=args(1).matrix_value();
  Matrix band_types=args(2).matrix_value();
  mbNoiseVectorStructBands *noise=mytod->band_noise;

  double *bb=band_types.fortran_vec();

  mbNoiseType *types=(mbNoiseType *)malloc(sizeof(mbNoiseType)*noise->nband);
  bool *do_rots=(bool *)malloc(sizeof(bool)*noise->nband);

  for (int i=0;i<noise->nband;i++) {
    //printf("band limits on %d are %5d %5d\n",i,noise->ibands[i],noise->ibands[i+1]);
    if (do_rots_octave(i))
      do_rots[i]=true;
    else
      do_rots[i]=false;
    int bbi=bb[i];
    switch(bbi) {
    case 0:
      types[i]=MBNOISE_CONSTANT;
      break;
    case 1:
      types[i]=MBNOISE_FULL;
      break;
    case 2:
      types[i]=MBNOISE_INTERP;
      break;
    default:
      printf("unrecognized type of noise - %d.  Assuming constant.\n",bbi);
      types[i]=MBNOISE_CONSTANT;
      break;      
    }
    
  }
  
  get_simple_banded_noise_model(mytod,do_rots,types);
  //get_simple_banded_noise_model_onerotmat(mytod,do_rots,types);

  dim_vector dm(mytod->ndet,1);
  Matrix facs(dm);
  for (int i=0;i<mytod->ndet;i++)
    facs(i,0)=mytod->band_noise->noise_params[0][i].noise_data[0];
  
  free(do_rots);
  free(types);
  //return octave_value_list();
  return octave_value(facs);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (test_apply_banded_rotations_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actComplex **matft=fft_all_data(mytod);
  int nn=get_nn(mytod->ndata);



  actComplex **mat2=apply_banded_rotations(mytod,matft,true);

  free(matft[0]);
  free(matft);
  

  actComplex **mat3=apply_banded_rotations(mytod,mat2,false);
  free(mat2[0]);
  free(mat2);
  ifft_all_data(mytod,mat3);
  free(mat3[0]);
  free(mat3);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_banded_data_corr_c, args, nargout, "Debutterworth the data.\n")

{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData nu_low=get_value(args(1));
  actData nu_high=get_value(args(2));


  actComplex **data_fft=fft_all_data(mytod);

  actData **mat=get_banded_correlation_matrix_from_fft(mytod,data_fft,nu_low,nu_high);
  
  dim_vector dm(mytod->ndet,mytod->ndet);
  Matrix mycorr(dm);
  memcpy(mycorr.fortran_vec(),mat[0],sizeof(actData)*mytod->ndet*mytod->ndet);
 
  free(mat[0]);
  free(mat);
  free(data_fft[0]);
  free(data_fft);

  return octave_value(mycorr);

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (apply_banded_noise_model_c, args, nargout, "Debutterworth the data.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  apply_banded_noise_model(mytod);
  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (apply_tod_noise_model_c, args, nargout, "Apply the noise.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  apply_noise(mytod);

  return octave_value_list();
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_noise_powlaw_c, args, nargout, "Apply the noise.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  Matrix amps=args(1).matrix_value();
  Matrix knees=args(2).matrix_value();
  Matrix pows=args(3).matrix_value();
  set_noise_powlaw(mytod,amps.fortran_vec(),knees.fortran_vec(),pows.fortran_vec());

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (multiply_all_data, args, nargout, "Multiply the data in a TOD by a constant.  Args are (tod,val).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData val=get_value(args(1));
  multiply_all_data(mytod,val);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (test_eig_c, args, nargout, "Multiply the data in a TOD by a constant.  Args are (tod,val).\n")
{
  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();

  double *vec=mat.fortran_vec();
  get_eigenvectors(&vec,dm(0));
  return octave_value(mat);
  

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (shift_tod_pointing_c, args, nargout, "Shift the TOD pointing in RA and Dec.  args are (tod,ra_shift,dec_shift).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  actData ra_shift=get_value(args(1));
  actData dec_shift=get_value(args(2));
  
  shift_tod_pointing(mytod,ra_shift,dec_shift);

  return octave_value_list();
  

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_radec_lims_c, args, nargout, "Return the tod pointing limits.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  dim_vector dm(1,4);
  Matrix mat(dm);
  mat(0,0)=mytod->ramin;
  mat(0,1)=mytod->ramax;
  mat(0,2)=mytod->decmin;
  mat(0,3)=mytod->decmax;

  return octave_value(mat);

}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (get_tod_altaz_lims_c, args, nargout, "Return the tod pointing limits.\n")
{

  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  dim_vector dm(1,4);
  Matrix mat(dm);

  actData altmin=mytod->alt[0];
  actData altmax=mytod->alt[0];

  actData azmin=mytod->az[0];
  actData azmax=mytod->az[0];

  for (int i=0;i<mytod->ndata;i++) {
    if (mytod->alt[i]>altmax)
      altmax=mytod->alt[i];
    if (mytod->az[i]>azmax)
      azmax=mytod->az[i];

    if (mytod->alt[i]<altmin)
      altmin=mytod->alt[i];
    if (mytod->az[i]<azmin)
      azmin=mytod->az[i];


  }


  mat(0,0)=altmin;
  mat(0,1)=altmax;
  mat(0,2)=azmin;
  mat(0,3)=azmax;

  return octave_value(mat);

}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (erase_tod_c, args, nargout, "Free as much of a TOD as we can ")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (mytod->have_data) {
    free(mytod->data[0]);
    free(mytod->data);
  }
  free(mytod->az);
  free(mytod->alt);
  if (mytod->ra)
    free(mytod->ra);
  if (mytod->dec)
    free(mytod->dec);
  if (mytod->corrs) {
    free(mytod->corrs[0]);
    free(mytod->corrs);
  }
  if (mytod->rotmat) {
    free(mytod->rotmat[0]);
    free(mytod->rotmat);
  }
  if (mytod->cuts)
    CutsFree(mytod->cuts);


  //don't actually free the TOD since it may have been allocated in a vector.
  //free(mytod);

  return octave_value_list();

    


}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_radec_from_altaz_ctime_c, args, nargout, "Get RA/Dec of a set of points.\n")
{

  Matrix alt_mat=args(0).matrix_value();
  Matrix az_mat=args(1).matrix_value();
  Matrix ctime_mat=args(2).matrix_value();
  double *alt=alt_mat.fortran_vec();
  double *az=az_mat.fortran_vec();
  double *ctime=ctime_mat.fortran_vec();
  dim_vector dm=alt_mat.dims();
  
  Matrix ra_mat(dm);
  Matrix dec_mat(dm);
  double *ra=ra_mat.fortran_vec();
  double *dec=dec_mat.fortran_vec();  

  int ndata=dm(0)*dm(1);

  if (ndata<1)
    return octave_value_list();
  

  Site site;
  ACTSite(&site);
  //printf("ctime, alt, az of first sample are %14.4g %14.6f %14.6f\n",ctime[0],alt[0],az[0]);
  act_observed_altaz_to_mean_radec(&site,150.0,ndata,ctime,alt,az,ra,dec);  
  //printf("ra/dec of first sample are %14.6f %14.6f\n",ra[0],dec[0]);
  octave_value_list retval;
  retval(0)=octave_value(ra_mat);
  retval(1)=octave_value(dec_mat);
  return retval;
}
 
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_detector_radec_c, args, nargout, "Get RA/Dec of a detector.\n")
{
  int nargin = args.length();
  if (nargin<2) {
    printf("Only have %d arguments in get_detector_radec_c\n",nargin);
    return octave_value_list();
  }

  bool do_exact=false;
  if (nargin>2)
    do_exact=true;
	      
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int det=(int)get_value(args(1));
  PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
  if (do_exact)
    get_radec_from_altaz_exact_1det(tod,det,scratch);
  else {
    //printf("Getting detector fit coarse.\n");
    get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
  }  

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
DEFUN_DLD (get_tod_coarse_ind_c, args, nargout, "Get indices of pivots for the coarse pointing.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix pivots(tod->pointing_fit->ncoarse,1);
  for (int i=0;i<tod->pointing_fit->ncoarse;i++) {
    pivots(i,0)=tod->pointing_fit->coarse_ind[i];
  }
  return octave_value(pivots);
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_radec_from_altaz_fit_tiled_c, args, nargout, "Evaluate pointing from a tiled fit.\n")
{
  int nargin = args.length();
  if (nargin<4) {
    printf("need 4 arguments in get_radec_from_altaz_fit_tiled_c (tod, alt_vec, az_vec, time_vec)\n");
    return octave_value_list();
  }
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix alt_vec=args(1).matrix_value();
  Matrix az_vec=args(2).matrix_value();
  Matrix time_vec=args(3).matrix_value();
  dim_vector dm=alt_vec.dims();
  int nelem=matrix_nelem(alt_vec);
  printf("applying pointing to %d elements.\n",nelem);
  double *alt=alt_vec.fortran_vec();
  double *az=az_vec.fortran_vec();
  double *mytime=time_vec.fortran_vec();
  Matrix ra_vec(dm(0),dm(1));
  Matrix dec_vec(dm(0),dm(1));
  double *ra=ra_vec.fortran_vec();
  double *dec=dec_vec.fortran_vec();
  get_radec_from_altaz_fit_tiled(tod->pointing_fit->tiled_fit,alt,az,mytime,ra,dec,nelem);
  octave_value_list retval;
  retval(0)=octave_value(ra_vec);
  retval(1)=octave_value(dec_vec);
  return retval;
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (set_tod_pointing_tiled_c, args, nargout, "Set pointing to a tiled fit.\n")
{
  int nargin = args.length();
  if (nargin<7) {
    printf("need 7 arguments in set_tod_pointing_tiled_c (tod, alt_vec, az_vec, ra_mat, dec_mat, ra_clock, dec_clock)\n");
    return octave_value_list();
  }

  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix alt_vec=args(1).matrix_value();
  Matrix az_vec=args(2).matrix_value();
  Matrix ra_mat=args(3).matrix_value();
  Matrix dec_mat=args(4).matrix_value();
  Matrix ra_clock=args(5).matrix_value();
  Matrix dec_clock=args(6).matrix_value();
  
  double *alt_ptr=alt_vec.fortran_vec();
  int nalt=matrix_nelem(alt_vec);
  double *az_ptr=az_vec.fortran_vec();
  int naz=matrix_nelem(az_vec);
  double *ra_ptr=ra_mat.fortran_vec();
  double *dec_ptr=dec_mat.fortran_vec();
  assert(matrix_nelem(ra_mat)==nalt*naz);
  assert(matrix_nelem(dec_mat)==nalt*naz);
  
  
  double *ra_clock_ptr=ra_clock.fortran_vec();
  double *dec_clock_ptr=dec_clock.fortran_vec();
  set_tod_pointing_tiled(tod,az_ptr,naz,alt_ptr,nalt,&ra_ptr,&dec_ptr, ra_clock_ptr, dec_clock_ptr, matrix_nelem(ra_clock));
  

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (get_all_detector_radec_c, args, nargout, "Get RA/Dec of a detector.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix ra_mat(tod->ndata,tod->ndet);
  Matrix dec_mat(tod->ndata,tod->ndet);
  
  
  //printf("tod->ra_saved is %ld, tod->dec_saved is %ld\n",(long)tod.ra_saved-(long)tod,(long)tod.dec_saved-(long)tod);
  //printf("tod->ra_saved is %ld, tod->dec_saved is %ld\n",(long)(&(tod->ra_saved))-(long)tod,(long)(&(tod->dec_saved))-(long)tod);

  if ((tod->ra_saved!=NULL)&&(tod->dec_saved!=NULL)) {
    //printf("returning saved pointing.\n");
    //printf("ndata and ndet are %d %d\n",tod->ndata,tod->ndet);
    actData *ra_ptr=ra_mat.fortran_vec();
    actData *dec_ptr=dec_mat.fortran_vec();
    memcpy(ra_ptr,tod->ra_saved[0],tod->ndata*tod->ndet*sizeof(actData));
    memcpy(dec_ptr,tod->dec_saved[0],tod->ndata*tod->ndet*sizeof(actData));
  
    octave_value_list retval;
    retval(0)=octave_value(ra_mat);
    retval(1)=octave_value(dec_mat);
    
    return retval;

  }
  


  int nargin = args.length();
  if (nargin<2) {
    printf("Only have %d arguments in get_detector_radec_c\n",nargin);
    return octave_value_list();
  }

  bool do_exact=false;
  if (nargin>2)
    do_exact=true;	      
  
#pragma omp parallel shared(ra_mat,dec_mat,tod,do_exact) default(none)
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    double *ra=ra_mat.fortran_vec();
    double *dec=dec_mat.fortran_vec();

#pragma omp for 
    for (int det=0;det<tod->ndet;det++) {
      if (do_exact)
	get_radec_from_altaz_exact_1det(tod,det,scratch);
      else
	get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
      
      for (int i=0;i<tod->ndata;i++) {
	int ii=i+det*tod->ndata;
	ra[ii]=scratch->ra[i];
	dec[ii]=scratch->dec[i];
      }
    }
    destroy_pointing_fit_scratch(scratch);
  }
  
  
  octave_value_list retval;
  retval(0)=octave_value(ra_mat);
  retval(1)=octave_value(dec_mat);
 
  return retval;

}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (find_max_tod_pointing_err_c, args, nargout, "Find maximum pointing error of a TOD.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  return octave_value(find_max_pointing_err(tod));
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_detector_altaz_c, args, nargout, "Get alt/az of a detector.\n")
{
  int nargin = args.length();
  if (nargin<2) {
    printf("Only have %d arguments in get_detector_radec_c\n",nargin);
    return octave_value_list();
  }

  bool do_exact=false;
  if (nargin>2)
    do_exact=true;
	      
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  int det=(int)get_value(args(1));
  //printf("working on detector %d with row %d and col %d\n",det,tod->rows[det],tod->cols[det]);
  PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
  if (do_exact)
    get_radec_from_altaz_exact_1det(tod,det,scratch);
  else {
    //get_radec_from_altaz_fit_1det_coarse(tod,det,scratch);
    get_radec_from_altaz_fit_1det(tod,det,scratch);
  }  

  Matrix altaz(tod->ndata,2);
  double *dat=altaz.fortran_vec();
  for (int i=0;i<tod->ndata;i++) {
    dat[i]=scratch->alt[i];
    dat[i+tod->ndata]=scratch->az[i];
  }

  destroy_pointing_fit_scratch(scratch);
 
  return octave_value(altaz);

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (convert_saved_pointing_to_pixellization, args, nargout, "Convert a TOD's saved RA/Dec to a map pixellization, freeing the RA/Dec storage..\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  MAP *mymap=(MAP *)get_pointer(args(1));
  convert_saved_pointing_to_pixellization(mytod,mymap);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (free_saved_pixellization, args, nargout, "Free a TOD's saved pixellization.\n")
{
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  if (tod->pixelization_saved) {
    free(tod->pixelization_saved[0]);
    free(tod->pixelization_saved);
    tod->pixelization_saved=NULL;
  }
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_alldet_altaz_lims_c, args, nargout, "Get alt/az limits spanned by detectors in a TOD.\n")
{
	      
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix lims_mat(4,1);
  lims_mat(0)=1000;
  lims_mat(1)=-1000;
  lims_mat(2)=1000;
  lims_mat(3)=-1000;
#pragma omp parallel shared(lims_mat,tod) default(none)
  {
    PointingFitScratch *scratch=allocate_pointing_fit_scratch(tod);
    double *lims=lims_mat.fortran_vec();
    double altmin=lims[0];
    double altmax=lims[1];
    double azmin=lims[2];
    double azmax=lims[3];
    
#pragma omp for
    for (int det=0;det<tod->ndet;det++) {
      get_radec_from_altaz_fit_1det(tod,det,scratch);
      for (int i=0;i<tod->ndata;i++) {
	if (scratch->alt[i]<altmin)
	  altmin=scratch->alt[i];
	if (scratch->alt[i]>altmax)
	  altmax=scratch->alt[i];
	if (scratch->az[i]<azmin)
	  azmin=scratch->az[i];
	if (scratch->az[i]>azmax)
	  azmax=scratch->az[i];
      }
    }
#pragma omp critical
    {
      if (altmin<lims[0])
	lims[0]=altmin;
      if (altmax>lims[1])
	lims[1]=altmax;
      if (azmin<lims[2])
	lims[2]=azmin;
      if (azmax>lims[3])
	lims[3]=azmax;
    }
    
    destroy_pointing_fit_scratch(scratch);
  }
  return octave_value(lims_mat);

}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_median_altaz_c, args, nargout, "Get median alt/az of a tod name.\n")
{

  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();
  
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  double azmed=(double)compute_median_inplace(tod->ndata,tod->az);
  double altmed=(double)compute_median_inplace(tod->ndata,tod->alt);
  double ctime_med=tod->ctime+0.5*tod->deltat*(double)(tod->ndata/2);

  octave_value_list retval;
  retval(0)=azmed;
  retval(1)=altmed;
  retval(2)=ctime_med;
  return retval;
  
  
  
}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (set_tod_pointing_c, args, nargout, "Set the pointing for a TOD.  Return TOD limits\n")
{
  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();
  
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  Matrix dalt=args(1).matrix_value();
  Matrix daz=args(2).matrix_value();
  dim_vector dm=dalt.dims();
  mbPointingOffset *offset;
  //printf("dimenstions are %d %d\n",dm(0),dm(1));
  if (tod->pointingOffset)
    offset=tod->pointingOffset; //assuming this has been allocated correctly
  else {
    offset = nkPointingOffsetAlloc(dm(0),dm(1),0);
    tod->pointingOffset=offset; 
  }
  //printf("allocated.\n");
  for (int i=0;i<dm(0);i++)
    for (int j=0;j<dm(1);j++) {
      offset->offsetAlt[i][j]=dalt(i,j);
      offset->offsetAzCosAlt[i][j]=daz(i,j);      
    }
  //printf("assigned.\n");
  //printf("Offsets of (0,0) are %14.7f %14.7f\n",offset->offsetAlt[0][0],offset->offsetAzCosAlt[0][0]);

  cut_mispointed_detectors(tod);
  //printf("cut mispointed detectors.\n");
  assign_tod_ra_dec(tod);
  //printf("assigned ra_dec.\n");
  
  //printf("Offsets of (0,0) are %14.7f %14.7f\n",offset->offsetAlt[0][0],offset->offsetAzCosAlt[0][0]);

  find_pointing_pivots(tod,0.5);
  //printf("Offsets of (0,0) are %14.7f %14.7f\n",offset->offsetAlt[0][0],offset->offsetAzCosAlt[0][0]);
  find_tod_radec_lims(tod);
  //printf("Offsets of (0,0) are %14.7f %14.7f\n",tod->pointingOffset->offsetAlt[0][0],tod->pointingOffset->offsetAzCosAlt[0][0]);

  Matrix lims(4,1);
  lims(0,0)=tod->ramin;
  lims(1,0)=tod->ramax;
  lims(2,0)=tod->decmin;
  lims(3,0)=tod->decmax;

  return octave_value(lims);
   
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_tod_pointing_offsets_c, args, nargout, "Make a matrix of the TOD pointing offsets.\n")
{

  int nargin = args.length();
  if (nargin==0)
    return octave_value_list();
  
  mbTOD  *tod=(mbTOD *)get_pointer(args(0));
  printf("Offsets of (0,0) are %14.7f %14.7f\n",tod->pointingOffset->offsetAlt[0][0],tod->pointingOffset->offsetAzCosAlt[0][0]);
  Matrix dalt(33,32);
  Matrix daz(33,32);
  for (int i=0;i<33;i++)
    for (int j=0;j<32;j++) {
      dalt(i,j)=-100;
      daz(i,j)=-100;
    }
  
  for (int i=0;i<tod->ndet;i++) {
    int row=tod->rows[i];
    int col=tod->cols[i];
    dalt(row,col)=tod->pointingOffset->offsetAlt[row][col];
    daz(row,col)=tod->pointingOffset->offsetAzCosAlt[row][col];
  }
  octave_value_list retval;
  retval(0)=octave_value(dalt);
  retval(1)=octave_value(daz);
  return retval;
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_alldata_fft_plans_c, args, nargout, "Make FFTW plans for a TOD using FFTW measures  Optionally, import/save them in a wisdom file.\n")
{

  //will check and see if data is allocated coming in.  If it is, don't free it at the end.
  bool had_data=false;
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->have_data)
    allocate_tod_storage(mytod);
  else
    had_data=true;
  
  char *wisdom_name=NULL;
  if (args.length()>1)
    wisdom_name=get_char_from_arg(args(1).char_matrix_value());
  
  if (wisdom_name)
    {
      FILE *wisdomfile=fopen(wisdom_name,"r");
      if (wisdomfile)
	fftw_import_wisdom_from_file(wisdomfile);
      else
	fprintf(stderr,"unable to read wisdom from %s\n",wisdom_name);
    }
  actComplex **data_fft=fft_all_data_flag(mytod,FFTW_MEASURE);
  //printf("did first  bit.\n");
  
  ifft_all_data_flag(mytod,data_fft,FFTW_MEASURE);

  free(data_fft[0]);
  free(data_fft);
  if (!had_data)
    free_tod_storage(mytod);
  
  if (wisdom_name)
    {
      //first, re-read wisdom in case someone else has updated in the meantime
      FILE *wisdomfile=fopen(wisdom_name,"r");
      if (wisdomfile) {
	fftw_import_wisdom_from_file(wisdomfile);
	fclose(wisdomfile);
      }
      wisdomfile=fopen(wisdom_name,"w");
      if (wisdomfile) {
	fftw_export_wisdom_to_file(wisdomfile);
	fclose(wisdomfile);
      }
      else
	fprintf(stderr,"Error writing to wisdomfile %s\n",wisdom_name);
      
    }
  return octave_value_list();
}



/*--------------------------------------------------------------------------------*/
#if 0
DEFUN_DLD (simple_test_diag_proj_noise_inv, args, nargout, "Test some C-version noise inverses.\n")
{
  Matrix data=args(0).matrix_value();
  Matrix noise=args(1).matrix_value();
  Matrix vecs=args(2).matrix_value();



  double **dataptr=matrix2cmat(data.fortran_vec(),data.dims());
  dim_vector dm=data.dims();
  Matrix dinv(dm);
  double **dinvptr=matrix2cmat(dinv.fortran_vec(),dinv.dims());
  double *noiseptr=noise.fortran_vec();
  double *vecfptr=vecs.fortran_vec();
  double **vecptr=matrix2cmat(vecfptr,vecs.dims());
  
  dim_vector vecs_dm=vecs.dims();

  int ndata=dm(0);
  int ndet=dm(1);
  int nvecs=vecs_dm(1);  


  if (args.length()>3) {
    int imin=(int)get_value(args(3));
    int imax=(int)get_value(args(4));
    for (int i=0;i<ndet;i++)
      noiseptr[i]=1.0/noiseptr[i];
    apply_diag_proj_noise_inv_bands(dataptr,dinvptr,noiseptr,vecptr,ndata,ndet,nvecs,imin,imax);
  }
  else
    simple_test_diag_proj_noise_inv(dataptr,dinvptr,noiseptr,vecptr,ndata,ndet,nvecs);


  free(dataptr);
  free(dinvptr);
  free(vecptr);
  return octave_value(dinv);
}
#endif
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (pull_oneband_tod_noise_banded_projvec,args,nargout,"Pull parameters for a single band in a TOD noise + projection vector model.\n")
{
  if (args.length()<2) {
    printf("Error in pull_oneband_tod_noise_banded_projvec, need at least 4 arguments.\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  if (!mytod->band_vecs_noise) {
    fprintf(stderr,"missing band_vecs_noise in tod.\n");
    return octave_value_list();
  }
  int myband=(int)get_value(args(1))-1;
  Matrix det_noise(mytod->ndet,1);
  for (int i=0;i<mytod->ndet;i++)
    det_noise(i,0)=mytod->band_vecs_noise->noises[myband][i];
  Matrix vecs(mytod->band_vecs_noise->nvecs[myband],mytod->ndet);
  for (int i=0;i<mytod->band_vecs_noise->nvecs[myband];i++)
    for (int j=0;j<mytod->ndet;j++)
      vecs(i,j)=mytod->band_vecs_noise->vecs[myband][i][j];
  octave_value_list retval;
  retval(0)=det_noise;
  retval(1)=vecs;
  return retval;

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (set_oneband_tod_noise_banded_projvec,args,nargout,"Set up a single band in a TOD noise + projection vector model.\n")
{
  if (args.length()<4) {
    printf("Error in set_oneband_tod_noise_banded_projvec, need at least 4 arguments.\n");
    return octave_value_list();
  }
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  int myband=(int)get_value(args(1))-1;
  Matrix noises=args(2).matrix_value();
  Matrix vecs=args(3).matrix_value();
  
  dim_vector dm=vecs.dims();
  //printf("dims are %d %d\n",dm(0),dm(1));
  
  actData *nvec=noises.fortran_vec();
  for (int i=0;i<mytod->ndet;i++)
    mytod->band_vecs_noise->noises[myband][i]=nvec[i];

  double *vecptr=vecs.fortran_vec();
  mytod->band_vecs_noise->nvecs[myband]=dm(1);

  mytod->band_vecs_noise->vecs[myband]=matrix(mytod->band_vecs_noise->nvecs[myband],mytod->ndet);

  
  for (int i=0;i<mytod->band_vecs_noise->nvecs[myband];i++)
    for (int j=0;j<mytod->ndet;j++)
      mytod->band_vecs_noise->vecs[myband][i][j]=vecptr[i*mytod->ndet+j];
  

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (save_tod_noise_banded_projvec,args,nargout,"Save the banded_projvec noise from a TOD.  Args are (tod,filename).\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  char *fname=get_char_from_arg(args(1).char_matrix_value());
  FILE *outfile=fopen(fname,"w");
  if (!outfile) {
    fprintf(stderr,"Unable to open %s for writing in save_tod_noise_banded_projvec.\n",fname);
    return octave_value_list();
  }
  mbNoiseStructBandsVecs *noise=mytod->band_vecs_noise;
  fwrite(&(noise->ndet),sizeof(noise->ndet),1,outfile);
  fwrite(&(noise->nband),sizeof(noise->nband),1,outfile);
  fwrite((noise->band_edges),sizeof(*(noise->band_edges)),1+noise->nband,outfile);
  fwrite(noise->nvecs,sizeof(*(noise->nvecs)),noise->nband,outfile);
  fwrite(noise->noises[0],sizeof(actData),noise->nband*noise->ndet,outfile);
  for (int i=0;i<noise->nband;i++) {
    fwrite(noise->vecs[i][0],sizeof(actData),noise->nvecs[i]*noise->ndet,outfile);
  }
  
  //mytod->band_vecs_noise->noises=matrix(nband,mytod->ndet);
  //mytod->band_vecs_noise->vecs=(actData ***)malloc(sizeof(actData **)*nband);
  
  
  fclose(outfile);
  return octave_value_list();
  

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (read_tod_noise_banded_projvec,args,nargout,"Read the banded_projvec noise into a TOD.  Args are (tod,filename).\n")
{
  //fprintf(stderr,"hello!\n");
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  //fprintf(stderr,"god tod pointer.\n");
  char *fname=get_char_from_arg(args(1).char_matrix_value());
  //fprintf(stderr,"fname is %s\n",fname);
  FILE *infile=fopen(fname,"r");
  //fprintf(stderr,"opened infile.\n");
  if (!infile) {
    fprintf(stderr,"Could not find file %s for reading in read_tod_noise_banded_projvec\n",fname);
    return octave_value_list();
  }
    
  mbNoiseStructBandsVecs *noise=(mbNoiseStructBandsVecs *)malloc(sizeof(mbNoiseStructBandsVecs));
  //printf("here 1.\n");
  fread(&(noise->ndet),sizeof(int),1,infile);
  fread(&(noise->nband),sizeof(int),1,infile);
  //printf("ndet and nband are %d %d\n",noise->ndet,noise->nband);
  noise->band_edges=(int *)malloc(sizeof(int)*(1+noise->nband));
  fread((noise->band_edges),sizeof(*(noise->band_edges)),1+noise->nband,infile);
  //printf("read edges.\n");

  noise->nvecs=(int *)malloc(sizeof(int)*noise->nband);
  fread(noise->nvecs,sizeof(*(noise->nvecs)),noise->nband,infile);
  //printf("read nvecs.\n");

  noise->noises=matrix(noise->nband,noise->ndet);
  fread(noise->noises[0],sizeof(actData),noise->nband*noise->ndet,infile);
  //printf("read detector noises.\n");
  
  noise->vecs=(actData ***)malloc(sizeof(actData **)*noise->nband);
  for (int i=0;i<noise->nband;i++) {
    noise->vecs[i]=matrix(noise->nvecs[i],noise->ndet);
    fread(noise->vecs[i][0],sizeof(actData),noise->nvecs[i]*noise->ndet,infile);
  }
  //printf("read vecs.\n");
  fclose(infile);
  mytod->band_vecs_noise=noise;
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (allocate_tod_noise_banded_projvec,args,nargout,"Set the TOD noise model to be in bands with non-completely projected vectors.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  
  Matrix bands=args(1).matrix_value();
  dim_vector dm=bands.dims();
  int nband=dm(0);
  if (dm(1)>nband)
    nband=dm(1);
  nband--;

  mytod->band_vecs_noise=(mbNoiseStructBandsVecs *)malloc(sizeof(mbNoiseStructBandsVecs));

  mytod->band_vecs_noise->ndet=mytod->ndet;
  mytod->band_vecs_noise->nband=nband;
  mytod->band_vecs_noise->band_edges=(int *)malloc(sizeof(int)*(nband+1));
  double *bb=bands.fortran_vec();
  for (int i=0;i<nband+1;i++) {
    mytod->band_vecs_noise->band_edges[i]=(int)bb[i];
    //printf("assigned band edge %2d %5d\n",i,mytod->band_vecs_noise->band_edges[i]);
  }
  mytod->band_vecs_noise->nvecs=(int *)malloc(sizeof(int)*(nband));
  

  mytod->band_vecs_noise->noises=matrix(nband,mytod->ndet);
  mytod->band_vecs_noise->vecs=(actData ***)malloc(sizeof(actData **)*nband);
  
  
  return octave_value_list();

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (query_tod_type,args,nargout,"Query what type of TOD you have.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));
  return octave_value(mytod->todtype);
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (insert_dataset_into_tod,args,nargout,"Put a pointer into a TOD structure.  Args are (structure, todtype, tod)\n")
{
  void  *myptr=(void *)get_pointer(args(0));
  int mytype=get_value(args(1));
  mbTOD  *mytod=(mbTOD *)get_pointer(args(2));
  mytod->todtype=mytype;
  mytod->generic=myptr;

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (get_generic_tod_pointer,args,nargout,"Pull the generic pointer out of a TOD structure.\n")
{
  mbTOD  *mytod=(mbTOD *)get_pointer(args(0));

  int64NDArray myptr(1);
  long myptr_asint=(long)mytod->generic;
  myptr(0)=myptr_asint;

  octave_value_list retval;
  retval(0)=myptr;
  retval(1)=mytod->todtype;

  return retval;
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (test_legendre_fit,args,nargout,"Do a Legendre polynomial fit from ninkasi.  args are (vec,order)")
{
  Matrix data=args(0).matrix_value();
  dim_vector dm=data.dims();
  int order=get_value(args(1));
  Matrix fitp(order,1);
  Matrix projp(order,1);
  int ndata=dm(0)*dm(1);
  actData *tmp=legendre_fit(data.fortran_vec(),ndata,order);
#if 0
  int i;
  for (i=0;i<order;i++)
    tmp[i]=0;
  tmp[1]=1;
#endif
  memcpy(fitp.fortran_vec(),tmp,order*sizeof(actData));
  free(tmp);
  Matrix fitdata(dm(0),dm(1));
  legendre_eval(fitdata.fortran_vec(),ndata,fitp.fortran_vec(),order);
  memset(projp.fortran_vec(),0,sizeof(double)*order);
  legendre_project(fitdata.fortran_vec(),ndata,projp.fortran_vec(),order);
  octave_value_list retval;
  retval(0)=octave_value(fitp);
  retval(1)=octave_value(fitdata);
  retval(2)=octave_value(projp);
  return retval;
}
