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
#include <qpoint/qpoint.h>

#ifdef __cplusplus
}  /* end extern "C" */
#endif


//from sasha on how to get boresight stuff set up
//[2015-05-21, 11:05:22 PM] Sasha Rahlin: so either you read RA/DEC/PA from the point00 product
//[2015-05-21, 11:05:25 PM] Sasha Rahlin: and call radecpa2quat
//[2015-05-21, 11:05:34 PM] Sasha Rahlin: or read QBORE[0-3] from the quat00 product
//to get a detector offset as a quaternion:  theres a qp_det_offset() function, or qp_det_offsetn for N detectors
//options can get set via qp_set_options, after calling qp_init_memory
//qp_set_opt_fast_math(mem, 1) to use poly approximations


/*--------------------------------------------------------------------------------*/

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=(double)myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/


DEFUN_DLD (qp_init_c, args, nargout, "Initialize qpoint memory.\n")
{
  qp_memory_t *qpt=  qp_init_memory();
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
DEFUN_DLD (qp_boresight_pointing_c, args, nargout, "Get boresight pointing from alt/az.\n")
{
  qp_memory_t *mem=  qp_init_memory();

  ColumnVector az=args(0).column_vector_value();
  int n=az.length();
  //printf("going to process %d samples.\n",n);
  ColumnVector el=args(1).column_vector_value();
  ColumnVector pitch=args(2).column_vector_value();
  ColumnVector roll=args(3).column_vector_value();
  ColumnVector lon=args(4).column_vector_value();
  ColumnVector lat=args(5).column_vector_value();
  ColumnVector ct=args(6).column_vector_value();
  
  

  //for (int i=0;i<n;i++) {
  //  printf("%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %15.3f\n",az(i),el(i),pitch(i),roll(i),lat(i),lon(i),ct(i));
  // }

  quat_t *q_bore = (quat_t *)malloc(sizeof(quat_t)*n);
  //memset(q_bore,0,sizeof(quat_t)*n);
  //after discussions with sasha, q_bore needs to be initialized, since one can put offsets in it.  default is (1,0,0,0) for no offsets.
  for (int i=0;i<n;i++) {
    q_bore [i][0]=1.0;
    q_bore[i][1]=0.0;
    q_bore[i][2]=0.0;
    q_bore[i][3]=0.0;
  }
  
  //printf("sizes are %d %d\n",sizeof(double),sizeof(quat_t));

  qp_azel2bore(mem, az.fortran_vec(), el.fortran_vec(), pitch.fortran_vec(), roll.fortran_vec(), lon.fortran_vec(), lat.fortran_vec(), ct.fortran_vec(), q_bore, n);

  quat_t qoff={1.0,0,0,0};
  ColumnVector ra(n);
  ColumnVector dec(n);
  ColumnVector psi(n);
  memset(ra.fortran_vec(),0,sizeof(double)*n);
  memset(dec.fortran_vec(),0,sizeof(double)*n);
  memset(psi.fortran_vec(),0,sizeof(double)*n);
  


  //qp_bore2radec(mem, delta_az, delta_el, delta_psi, ct.fortran_vec(), q_bore,ra.fortran_vec(), dec.fortran_vec(), sin2psi.fortran_vec(), cos2psi.fortran_vec(), n);
  //qp_bore2radec(mem, qoff, ct.fortran_vec(), q_bore,ra.fortran_vec(), dec.fortran_vec(), sin2psi.fortran_vec(), cos2psi.fortran_vec(), n);
  //void qp_quat2radecpan(qp_memory_t *mem, quat_t *q, double *ra, double *dec,double *pa, int n);
  qp_quat2radecpan(mem,q_bore,ra.fortran_vec(),dec.fortran_vec(),psi.fortran_vec(),n);


  free(q_bore);

  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;
  retval(2)=psi;
  
  qp_free_memory(mem);

  return retval;
  
}
/*--------------------------------------------------------------------------------*/

DEFUN_DLD(qp_detector_pointing_c, args, nargout, "Get a detector pointing offset from boresight RA/Dec.\n")
{
  ColumnVector ra=args(0).column_vector_value();
  ColumnVector dec=args(1).column_vector_value();
  ColumnVector pa=args(2).column_vector_value();
  ColumnVector ct=args(3).column_vector_value();
  double daz=get_value(args(4));
  double del=get_value(args(5));
  
  int n=ra.length();
  if (dec.length()!=n) {
    printf("size mismatch between ra and dec in qp_detector_pointing_c.\n");
    return octave_value_list();
  }
  if (pa.length()!=n) {
    printf("size mismatch between ra and pa in qp_detector_pointing_c.\n");
    return octave_value_list();
  }

  qp_memory_t *qpt=  qp_init_memory();  
  quat_t *q=(quat_t *)malloc(sizeof(quat_t)*n);
  qp_radecpa2quatn(qpt, ra.fortran_vec(), dec.fortran_vec(),pa.fortran_vec(),q,n);
  
  quat_t detq;
  qp_det_offset(daz,del,0.0,detq);
  ColumnVector sin2psi(n);
  ColumnVector cos2psi(n);
  //printf("calling bore2radec.\n");
  qp_bore2radec(qpt,detq,ct.fortran_vec(),q,ra.fortran_vec(),dec.fortran_vec(),sin2psi.fortran_vec(),cos2psi.fortran_vec(),n);
  //printf("back.\n");
  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;
  retval(2)=sin2psi;
  retval(3)=cos2psi;

  return retval;

}
