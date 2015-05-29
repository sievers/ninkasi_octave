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
  printf("calling bore2radec.\n");
  qp_bore2radec(qpt,detq,ct.fortran_vec(),q,ra.fortran_vec(),dec.fortran_vec(),sin2psi.fortran_vec(),cos2psi.fortran_vec(),n);
  printf("back.\n");
  octave_value_list retval;
  retval(0)=ra;
  retval(1)=dec;
  retval(2)=sin2psi;
  retval(3)=cos2psi;

  return retval;

}
