#include <octave/oct.h>

#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#ifdef _MKL
#include <mkl.h>
#include <omp.h>
#endif

using namespace std;


DEFUN_DLD (mpi_init, args, nargout,
	   "Fire up MPI.  Please only call this once...\n")
{
  int argc=0;



#ifdef _MKL
  int nthread;
#pragma omp parallel
#pragma omp single
  nthread=omp_get_num_threads();
  //printf("nthread is %d\n",nthread);
#endif

  //Make sure we haven't already called mpi_init!
  int flag=0;
  MPI_Initialized(&flag);
  

  
  if (!flag)
    MPI_Init (&argc,NULL);

#ifdef _MKL
  //mkl_set_num_threads(8);
  omp_set_num_threads(8);
  //printf("set num_threads to %d\n",nthread);
#endif

  return octave_value_list ();
}

