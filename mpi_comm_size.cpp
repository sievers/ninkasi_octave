#include <mpi.h>
#include <iostream>
#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
#include <dMatrix.h>
using namespace std;


DEFUN_DLD (mpi_comm_size, args, nargout,
	   "Return total number of proceses.\n")
{
  int size;
  MPI_Comm_size (MPI_COMM_WORLD,&size);
  //printf("There are a total of  %d\n",size);
  return octave_value(size);
}

