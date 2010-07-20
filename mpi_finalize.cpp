#include <mpi.h>
#include <iostream>
#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;


DEFUN_DLD (mpi_finalize, args, nargout,
	   "Finish up MPI.  Please only call this once...\n")
{
  int argc=0;
  MPI_Finalize ();
  return octave_value_list ();
}

