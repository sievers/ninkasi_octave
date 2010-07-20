#include <mpi.h>
#include <iostream>
#include <octave/oct.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;


DEFUN_DLD (mpi_comm_rank, args, nargout,
	   "Return which process you are.\n")
{
  int flag=0;
  MPI_Initialized(&flag);


  int rank=0;
  if (flag)
    MPI_Comm_rank (MPI_COMM_WORLD,&rank);
  //printf("I am %d\n",rank);
  return octave_value(rank);

  //return octave_value_list ();
}

