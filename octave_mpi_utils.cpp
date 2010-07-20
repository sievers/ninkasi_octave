#include <octave/oct.h>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include <string.h>
#include <ctype.h>
#ifdef __cplusplus
}  /* end extern "C" */
#endif

using namespace std;


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

double get_value(octave_value val)
{
  NDArray myptr=val.array_value();
  double myval=myptr(0,0);
  return myval;

}

/*--------------------------------------------------------------------------------*/


MPI_Op get_op(octave_value mystr)
{
  charMatrix ch=mystr.char_matrix_value();
  char *myc=get_char_from_arg(ch);

  int len=strlen(myc);
  if (len<=0)
    return MPI_SUM;
  for (int i=0;i<len;i++)
    myc[i]=tolower(myc[i]);

  MPI_Op op=MPI_SUM;
  if ((strcmp(myc,"sum")==0)||(strcmp(myc,"mpi_sum")==0)) 
    op=MPI_SUM;
  
  if ((strcmp(myc,"max")==0)||(strcmp(myc,"mpi_max")==0)) 
    op=MPI_MAX;
  
  if ((strcmp(myc,"min")==0)||(strcmp(myc,"mpi_min")==0))
    op=MPI_MIN;
  
  if ((strcmp(myc,"prod")==0)||(strcmp(myc,"mpi_prod")==0))
    op=MPI_PROD;
  
  free(myc);
  return op;
  
}


/*--------------------------------------------------------------------------------*/

#if 1

DEFUN_DLD (mpi_allreduce, args, nargout, "Reduce .\n")
{
  
  NDArray  val=args(0).array_value();
  dim_vector dm=val.dims();
  NDArray reduced(dm);
  dim_vector dm2=reduced.dims();
  double *valptr=val.fortran_vec();
  double *reducedptr=reduced.fortran_vec();
  
  long numel=1;
  long numel2=1;
  assert(dm.length()==dm2.length());
  for (int i=0;i<dm.length();i++) {
      numel*=dm(i);
      numel2*=dm2(i);
  }
  assert(numel==numel2);
  //printf("Have %ld elements.\n",numel);
  long numel_min,numel_max;
  MPI_Allreduce(&numel,&numel_min,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&numel,&numel_max,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  assert(numel_min==numel);
  assert(numel_max==numel);
  
  

  MPI_Op op=MPI_SUM;
  if (args.length()>1)
    op=get_op(args(1));
  //if (MPI_Allreduce (valptr,reducedptr,numel,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD)) {
  if (MPI_Allreduce (valptr,reducedptr,numel,MPI_DOUBLE,op,MPI_COMM_WORLD)) {
    fprintf(stderr,"Error in MPI_Allreduce.\n");
    return octave_value_list();
  }
  return octave_value(reduced);
  
}

#endif

/*--------------------------------------------------------------------------------*/

DEFUN_DLD (mpi_barrier, args, nargout, "Barrier .\n")
{
  MPI_Barrier(MPI_COMM_WORLD);
  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/
#if 1
DEFUN_DLD (mpi_bcast_string, args, nargout, "Broadcast .\n")
{
  int narg=args.length();

  int from_whom=0;
  if (narg>1)
    from_whom=(int)get_value(args(1));
  


  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD,&rank);

  int len;
  charMatrix ch;
  if (rank==from_whom) {
    ch = args(0).char_matrix_value ();
    len=ch.length();
  }

  MPI_Bcast(&len,1,MPI_INT,from_whom,MPI_COMM_WORLD);
  
  dim_vector dm(1,len);
  charMatrix ch2(dm);
  
  char  *tosend;
  if (rank==from_whom)
    tosend=ch.fortran_vec();
  else
    tosend=ch2.fortran_vec();
  MPI_Bcast(tosend,len,MPI_CHAR,from_whom,MPI_COMM_WORLD);
  
  if (rank==from_whom)
    return octave_value(ch,true);
  else
    return octave_value(ch2,true);
  
  
  return octave_value_list();  //never get here. 
}
#endif
/*--------------------------------------------------------------------------------*/

DEFUN_DLD (copy_string, args, nargout, "Broadcast .\n")
{
  int narg=args.length();

  charMatrix ch = args(0).char_matrix_value ();  
  dim_vector dm=ch.dims();
  charMatrix ch2(dm);
  ch2(ch.length()-1)='a';
  return octave_value(ch2,true);

}
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (mpi_bcast_array, args, nargout, "Broadcast .\n")
{
  int narg=args.length();
  int from_whom=0;
  if (narg>1)
    from_whom=(int)get_value(args(1));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank==from_whom)  { //I am the broadcaster
    NDArray myptr=args(0).array_value();
    dim_vector dm=myptr.dims();
    int ndim=dm.length();
    MPI_Bcast(&ndim,1,MPI_INT,from_whom,MPI_COMM_WORLD);



    //would like to do the following, but no fortran_vec
    //MPI_Bcast(dm.fortran_vec(),ndim,MPI_DOUBLE,from_whom,MPI_COMM_WORLD);
    

    for (int i=0;i<ndim;i++) {
      int cur_dim=dm(i);
      MPI_Bcast(&cur_dim,1,MPI_INT,from_whom,MPI_COMM_WORLD);
    }

    long numel=1;
    for (int i=0;i<ndim;i++)
      numel*=dm(i);
    MPI_Bcast(myptr.fortran_vec(),numel,MPI_DOUBLE,from_whom,MPI_COMM_WORLD);
    return octave_value(myptr);
  }
  else  {
    int ndim;
    MPI_Bcast(&ndim,1,MPI_INT,from_whom,MPI_COMM_WORLD);
    dim_vector dm(ndim);

    //would like to do the following, but no fortran_vec
    // MPI_Bcast(dm.fortran_vec(),ndim,MPI_DOUBLE,from_whom,MPI_COMM_WORLD);
    for (int i=0;i<ndim;i++) {
      int cur_dim;
      MPI_Bcast(&cur_dim,1,MPI_INT,from_whom,MPI_COMM_WORLD);
      dm(i)=cur_dim;
    }
    NDArray arr(dm);
    long numel=1;
    for (int i=0;i<ndim;i++)
      numel*=dm(i);
    MPI_Bcast(arr.fortran_vec(),numel,MPI_DOUBLE,from_whom,MPI_COMM_WORLD);
    
    return octave_value(arr);
  }
  return octave_value_list(); //never get here;
}
