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
void *get_pointer(octave_value val)
{
  int64NDArray myptr=val.array_value();
  long myptr2=myptr(0,0);
  return (void *)myptr2;

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

#if 1

DEFUN_DLD (mpi_reduce, args, nargout, "Reduce.  Args are (array, destination, operator, communicator) (\n")
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
  
  int mpi_target=0;
  //find the target - octave is unit-offset, so shift down by 1
  if (args.length()>1)
    mpi_target=(int)get_value(args(1))-1;
  
  
  
  MPI_Op op=MPI_SUM;
  if (args.length()>2)
    op=get_op(args(2));
  MPI_Comm comm=MPI_COMM_WORLD;
  if (args.length()>3) {
    MPI_Comm *myptr=(MPI_Comm *)get_pointer(args(3));
    comm=*myptr;
  }
  
  

  if (MPI_Reduce (valptr,reducedptr,numel,MPI_DOUBLE,op,mpi_target,comm)) {
    fprintf(stderr,"Error in MPI_Reduce.\n");
    return octave_value_list();
  }
  
  int myrank=-1;
  MPI_Comm_rank(comm,&myrank);
  if (myrank==mpi_target)
    return octave_value(reduced);
  else {
    Matrix my_empty;
    //return octave_value_list();
    return octave_value(my_empty);
  }
  
}

#endif






/*--------------------------------------------------------------------------------*/
//int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,	     int tag, MPI_Comm comm)
DEFUN_DLD(mpi_send,args,nargout,"MPI Send.  args are data, destination, <tag>, <communicator>\n")
{
  
  if (args.length()<2) {
    printf("Need at least 2 arguments in mpi_send.  only got %d\n",args.length());
    return octave_value_list();
  }
  int dest=(int)get_value(args(1))-1;
  int tag=0;
  if (args.length()>2)
    tag=get_value(args(2));

  MPI_Comm comm=MPI_COMM_WORLD;
  if (args.length()>3) {
    MPI_Comm *myptr=(MPI_Comm *)get_pointer(args(3));
    comm=*myptr;
  }
   
  NDArray  val=args(0).array_value();
  dim_vector dm=val.dims();
  int ndim=dm.length();

  int dims[ndim];
  int nelem=1;
  for (int i=0;i<ndim;i++) {
    nelem*=dm(i);
    dims[i]=dm(i);
  }
  
  MPI_Send(&ndim,1,MPI_INT,dest,tag,comm);
  MPI_Send(dims,ndim,MPI_INT,dest,tag,comm);

  MPI_Send(val.fortran_vec(),nelem,MPI_DOUBLE,dest,tag,comm);

  return octave_value_list();
}

/*--------------------------------------------------------------------------------*/

//int MPI_Recv(void *buf, int count, MPI_Datatype datatype,	     int source, int tag, MPI_Comm comm, MPI_Status *status)


DEFUN_DLD(mpi_recv,args,nargout,"MPI Receive.  args are source, <tag>, <communicator>\n")
{
  
  if (args.length()<1) {
    printf("Need at least 1 argument in mpi_recv (source).  only got %d\n",args.length());
    return octave_value_list();
  }
  MPI_Status mystatus;
  int source=(int)get_value(args(0))-1;
  int tag=0;
  if (args.length()>1)
    tag=get_value(args(1));
  
  MPI_Comm comm=MPI_COMM_WORLD;
  if (args.length()>2) {
    MPI_Comm *myptr=(MPI_Comm *)get_pointer(args(2));
    comm=*myptr;
  }
   
  int ndim;
  MPI_Recv(&ndim,1,MPI_INT,source,tag,comm,&mystatus);
  int dims[ndim];
  MPI_Recv(dims,ndim,MPI_INT,source,tag,comm,&mystatus);
  //according to the header file, this is how to get n-dimensional dim_vectors assembled.
  dim_vector dm;
  if (ndim==1) {
    dim_vector dm_tmp(dims[0]);
    dm=dm_tmp;
  }

  if (ndim==2) {
    dim_vector dm_tmp(dims[0],dims[1]);
    dm=dm_tmp;
  }

  if (ndim==3) {
    dim_vector dm_tmp(dims[0],dims[1],dims[2]);
    dm=dm_tmp;
  }

  if (ndim==4) {
    dim_vector dm_tmp(dims[0],dims[1],dims[2],dims[3]);
    dm=dm_tmp;
  }
  int nelem=1;
  for (int i=0;i<ndim;i++) {
    dm(i)=dims[i];
    nelem*=dm(i);
  }
  NDArray val(dm);
  dim_vector dm2=val.dims();

  double *myptr=val.fortran_vec();
  //memset(myptr,0,sizeof(double)*nelem);

  MPI_Recv(val.fortran_vec(),nelem,MPI_DOUBLE,source,tag,comm,&mystatus);
  return octave_value(val);
}


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
/*--------------------------------------------------------------------------------*/
DEFUN_DLD (mpi_allgather_c, args, nargout, "Allgather(v).  Note that this will turn everything into a column vector.\n")
{
  Matrix mat=args(0).matrix_value();
  dim_vector dm=mat.dims();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  int nelem=dm(0)*dm(1);
  int *all_nelem=(int *)malloc(sizeof(int)*nproc);
  int *displs=(int *)malloc(sizeof(int)*nproc);
  displs[0]=0;
  memset(all_nelem,0,sizeof(int)*nproc);
  MPI_Allgather(&nelem,1,MPI_INT,all_nelem,1,MPI_INT,MPI_COMM_WORLD);
  int nelem_tot=0;
  for (int i=0;i<nproc;i++) {
    nelem_tot+=all_nelem[i];
    if (i>0)
      displs[i]=displs[i-1]+all_nelem[i-1];
  }

  Matrix allmat(nelem_tot,1);
  dim_vector dd=allmat.dims();
  double *allptr=allmat.fortran_vec();
  double *ptr=mat.fortran_vec();

  MPI_Allgatherv(ptr,nelem,MPI_DOUBLE,allptr,all_nelem,displs,MPI_DOUBLE,MPI_COMM_WORLD);
		 
  free(all_nelem);
  free(displs);
  return octave_value(allmat);
}

