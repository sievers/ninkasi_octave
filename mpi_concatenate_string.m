function[big_str]=mpi_concatenate_string(str)



myid=mpi_comm_rank;
nproc=mpi_comm_size;
nn=mpi_allreduce(1);assert(nn==nproc);
big_str='';
for j=0:nproc-1,
  if (myid==j)
    mystr=mpi_bcast_string(str,j);
    big_str=[big_str str];
  else
    mystr=mpi_bcast_string('',j);
    big_str=[big_str mystr];
  end

end

