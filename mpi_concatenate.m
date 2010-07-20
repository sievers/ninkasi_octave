function[big_mat]=mpi_concatenate(mat,which_dim)

if ~exist('which_dim')
  which_dim=1;
end


myid=mpi_comm_rank;
nproc=mpi_comm_size;
big_mat=[];
for j=0:nproc-1,
  if (myid==j)
    mymat=mpi_bcast_array(mat,j);
  else
    mymat=mpi_bcast_array('',j);
  end
  if (which_dim==1)
    big_mat=[big_mat;mymat];
  else
    big_mat=[big_mat mymat];
  end
end
