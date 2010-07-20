function[seed_in]=mpi_setup_seeds(seed_in)
if ~exist('seed_in')
  seed_in=round(1000*rand);
end

seed_in=mpi_bcast_array(seed_in,0);
nproc=mpi_comm_size;
rand ('state', seed_in)

myid=mpi_comm_rank;
myseed=round(1e12*rand(1))+myid;

%printf('starting seed is %d on %d\n',myseed,myid)
rand('state',myseed);

