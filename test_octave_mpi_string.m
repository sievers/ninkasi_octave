mpi_init;
myid=mpi_comm_rank;

if ( myid == 0 )
  ss='abcdef';
  mpi_bcast_string(ss);
else
  ss=mpi_bcast_string;
end

disp(['string is ' ss ' on ' num2str(myid)]);


mpi_finalize;