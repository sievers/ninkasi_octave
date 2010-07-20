mpi_init;

n=mpi_comm_size;
 myid=mpi_comm_rank;
disp(['I am ' num2str(myid) ' of ' num2str(n)]);

x=randn(3);

tot=mpi_allreduce(x);
disp('private x is :')
disp(x);

disp('sum is: ')
disp(tot);

mpi_finalize;