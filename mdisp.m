function[]=mdisp(mystr)
myid=mpi_comm_rank();
if (myid==0)
  disp(mystr);
end
