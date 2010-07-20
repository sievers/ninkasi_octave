function[value]=odisp(mystr,pauselen)
%do an ordered print to screen from mpi processes.  Note - could be very slow.
%pauselen tells you how long to wait between processes.

myid=mpi_comm_rank;
if ~exist('pauselen'),
  pauselen=0.05;
end
np=mpi_allreduce(1);
pause(pauselen*myid)
disp(mystr)
