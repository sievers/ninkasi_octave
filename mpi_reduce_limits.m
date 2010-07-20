function[lims]=mpi_reduce_limits(lims)
lims(1)=mpi_allreduce(lims(1),'min');
lims(2)=mpi_allreduce(lims(2),'max');
lims(3)=mpi_allreduce(lims(3),'min');
lims(4)=mpi_allreduce(lims(4),'max');

