addpath /home/sievers/matlab

more off
crash_dumps_octave_core(false);



fftw_init_threaded

mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

if nproc>1
  ignore_function_time_stamp ("all")
end
