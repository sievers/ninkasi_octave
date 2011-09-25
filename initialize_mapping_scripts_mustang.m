%addpath /home/sievers/matlab
format short g
more off
crash_dumps_octave_core(false);

%fftw_init_threaded

mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

if nproc>1
  ignore_function_time_stamp ("all")
end

test_omp_c(1);  %set to single-threaded
