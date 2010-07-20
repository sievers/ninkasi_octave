more off
addpath /cita/d/raid-sievers/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;


%[tods_org,lims]=read_all_tod_headers(); pixsize=30/60/60*pi/180;
[tods,lims]=read_all_tod_headers(1,2); pixsize=30/60/60*pi/180;



ntod=length(tods);
for j=1:ntod,
  purge_cut_detectors(tods(j));
end




allocate_tod_storage(tods(1));

tic;read_tod_data(tods(1));toc


tic;corrs=get_data_correlation(tods(1));toc;

tic;a=get_band_fft(tods(1),0.002,4);toc;
