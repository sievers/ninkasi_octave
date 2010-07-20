fftw_init_threaded

more off
addpath /home/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;

do_calib=true;

nproc=mpi_comm_size;
format short g
dr='/home/sievers/act/data/';
tod_names={
'1221194194.1235589444.ar1'
'1221993788.1221993903.ar1'
'1221390748.1221390761.ar1'
'1222252186.1222252283.ar1'
'1223630848.1223630964.ar1'
'1223700209.1234326370.ar1'
'1223717251.1234327922.ar1'
'1223975543.1223975689.ar1'
'1224130828.1224130940.ar1'
'1224217241.1224217412.ar1'
'1224475437.1224475573.ar1'
'1224578536.1224578653.ar1'
%%%%'1224923283.1224923419.ar1'
'1225009631.1225009755.ar1'
'1225078521.1225078640.ar1'
'1225267927.1225268037.ar1'
'1225354329.1225354462.ar1'
'1225699022.1234639633.ar1'
'1225768004.1234649611.ar1'
'1226112714.1226112873.ar1'
'1226285106.1226357005.ar1'
'1226302013.1226348217.ar1'
'1226560306.1226619322.ar1'
'1226646707.1226646820.ar1'
'1226991401.1226991516.ar1'
'1227491393.1227491515.ar1'
'1227939086.1227939191.ar1'
'1228025486.1228025589.ar1'
'1228283781.1228283887.ar1'
'1228783817.1228783833.ar1'
};


for j=1:length(tod_names),
  tod_names{j}=[dr tod_names{j}];
end
if (nproc==1)
  tod_names=tod_names(1);
else
  tod_names=tod_names(myid:nproc:end);
end

tod_names=guess_tod_name(tod_names);

%tod_names=tod_names(1:2:end);


decimate=1;
[tods,lims]=read_all_tod_headers(tod_names,decimate); pixsize=30/60/60*pi/180;
mdisp(['lims are ' num2str(lims)]);
lims(1)=mpi_allreduce(lims(1),'min');
lims(2)=mpi_allreduce(lims(2),'max');
lims(3)=mpi_allreduce(lims(3),'min');
lims(4)=mpi_allreduce(lims(4),'max');
mdisp(['lims are now ' num2str(lims)]);

ntod=length(tods);

do_noise=true;nbad=10;
do_rots=[0 0];

read_cuts(tods);
crap=cut_tod_ends(tods,20);

for j=1:ntod,
  if do_calib,
    get_tod_calib_factors(tods(j));  %this also cuts uncalibrated detectors.
  end
  purge_cut_detectors(tods(j));
end



get_tod_uncut_regions(tods);


allocate_tod_storage(tods);read_tod_data(tods);
      calibrate_data(tods);
gapfill_data_c(tods);
detrend_data_c(tods);debutterworth_c(tods);
window_data(tods);

dat=get_tod_data(tods);
bands=[0.01*[-1 2.^(1:12)] 1000];



mytod=tods;
allocate_tod_noise_bands_c(mytod,bands);

%rots=[0 0];
%noise_types=[1 0];

rots=zeros(size(bands));rots=rots(1:end-1);
noise_types=zeros(size(bands));noise_types=noise_types(1:end-1);

get_simple_banded_noise_model_c(mytod,rots,noise_types);
tic;apply_banded_noise_model_c(tods);toc;
dat2=get_tod_data(tods);

%ii=500;a=abs(fft(dat(:,ii)));b=abs(fft(dat2(:,ii)));loglog(a(1:end/2).*b(1:end/2))
