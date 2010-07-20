%mpiexec --hostfile ~/f120n063.txt octave --eval "myset='0';" test_octave_mapping_season_strip_cuts3_4way_nofilt_025.m | tee split_4way_set0.out

aa=now;

more off
crash_dumps_octave_core(false);
fftw_init_threaded

addpath /home/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

%disp(['greetings from ' num2str(myid)])


do_calib=true;
highpass=0.0;

do_noise=true;nbad=10;

bands=[-1 0.5 0.7 1 1.4 2 2.8 4 5.6 8 11 16 22 32 45 64 1000];
noise_scale_facs=[0 0.25 0.5 ];
rots=zeros(size(bands));noise_types=zeros(size(bands));
rots=rots(1:end-1);
noise_types=noise_types(1:end-1);







format short g

tod_names=read_text_file('/home/sievers/act//tod_lists/bright_source.txt');


mdisp(['starting with ' num2str(length(tod_names)) ' TODs.']);

for j=1:length(tod_names),
  tt=tod_names{j};if tt(end)=='/', tod_names(j)={tt(1:end-1)};end;
end


if nproc==1,
   tod_names=tod_names(1);
   highpass=1.0;
end




tod_names=tod_names(myid:nproc:end);
tod_names=guess_tod_name(tod_names);



decimate=1;
[tods,lims]=read_all_tod_headers(tod_names,decimate); pixsize=30/60/60*pi/180;
mdisp(['lims are ' num2str(lims)]);
lims(1)=mpi_allreduce(lims(1),'min');
lims(2)=mpi_allreduce(lims(2),'max');
lims(3)=mpi_allreduce(lims(3),'min');
lims(4)=mpi_allreduce(lims(4),'max');
mdisp(['lims are now ' num2str(lims)]);



found_cuts=read_cuts(tods);
if sum(found_cuts)
  disp(['going to cut ' num2str(sum(found_cuts==false)) ' tods on ' num2str(myid)]);
end
tods=tods(found_cuts);
tod_names=tod_names(found_cuts);
crap=cut_tod_ends(tods,20);

ntod=length(tods);

for j=1:ntod,
  if do_calib,
    get_tod_calib_factors(tods(j));  %this also cuts uncalibrated detectors.
  end
  purge_cut_detectors(tods(j));
end

get_tod_uncut_regions(tods);


for j=1:length(tods),
  cut_mostly_cut_detectors(tods(j));
  purge_cut_detectors(tods(j));
end




[tods,tod_names]=cull_empty_tods(tods,tod_names,400);
if (length(tods)==0) 
   error(['no tods on process ' num2str(myid)]);
end

ntod=length(tods);


map=allocate_ninkasi_skymap(pixsize,lims(1)-1e-2,lims(2)+1e-2, lims(3)-1e-2,lims(4)+1e-2);
set_skymap_cea_simple_c(map);



mapset.skymap.mapptr=map;
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);



wt=make_map_copy(map);

make_weightmap_octave(tods,wt);
mdisp('made weightmap');





weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
octave2skymap(weight,wt);

destroy_map(wt);




read_tod_data(tods)
calibrate_data(tods);

dat=get_tod_data(tods);
meds=omp_median(dat,1);
%tic;for j=1:length(meds),dat(:,j)=dat(:,j)-meds(j);end;toc;
tic;crud=repmat(meds,[size(dat,1) 1]);dat=dat-crud;clear crud;toc;
cm=omp_median(dat,2);
facs=(cm'*dat)/(cm'*cm);

dat2=dat-cm*facs;

tic;crud=dat2';mycorr=crud*dat2;toc

mycorr=0.5*(mycorr+mycorr');

[vv,ee]=eig(mycorr);ee=diag(ee);ee=sqrt(ee);ind=ee>3*median(ee);
vecs=vv(:,ind);
mymodes=dat2*vecs;
fit_params=mymodes'*dat2;
dat_clean=dat2-mymodes*fit_params;

fid=fopen('dat_clean.out','w');fwrite(fid,size(dat_clean),'int');fwrite(fid,dat_clean,'double');fclose(fid);

bb=now;
disp(86400*(bb-aa))


clear_map(map);
tod2map(tods,map);
mm=skymap2octave(map);
weight_inv=1./weight;weight_inv(weight==0)=0;  
mm=mm.*weight_inv;
octave2skymap(mm,map);




return




precon.skymap.map=weight;
precon.corrnoise=mapset.corrnoise;

mapset2=clear_mapset(mapset,true);


[b,meds]=create_initial_mapset_octave(tods,mapset,'', 'scale_factor',-1.0,'find_modes',true,'nbadmode',nbad,'detrend',true,'highpass',highpass,'debutter',true,'do_noise',do_noise,'do_calib',do_calib,'bands',bands,'rots',rots,'noise_types',noise_types,'noise_scale_facs',noise_scale_facs);

weight_inv=weight;weight_inv(weight>0)=1./weight(weight>0);
x=clone_mapset(mapset);
x.corrnoise=b.corrnoise;
for j=1:length(x.corrnoise),
  x.corrnoise(j).vecs=b.corrnoise(j).vecs;
  x.corrnoise(j).map(:,1)=meds{j}/mean(x.corrnoise(j).vecs(1,:));
end
precon.corrnoise=x.corrnoise;

x=fit_mapset_tictoc(tods,x,b,weight_inv);
if (myid==1)  
  write_map(x.skymap,[outroot 'starting']);
end

  
if (0)

  mdisp('getting priors');
  for j=1:length(tods),
    x.corrnoise(j).prior=1./get_smoothed_ps(x.corrnoise(j).map).^2*fac^2;    
    b.corrnoise(j).prior=1./get_smoothed_ps(x.corrnoise(j).map).^2*fac^2;    
  end
  precon.corrnoise=x.corrnoise;
  mdisp('got em');
  
else
  outroot=[outroot 'noprior_'];
end

[x,best]=run_pcg_corrnoise_precon_octave(tods,x,b,precon,@apply_precon,@apply_prior_to_mapset,'maxiter',1000,'tol',1e-10,'save_interval',10,'save_tag',outroot,'do_noise',do_noise);

mpi_finalize;return
