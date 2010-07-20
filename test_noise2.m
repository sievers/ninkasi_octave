
more off

system('hostname');
fftw_init_threaded

addpath /home/sievers/matlab
mpi_init;
myid=mpi_comm_rank+1;
nproc=mpi_comm_size;

disp(['greetings from ' num2str(myid)])



do_calib=true;
highpass=0.5;



set='2';


format short g

tod_names={
%'20080914/1221369887.1221369955.ar1'
%'20080914/1221366269.1221366353.ar1'
'20080912/1221196003.1221230597.ar1'
%'20081104/1225797346.1234651384.ar1'
%'20080918/1221705782.1221705825.ar1'
%'20080914/1221382603.1221382663.ar1'
};

mdisp(['starting with ' num2str(length(tod_names)) ' TODs.']);

for j=1:length(tod_names),
  tt=tod_names{j};if tt(end)=='/', tod_names(j)={tt(1:end-1)};end;
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

do_noise=true;nbad=10;
bands=[0.01*[-1 2.^(1:12)] 1000];
rots=zeros(size(bands));
noise_types=zeros(size(bands));
rots=rots(1:end-1);
noise_types=noise_types(1:end-1);


found_cuts=read_cuts(tods);
if sum(found_cuts)
  disp(['going to cut ' num2str(sum(found_cuts==false)) ' tods on ' num2str(myid)]);
end
tods=tods(found_cuts);
tod_names=tod_names(found_cuts);


%crap=cut_tod_ends(tods,20);

ntod=length(tods);

for j=1:ntod,
  if do_calib,
    get_tod_calib_factors(tods(j));  %this also cuts uncalibrated detectors.
  end

  purge_cut_detectors(tods(j));
end

get_tod_uncut_regions(tods);


crap=cut_tod_ends(tods,20);

cut_mostly_cut_detectors(tods(j));	
purge_cut_detectors(tods(j));

[tods,tod_names]=cull_empty_tods(tods,tod_names,400);

if (length(tods)==0) 
   error(['no tods on process ' num2str(myid)]);
end



ntod=length(tods);

map=allocate_ninkasi_skymap(pixsize,lims(1)-1e-2,lims(2)+1e-2, lims(3)-1e-2,lims(4)+1e-2);
set_skymap_cea_simple_c(map);
mdisp('set cea');

mapset.skymap.mapptr=map;
mapset.skymap.map=skymap2octave(mapset.skymap.mapptr);


for j=1:ntod,
  mdisp(['crunching  tod ' num2str(j)]);
  [rows,cols]=get_tod_rowcol(tods(j));
  mapset.corrnoise(j).map=zeros(get_tod_ndata(tods(j)),3);
  mapset.corrnoise(j).vecs=[ones(1,how_many_dets_are_kept(tods(j)));rows'-mean(rows);cols'-mean(cols)];
end

mdisp('setup mapset.');


fac=0;
outroot=['temp_map_season_strip_set' set '_small_cuts3_filt_' num2str(highpass) '_prior_x' num2str(fac) '_find_' num2str(nbad) '_modes_'];
if (do_noise),
   outroot=[outroot 'noise_'];
end

wt=make_map_copy(map);

make_weightmap_octave(tods,wt);
mdisp('made weightmap');

weight=skymap2octave(wt);
weight=mpi_allreduce(weight);
octave2skymap(weight,wt);
if (myid==1)
  write_map(wt,[outroot 'weights']);
end

destroy_map(wt);
precon.skymap.map=weight;
precon.corrnoise=mapset.corrnoise;

mapset2=clear_mapset(mapset,true);


[b,meds]=create_initial_mapset_octave(tods,mapset,'', 'scale_factor',1.0,'find_modes',true,'nbadmode',nbad,'detrend',true,'highpass',highpass,'debutter',true,'do_noise',do_noise,'do_calib',do_calib,'bands',bands,'rots',rots,'noise_types',noise_types);


return



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
