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

%tod_names=read_text_file('/home/sievers/act//tod_lists/bright_source.txt');
tod_names={'/scr/queequeg1/colossus/season2/merlin/20080906/1220698325.1235547164.ar1'};

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





[rr,cc]=get_tod_rowcol(tods);
ind=(rr==6)&(cc==6);
assert(sum(ind)==1);
[a,myind]=max(ind);



read_tod_data(tods);
dat=get_tod_data(tods);
dat_raw=dat(:,myind);

debutterworth_c(tods);
dat=get_tod_data(tods);
dat_debutter=dat(:,myind);

calibrate_data(tods);
dat=get_tod_data(tods);
dat_debutter_calib=dat(:,myind);

fid=fopen(['for_toby_det66_decimate_' num2str(decimate) '.dat'],'w');
fwrite(fid,length(dat_raw),'int');
fwrite(fid,dat_raw,'double');
fwrite(fid,dat_debutter,'double');
fwrite(fid,dat_debutter_calib,'double');
fclose(fid);



fid=fopen(['for_toby_det66_decimate_' num2str(decimate) '_raw.txt'],'w');
fprintf(fid,'%16.7e\n',dat_raw);
fclose(fid);

fid=fopen(['for_toby_det66_decimate_' num2str(decimate) '_debutter.txt'],'w');
fprintf(fid,'%16.7e\n',dat_debutter);
fclose(fid);

fid=fopen(['for_toby_det66_decimate_' num2str(decimate) '_debutter_calib.txt'],'w');
fprintf(fid,'%16.7e\n',dat_debutter_calib);
fclose(fid);



return




calibrate_data(tods);

dat=get_tod_data(tods);
