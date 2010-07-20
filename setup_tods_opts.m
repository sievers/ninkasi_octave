function[tods]=setup_tods_opts(tods,opts)

ntod_org=length(tods);
myid=mpi_comm_rank+1;


found_cuts=read_cuts(tods);
if sum(~found_cuts)
  mdisp(['going to cut ' num2str(sum(found_cuts==false)) ' tods on ' num2str(myid)]);
end
tods=tods(found_cuts);

if isfield(opts,'endcut')
  crap=cut_tod_ends(tods,opts.endcut);
else
  disp('skipping chopping of tod ends.');
end


ntod=length(tods);
for j=1:ntod,
  if opts.do_calib,
    get_tod_calib_factors(tods(j));  %this also cuts uncalibrated detectors.
  end
  purge_cut_detectors(tods(j));
end


if isfield(opts,'cut_magic_carpets')
  if opts.cut_magic_carpets,
    mdisp('cutting magic carpets.');
    cut_magic_carpets(tods);
  end
end

get_tod_uncut_regions(tods);

for j=1:length(tods),
  cut_mostly_cut_detectors(tods(j));
  purge_cut_detectors(tods(j));
end

tod_names=get_tod_names(tods);
if isfield(opts,'mindets')
  tods=cull_empty_tods(tods,tod_names,opts.mindets);
end
if isfield(opts,'minsamps')
  tods=cull_short_tods(tods,tod_names,opts.minsamps);
end

if isfield(opts,'distthresh')
  total_tods=mpi_allreduce(length(tods));
  tods=cull_distant_tods(tods,opts);
  new_total=mpi_allreduce(length(tods));
  mdisp(['went from ' num2str(total_tods) ' to ' num2str(new_total) ' tods from distance cut.']);
else
  error('not cutting');
end



if length(tods)~=ntod_org,
  mdisp(['cutting ' num2str(ntod_org-length(tods)) ' tods on ' num2str(myid)]);
end

