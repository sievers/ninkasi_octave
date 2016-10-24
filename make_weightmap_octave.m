function[map]=make_weightmap_octave(tods,map,do_window,varargin)

do_new_pointing=get_keyval_default('do_new_pointing',false,varargin{:})
do_reduce=get_keyval_default('mpi_reduce',false,varargin{:});
free_2gamma=get_keyval_default('free_2gamma',true,varargin{:});
do_noise=get_keyval_default('do_noise',false,varargin{:});  %added ability to do noise instead of hitcounts.


if isstruct(map)
  if isfield(map,'partition')
    clear_map(map.mapptr);
    make_weightmap_octave(tods,map.mapptr,do_window,'mpi_reduce',false,varargin{:});
    if do_reduce
      map=mpi_reduce_map(map);
    end
    return;
  end
end



clear_map(map);

if isempty('do_window')
  clear do_window;
end


if ~exist('do_window')
  do_window=true;
end

for j=1:length(tods),    
  mytod=tods(j);
  allocate_tod_storage(mytod);
  if (do_noise)
    mynoise=get_detector_average_noise_banded_projvecs(mytod);
    nn=get_tod_ndata(mytod);
    wt_per_samp=1./mynoise.^2;
    dat=repmat(wt_per_samp',[nn 1]);
    push_tod_data(dat,mytod);
    clear dat;
    array_noise=1/sqrt(sum(wt_per_samp))*sqrt(get_tod_dt(mytod));
    disp(['TOD ' get_tod_name(mytod) ' has array noise ' num2str(array_noise)]);

  else
    assign_tod_value(mytod,1.0);
  end
  if do_window,
    window_data(mytod);
    window_data(mytod);
  end

  
  if (do_new_pointing)
    mdisp('doing pointing')
    precalc_actpol_pointing_exact(mytod);
    convert_saved_pointing_to_pixellization(mytod,map);
    free_tod_pointing_saved(mytod,free_2gamma);
  end
  
  mdisp('assigned value');
  if (is_map_polarized(map))
    tod2polmap(mytod,map);
  else
    tod2map(mytod,map);
  end
  mdisp('projected');
  
  if (do_new_pointing)
    %free_tod_pointing_saved(mytod);
    free_saved_pixellization(mytod);
    
  end


  free_tod_storage(mytod);    
  %disp('freed');
end


if (do_reduce)
  %disp(['at barrier ' num2str(mpi_comm_rank+1)]);
  mpi_barrier;
  mdisp('reducing')
  weight=skymap2octave(map);
  %whos weight
  %scinet mpi_allreduce started to act up for large matrices, so split up polarization ones
  if length(size(weight))==3,
    for j=1:size(weight,1),
      weight(j,:,:)=mpi_allreduce(weight(j,:,:));
    end
  else
    weight=mpi_allreduce(weight);
  end
  octave2skymap(weight,map);
  clear weight;
end

  
if (do_noise)
  weight=skymap2octave(map);
  if numel(size(weight))==3,
    asdf=weight(1,:,:);
    asdf(asdf>0)=1./sqrt(asdf(asdf>0));
    weight(1,:,:)=asdf;
  else
    weight(weight>0)=1./sqrt(weight(weight>0));
  end
  octave2skymap(weight,map);
end