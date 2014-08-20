function[new_mapset,tod_times,node_time,mpi_time]=mapset2mapset_corrnoise_octave(tods,mapset,varargin)

if (1)
  if numel(varargin)>1,
    clear myopts;
    for j=1:2:numel(varargin)
      eval(['myopts.' varargin{j} ' = varargin{j+1};']);
    end
  else
    myopts=varargin{1};
  end

  do_noise=get_struct_mem(myopts,'do_noise',false);
  window_symmetric=get_struct_mem(myopts,'window_symmetric',false);
  skip_window=get_struct_mem(myopts,'skip_window',false);
  barriers=get_struct_mem(myopts,'do_barriers',false);
  skip_mpi=get_struct_mem(myopts,'skip_mpi',false);
  check_empty=get_struct_mem(myopts,'check_empty',false);
  new_mapptr=get_struct_mem(myopts,'new_mapptr',[]);
  do_actpol_pointing=get_struct_mem(myopts,'do_actpol_pointing',false);
else
  
  do_noise=get_keyval_default('do_noise',false,varargin{:});
  window_symmetric=get_keyval_default('window_symmetric',false,varargin{:});
  skip_window=get_keyval_default('skip_window',false,varargin{:});  %don't window at all (e.g. if there's no noise)
  barriers=get_keyval_default('do_barriers',false,varargin{:});
  skip_mpi=get_keyval_default('skip_mpi',false,varargin{:});  %for some things, we can do reduces later, particularly useful when doing source fits.
  check_empty=get_keyval_default('check_empty',false,varargin{:});  %in some cases, it may be worth seeing if the TOD is empty before applying noise
end

new_mapset=clear_mapset(mapset,'true');
if ~isempty(new_mapptr)
  new_mapset.skymap.mapptr=make_map_copy(new_mapptr);
  clear_map(new_mapset.skymap.mapptr);
  new_mapset.skymap.map=skymap2octave(new_mapset.skymap.mapptr);
end

if isfield(mapset,'skymap')
  octave2skymap(mapset.skymap);
end
if isfield(mapset,'ground')
  octave2skymap(mapset.ground.ground,mapset.ground.groundptr);
end



%tod_times=zeros(size(tods));
tod_times=zeros(length(tods),4);

if barriers,
  min_ntod=mpi_allreduce(numel(tods),'min');
end


node_start=now;
for j=1:length(tods),
  t1=now;
  mtic
  mytod=tods(j);
  
  if (do_actpol_pointing)
    precalc_actpol_pointing_exact(mytod);
    if isfield(mapset,'skymap')
      if isfield(mapset.skymap,'mapptr')
        convert_saved_pointing_to_pixellization(mytod,mapset.skymap.mapptr)
      end
    end
    free_tod_pointing_saved(mytod);
  end
  
  allocate_tod_storage(mytod);
  assign_tod_value(mytod,0.0);

  %mtoc
  t1a=now;
  mapset2tod_octave(mapset,mytod,j);
  t1b=now;
  

  mtoc
  if (barriers)
    if j<=min_ntod
      mpi_barrier;
      mdisp(['completed mapset2tod on ' num2str(j)]);
    end
  end



  if (0)
    crap=get_tod_data(mytod);
    disp(sprintf('tod %d data %15.6e',j,sum(sum(crap.^2))))
    clear crap
  end

  if (~skip_window)
    window_data(mytod);
  end



  %mtoc
  %do noise in here
  t2a=now;
  %if (do_noise)
  %apply_banded_noise_model_c(mytod);
  %end	
  
  skip_noise=false;
  if check_empty,
    skip_noise=(sum_tod_data(mytod)==0);
  end
  if ~skip_noise,
    %apply_tod_noise_model_c(mytod);
    apply_tod_noise_model(mytod);
  end


  t2b=now;
  mtoc
  if (barriers)
    if j<=min_ntod
      mpi_barrier;
      mdisp(['completed noise on ' num2str(j)]);
    end
  end


  if (window_symmetric)
    window_data(mytod);
  end
  t3a=now;
  new_mapset=tod2mapset_octave(new_mapset,mytod,j);
  t3b=now;
  if (barriers)
    if j<=min_ntod
      mpi_barrier;
      mdisp(['completed tod2mapset on ' num2str(j)]);
    end
  end

  %mtoc
  if (length(tods>1))
    free_tod_storage(mytod);
  else
    assign_tod_value(mytod,0.0);
  end
  
  if (do_actpol_pointing)
    %free_tod_pointing_saved(mytod);
    free_saved_pixellization(mytod);
  end
  
  
  mtoc
  t2=now;
  tod_times(j,1)=86400*(t2-t1);
  tod_times(j,2)=86400*(t1b-t1a);
  tod_times(j,3)=86400*(t2b-t2a);
  tod_times(j,4)=86400*(t3b-t3a);
  
end
node_time=86400*(now-node_start);
%mdisp(tod_times)

mpi_time=0;
%need to make sure we don't memory leak...
if isfield(mapset,'skymap')
  new_mapset.skymap.map=skymap2octave(new_mapset.skymap.mapptr);
  aa=now;

  if (~skip_mpi)
    global weight_template
    if ~isempty(weight_template)
      fwee=mpi_allreduce(new_mapset.skymap.map(weight_template));
      new_mapset.skymap.map(weight_template)=fwee;
    else
      new_mapset.skymap.map=mpi_allreduce(new_mapset.skymap.map);
    end
  end
  mpi_time=86400*(now-aa);

  %if we are projecting into a new map pointer, don't free stuff
  if isempty(new_mapptr)
    destroy_map(new_mapset.skymap.mapptr);
    new_mapset.skymap.mapptr=mapset.skymap.mapptr;
  end
end

if isfield(mapset,'ground')
  %mdisp(['current ground tot is ' num2str(sum(sum(sum(abs(skymap2octave(new_mapset.ground.groundptr))))))]);
  if ~skip_mpi
    mpi_reduce_map(new_mapset.ground.groundptr);
  end
  new_mapset.ground.ground=skymap2octave(new_mapset.ground.groundptr);
  if isempty(new_mapptr)
    destroy_map(new_mapset.ground.groundptr);
    new_mapset.ground.groundptr=mapset.ground.groundptr;
  end

end


if isfield(mapset,'srccat')
  if ~skip_mpi
    if iscell(new_mapset.srccat)
      for ss=1:numel(new_mapset.srccat)
        new_mapset.srccat{ss}.amps=mpi_allreduce(new_mapset.srccat{ss}.amps);
      end
    else
      new_mapset.srccat.amps=mpi_allreduce(new_mapset.srccat.amps);
    end
  end
end





function[value]=mtic()
if (1)
  tic;
end

function[value]=mtoc()
if (0)
  if (mpi_comm_rank)==0
    toc;
  end
end
%toc;