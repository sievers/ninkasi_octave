function[new_mapset,tod_times,node_time]=mapset2mapset_corrnoise_octave(tods,mapset,varargin)

do_noise=get_keyval_default('do_noise',false,varargin{:});

new_mapset=clear_mapset(mapset,'true');

if isfield(mapset,'skymap')
  octave2skymap(mapset.skymap);
end


%tod_times=zeros(size(tods));
tod_times=zeros(length(tods),4);


node_start=now;
for j=1:length(tods),
  t1=now;
  mtic
  mytod=tods(j);
  allocate_tod_storage(mytod);
  assign_tod_value(mytod,0.0);

  %mtoc
  t1a=now;
  mapset2tod_octave(mapset,mytod,j);
  t1b=now;

  mtoc


  if (0)
    crap=get_tod_data(mytod);
    disp(sprintf('tod %d data %15.6e',j,sum(sum(crap.^2))))
    clear crap
  end

  window_data(mytod);


  %mtoc
  %do noise in here
  t2a=now;
  %if (do_noise)
  %apply_banded_noise_model_c(mytod);
  %end	
  apply_tod_noise_model_c(mytod);
  t2b=now;
  mtoc
  t3a=now;
  new_mapset=tod2mapset_octave(new_mapset,mytod,j);
  t3b=now;
  %mtoc
  if (length(tods>1))
    free_tod_storage(mytod);
  else
    assign_tod_value(mytod,0.0);
  end
  mtoc
  t2=now;
  tod_times(j,1)=86400*(t2-t1);
  tod_times(j,2)=86400*(t1b-t1a);
  tod_times(j,3)=86400*(t2b-t2a);
  tod_times(j,4)=86400*(t3b-t3a);
  
end
node_time=86400*(now-node_start);


%need to make sure we don't memory leak...
if isfield(mapset,'skymap')
  new_mapset.skymap.map=skymap2octave(new_mapset.skymap.mapptr);
  new_mapset.skymap.map=mpi_allreduce(new_mapset.skymap.map);
  destroy_map(new_mapset.skymap.mapptr);
  new_mapset.skymap.mapptr=mapset.skymap.mapptr;
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